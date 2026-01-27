import time
import os
import anndata as ad
from datetime import datetime, timedelta

import requests

from task_manager import Session, ProteinTask, get_s3_client, send_notification
from dgat_utils.predict_util import web_predict
from dgat_utils.utils.Preprocessing import preprocess_ST

from dgat_utils.downstream import _plot_leiden_clustering, _plot_spatial_expr, _plot_spatial_expr_mrna, _plot_tissue_only, _plot_image_placeholder, IMAGE_NA_PATH, _probe_spatial_meta, compute_moran_single, compute_bivariate_moran_single

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io




URL_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'
protein_names = requests.get(f"{URL_REPO}/common_protein_31.txt").text.strip().splitlines()

def save_plot_to_s3(fig, s3_client, bucket, key):
    img_data = io.BytesIO()
    fig.savefig(img_data, format='png', bbox_inches='tight', dpi=100)
    img_data.seek(0)
    s3_client.upload_fileobj(img_data, bucket, key, ExtraArgs={'ContentType': 'image/png'})
    plt.close(fig)


def save_plot_to_s3(fig, s3_client, bucket, key):
    img_data = io.BytesIO()
    fig.savefig(img_data, format='png', bbox_inches='tight', dpi=100)
    img_data.seek(0)
    s3_client.upload_fileobj(
        img_data,
        bucket,
        key,
        ExtraArgs={'ContentType': 'image/png'}
    )
    plt.close(fig)


def cleanup_old_tasks():
    print(f"[{datetime.now()}] Running cleanup: checking for tasks older than 48 hours...")
    session = Session()
    s3 = get_s3_client()
    bucket_name = os.getenv("BUCKET_NAME")

    threshold_time = datetime.utcnow() - timedelta(hours=48)

    old_tasks = session.query(ProteinTask).filter(
        ProteinTask.created_at < threshold_time,
        ProteinTask.status.in_(['completed', 'failed'])
    ).all()

    for task in old_tasks:
        try:
            print(f"Cleaning up S3 files for task: {task.feature_code}")

            prefix = f"task_{task.feature_code}/"
            objects_to_delete = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

            if 'Contents' in objects_to_delete:
                delete_keys = {'Objects': [{'Key': obj['Key']} for obj in objects_to_delete['Contents']]}
                s3.delete_objects(Bucket=bucket_name, Delete=delete_keys)
                print(f"Deleted files for task {task.feature_code} from bucket.")

            task.status = 'deleted'
            task.output_path = None

        except Exception as e:
            print(f"Error cleaning up task {task.feature_code}: {e}")


    session.commit()
    session.close()
    print("Cleanup finished.")


def run_worker():
    print("Worker started. Monitoring MySQL...")

    last_cleanup_time = 0
    CLEANUP_INTERVAL = 1800

    while True:
        current_time = time.time()
        if current_time - last_cleanup_time > CLEANUP_INTERVAL:
            try:
                cleanup_old_tasks()
            except Exception as e:
                print(f"Cleanup failed: {e}")
            last_cleanup_time = current_time

        session = Session()
        task = session.query(ProteinTask).filter_by(status='pending').first()

        if task:
            try:
                task.status = 'running'
                session.commit()
                print(f"Processing Task: {task.feature_code}")

                s3 = get_s3_client()
                bucket = os.getenv("BUCKET_NAME")

                local_in = f"tmp_{task.feature_code}.h5ad"
                s3.download_file(bucket, task.input_path, local_in)


                adata_in = ad.read_h5ad(local_in)
                adata_origin = ad.read_h5ad(local_in)
                preprocess_ST(adata_origin)
                adata_out, missing_genes = web_predict(URL_REPO, adata_in)

                plot_prefix = f"task_{task.feature_code}/spatial_plots"
                print(f"Generating plots for {task.feature_code}...")

                has_img, lib_id, img_k, meta = _probe_spatial_meta(adata_out)
                #print(f"[{task.feature_code}] Spatial detect: has_img={has_img}, lib_id={lib_id}, img_key={img_k}")
                try:
                    fig_tissue = _plot_tissue_only(adata_out, None, None)
                    save_plot_to_s3(fig_tissue, s3, bucket, f"{plot_prefix}/tissue.png")
                except Exception as e:
                    print(f"Tissue plot failed: {e}")

                for p_name in protein_names:
                    try:
                        p_name = p_name.strip()
                        fig_p = _plot_spatial_expr(adata_out, p_name, lib_id, img_k)
                        save_plot_to_s3(fig_p, s3, bucket, f"{plot_prefix}/protein_{p_name}.png")

                        if p_name in adata_origin.var_names:
                            fig_m = _plot_spatial_expr_mrna(adata_origin, p_name, lib_id, img_k)
                        else:
                            fig_m = _plot_image_placeholder(f"{p_name}\nNot in mRNA")
                        save_plot_to_s3(fig_m, s3, bucket, f"{plot_prefix}/mrna_{p_name}.png")
                    except Exception as e:
                        print(f"Plotting skipped for {p_name} due to: {e}")

                # Protein resolutions: 0.3, 0.4, 0.5, 0.6
                for res in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]:
                    fig, ax = plt.subplots(figsize=(4.8, 4.8))
                    _plot_leiden_clustering(adata_out, ax=ax, n_neighbors=15, resolution=res,
                                            title=f"Protein Res {res}")
                    save_plot_to_s3(fig, s3, bucket, f"{plot_prefix}/leiden_prot_{res}.png")

                # mRNA resolutions: 0.8, 0.9, 1.0, 1.1
                for res in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]:
                    fig, ax = plt.subplots(figsize=(4.8, 4.8))
                    _plot_leiden_clustering(adata_origin, ax=ax, n_neighbors=15, resolution=res, title=f"mRNA Res {res}")
                    save_plot_to_s3(fig, s3, bucket, f"{plot_prefix}/leiden_mrna_{res}.png")

                local_out = f"out_{task.feature_code}.h5ad"
                adata_out.write_h5ad(local_out)

                moran_df, moran_stats = moran_df, moran_stats = compute_moran_single(
                    adata_out,
                    adata_origin,
                    tissue_name="Current Sample",
                    coord_type="grid",
                    n_perms=10
                )
                csv_buf = io.StringIO()
                moran_df.to_csv(csv_buf, index=False)
                s3.put_object(Bucket=bucket, Key=f"task_{task.feature_code}/moran_statistics.csv",
                              Body=csv_buf.getvalue())

                bivariate_df = compute_bivariate_moran_single(
                    adata_out,
                    tissue_name="Current Sample",
                    output_dir=None,
                    coord_type="grid"
                )

                csv_buf = io.StringIO()
                bivariate_df.to_csv(csv_buf, index=False)
                s3.put_object(
                    Bucket=bucket,
                    Key=f"task_{task.feature_code}/bivariate_moran_colocalization.csv",
                    Body=csv_buf.getvalue()
                )

                local_in_pre = f"pre_{task.feature_code}.h5ad"

                adata_origin.write_h5ad(local_in_pre)

                s3.upload_file(local_out, bucket, f"task_{task.feature_code}/output.h5ad")

                s3.upload_file(local_in_pre, bucket, f"task_{task.feature_code}/input_preprocessed.h5ad")

                task.output_path = f"task_{task.feature_code}/output.h5ad"
                task.status = 'completed'
                session.commit()

                send_notification(task.email, task.feature_code, success=True, note=f"The number of missing genes: {len(missing_genes)}.\n Missing genes:\n {missing_genes}" if missing_genes else "All genes present.")

                for f in [local_in, local_out, local_in_pre]:
                    if os.path.exists(f): os.remove(f)

            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"Error processing {task.feature_code}: {e}")
                session.rollback()
                task.status = 'failed'
                session.commit()
                send_notification(task.email, task.feature_code, success=False, note=str(e))

            session.close()
            time.sleep(10 if not task else 1)

        else:
            session.close()
            time.sleep(10)




if __name__ == "__main__":
    run_worker()