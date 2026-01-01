import time
import os
import anndata as ad
from datetime import datetime, timedelta

import requests

from task_manager import Session, ProteinTask, get_s3_client, send_notification
from dgat_utils.predict_util import web_predict

from dgat_utils.downstream import _plot_leiden_clustering, _plot_spatial_expr, _plot_spatial_expr_mrna, _plot_tissue_only, _plot_image_placeholder, IMAGE_NA_PATH, _probe_spatial_meta, compute_moran_single, compute_bivariate_moran_single

import matplotlib
matplotlib.use('Agg') # 必须在导入 pyplot 之前
import matplotlib.pyplot as plt
import io

URL_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'
protein_names = requests.get(f"{URL_REPO}/common_protein_31.txt").text.strip().splitlines()

def save_plot_to_s3(fig, s3_client, bucket, key):
    """将 Matplotlib 图表直接上传至 S3，不产生本地文件"""
    img_data = io.BytesIO()
    fig.savefig(img_data, format='png', bbox_inches='tight', dpi=100)
    img_data.seek(0)
    s3_client.upload_fileobj(img_data, bucket, key, ExtraArgs={'ContentType': 'image/png'})
    plt.close(fig)


def save_plot_to_s3(fig, s3_client, bucket, key):
    img_data = io.BytesIO()
    fig.savefig(img_data, format='png', bbox_inches='tight', dpi=100)
    img_data.seek(0)
    # 关键点：使用 upload_fileobj 处理内存中的数据流
    s3_client.upload_fileobj(
        img_data,
        bucket,
        key,
        ExtraArgs={'ContentType': 'image/png'}
    )
    plt.close(fig)  # 释放内存

def cleanup_old_tasks():
    """清理48小时前的数据：包括存储桶文件和数据库记录"""
    print(f"[{datetime.now()}] Running cleanup: checking for tasks older than 48 hours...")
    session = Session()
    s3 = get_s3_client()
    bucket_name = os.getenv("BUCKET_NAME")

    # 1. 计算48小时前的时间点

    # !!!!!!!!现在是测试所以是1 ！！！！！！************************
    threshold_time = datetime.utcnow() - timedelta(hours=12)


    # 2. 查询过期的任务
    old_tasks = session.query(ProteinTask).filter(ProteinTask.created_at < threshold_time).all()

    for task in old_tasks:
        try:
            print(f"Cleaning up task: {task.feature_code}")

            # 3. 删除存储桶中该特征码对应的所有文件 (整个目录)
            # S3没有文件夹概念，需要列出所有该前缀的文件并删除
            prefix = f"task_{task.feature_code}/"
            objects_to_delete = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

            if 'Contents' in objects_to_delete:
                delete_keys = {'Objects': [{'Key': obj['Key']} for obj in objects_to_delete['Contents']]}
                s3.delete_objects(Bucket=bucket_name, Delete=delete_keys)
                print(f"Deleted files for task {task.feature_code} from bucket.")

            # 4. 从数据库中删除记录
            session.delete(task)

        except Exception as e:
            print(f"Error cleaning up task {task.feature_code}: {e}")

    session.commit()
    session.close()
    print("Cleanup finished.")


def run_worker():
    print("Worker started. Monitoring MySQL...")

    # 初始化上次清理时间
    last_cleanup_time = 0
    # 30分钟 = 1800秒
    CLEANUP_INTERVAL = 1800

    while True:
        # --- 定时清理逻辑 ---
        current_time = time.time()
        if current_time - last_cleanup_time > CLEANUP_INTERVAL:
            try:
                cleanup_old_tasks()
            except Exception as e:
                print(f"Cleanup failed: {e}")
            last_cleanup_time = current_time

        # --- 正常任务处理逻辑 ---
        session = Session()
        task = session.query(ProteinTask).filter_by(status='pending').first()

        if task:
            try:
                task.status = 'running'
                session.commit()
                print(f"Processing Task: {task.feature_code}")

                s3 = get_s3_client()
                bucket = os.getenv("BUCKET_NAME")

                # 2. 下载
                local_in = f"tmp_{task.feature_code}.h5ad"
                s3.download_file(bucket, task.input_path, local_in)

                # 3. 核心计算
                adata_in = ad.read_h5ad(local_in)
                adata_out = web_predict(URL_REPO, adata_in)

                plot_prefix = f"task_{task.feature_code}/spatial_plots"
                print(f"Generating plots for {task.feature_code}...")

                # 4.1 绘制 Tissue 基准图

                has_img, lib_id, img_k, meta = _probe_spatial_meta(adata_out)
                #print(f"[{task.feature_code}] Spatial detect: has_img={has_img}, lib_id={lib_id}, img_key={img_k}")
                try:
                    fig_tissue = _plot_tissue_only(adata_out, None, None)
                    save_plot_to_s3(fig_tissue, s3, bucket, f"{plot_prefix}/tissue.png")
                except Exception as e:
                    print(f"Tissue plot failed: {e}")

                # 4.2 遍历 31 个蛋白质绘图
                for p_name in protein_names:
                    try:
                        p_name = p_name.strip()
                        # 绘制 Protein
                        fig_p = _plot_spatial_expr(adata_out, p_name, lib_id, img_k)
                        save_plot_to_s3(fig_p, s3, bucket, f"{plot_prefix}/protein_{p_name}.png")

                        # 绘制 mRNA
                        if p_name in adata_in.var_names:
                            fig_m = _plot_spatial_expr_mrna(adata_in, p_name, lib_id, img_k)
                        else:
                            fig_m = _plot_image_placeholder(f"{p_name}\nNot in mRNA")
                        save_plot_to_s3(fig_m, s3, bucket, f"{plot_prefix}/mrna_{p_name}.png")
                    except Exception as e:
                        print(f"Plotting skipped for {p_name} due to: {e}")

                # 4.3 预渲染 Leiden 聚类图 (使用固定参数)
                # Protein resolutions: 0.3, 0.4, 0.5, 0.6
                for res in [0.3, 0.4, 0.5, 0.6]:
                    fig, ax = plt.subplots(figsize=(4.8, 4.8))
                    # n_neighbors 设置为 15 (Scanpy 默认，0 是无效的)
                    _plot_leiden_clustering(adata_out, ax=ax, n_neighbors=15, resolution=res,
                                            title=f"Protein Res {res}")
                    save_plot_to_s3(fig, s3, bucket, f"{plot_prefix}/leiden_prot_{res}.png")

                # mRNA resolutions: 0.8, 0.9, 1.0, 1.1
                for res in [0.8, 0.9, 1.0, 1.1]:
                    fig, ax = plt.subplots(figsize=(4.8, 4.8))
                    _plot_leiden_clustering(adata_in, ax=ax, n_neighbors=15, resolution=res, title=f"mRNA Res {res}")
                    save_plot_to_s3(fig, s3, bucket, f"{plot_prefix}/leiden_mrna_{res}.png")

                local_out = f"out_{task.feature_code}.h5ad"
                adata_out.write_h5ad(local_out)

                moran_df, moran_stats = moran_df, moran_stats = compute_moran_single(
                    adata_out,
                    adata_in,
                    tissue_name="Current Sample",
                    coord_type="grid",
                    n_perms=10
                )
                # 将结果保存并上传到 S3
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

                # --- 【关键修复点 1】: adata_in 也需要写入本地再上传 ---
                local_in_pre = f"pre_{task.feature_code}.h5ad"
                adata_in.write_h5ad(local_in_pre)

                s3.upload_file(local_out, bucket, f"task_{task.feature_code}/output.h5ad")
                s3.upload_file(local_in_pre, bucket, f"task_{task.feature_code}/input_preprocessed.h5ad")

                # 6. 更新状态
                task.output_path = f"task_{task.feature_code}/output.h5ad"
                task.status = 'completed'
                session.commit()

                # 发送通知
                send_notification(task.email, task.feature_code, success=True)

                # 7. 清理所有本地临时文件
                for f in [local_in, local_out, local_in_pre]:
                    if os.path.exists(f): os.remove(f)

            except Exception as e:
                import traceback
                traceback.print_exc()  # 打印详细堆栈到 Railway 日志
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