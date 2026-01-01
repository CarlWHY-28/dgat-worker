import time
import os
import anndata as ad
from datetime import datetime, timedelta

import requests

from task_manager import Session, ProteinTask, get_s3_client, send_notification
from dgat_utils.predict_util import web_predict

from dgat_utils.downstream import _plot_leiden_clustering, _plot_spatial_expr, _plot_spatial_expr_mrna, _plot_tissue_only, _plot_image_placeholder, IMAGE_NA_PATH

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
                # 1. 锁定任务
                task.status = 'running'
                session.commit()
                print(f"Processing Task: {task.feature_code}")

                # 2. 从 Bucket 下载
                s3 = get_s3_client()
                bucket = os.getenv("BUCKET_NAME")
                local_in = f"tmp_{task.feature_code}.h5ad"
                s3.download_file(bucket, task.input_path, local_in)

                # 3. 核心计算
                adata_in = ad.read_h5ad(local_in)
                adata_out = web_predict(URL_REPO, adata_in)

                s3 = get_s3_client()
                bucket = os.getenv("BUCKET_NAME")
                plot_prefix = f"task_{task.feature_code}/spatial_plots"

                print(f"Generating plots for {task.feature_code}...")

                # 4.1 绘制 Tissue 基准图
                fig_tissue = _plot_tissue_only(adata_out, None, None)  # 使用你之前的函数
                save_plot_to_s3(fig_tissue, s3, bucket, f"{plot_prefix}/tissue.png")

                # 4.2 遍历 31 个蛋白质绘图 (Protein + mRNA)
                for p_name in protein_names:
                    # 绘制 Protein (plasma)
                    fig_p = _plot_spatial_expr(adata_out, p_name, None, None)
                    save_plot_to_s3(fig_p, s3, bucket, f"{plot_prefix}/protein_{p_name}.png")

                    # 绘制 mRNA (viridis)
                    if p_name in adata_in.var_names:
                        fig_m = _plot_spatial_expr_mrna(adata_in, p_name, None, None)
                    else:
                        fig_m = _plot_image_placeholder(IMAGE_NA_PATH)
                    save_plot_to_s3(fig_m, s3, bucket, f"{plot_prefix}/mrna_{p_name}.png")

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

                # 4. 上传结果 (task_{code}/output.h5ad)
                local_out = f"out_{task.feature_code}.h5ad"
                adata_out.write_h5ad(local_out)
                output_key = f"task_{task.feature_code}/output.h5ad"
                s3.upload_file(local_out, bucket, output_key)
                s3.upload_file(adata_in, bucket, f"task_{task.feature_code}/input_preprocessed.h5ad")

                # 5. 更新状态并发送邮件
                task.output_path = output_key
                task.status = 'completed'
                session.commit()
                send_notification(task.email, task.feature_code, success=True)

                # 清理本地临时文件
                if os.path.exists(local_in): os.remove(local_in)
                if os.path.exists(local_out): os.remove(local_out)

            except Exception as e:
                print(f"Error processing {task.feature_code}: {e}")
                session.rollback()
                task.status = 'failed'
                session.commit()

                send_notification(task.email, task.feature_code, success=False, note=str(e))

        session.close()
        time.sleep(10 if not task else 1)


if __name__ == "__main__":
    run_worker()