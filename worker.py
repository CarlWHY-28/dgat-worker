import time
import os
import anndata as ad
from datetime import datetime, timedelta
from task_manager import Session, ProteinTask, get_s3_client, send_notification
from dgat_utils.predict_util import web_predict

URL_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'


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


                # 4. 上传结果 (task_{code}/output.h5ad)
                local_out = f"out_{task.feature_code}.h5ad"
                adata_out.write_h5ad(local_out)
                output_key = f"task_{task.feature_code}/output.h5ad"
                s3.upload_file(local_out, bucket, output_key)

                # 5. 更新状态并发送邮件
                task.output_path = output_key
                task.status = 'completed'
                session.commit()
                send_notification(task.email, task.feature_code, success=True)

                # 6. 把adata_in(这里是预处理过的)保存到s3备用
                local_in_pre = f"pre_{task.feature_code}.h5ad"
                s3.upload_file(local_in, bucket, f"task_{task.feature_code}/input_preprocessed.h5ad")



                # 清理本地临时文件
                if os.path.exists(local_in): os.remove(local_in)
                if os.path.exists(local_out): os.remove(local_out)

            except Exception as e:
                print(f"Error processing {task.feature_code}: {e}")
                # 标记失败并回滚其他可能的更改
                session.rollback()
                task.status = 'failed'
                session.commit()

                # 发送带有详细错误信息的通知
                send_notification(task.email, task.feature_code, success=False, note=str(e))

        session.close()
        # 轮询间隔：如果有任务则快速进入下一轮，没任务休眠10秒
        time.sleep(10 if not task else 1)


if __name__ == "__main__":
    run_worker()