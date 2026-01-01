import time
import os
import anndata as ad
from task_manager import Session, ProteinTask, get_s3_client, send_notification
from dgat_utils.predict_util import web_predict

URL_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'


def run_worker():
    print("Worker started. Monitoring MySQL...")
    while True:
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

                # 3. 模拟计算
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

                # 清理本地临时文件
                os.remove(local_in)
                os.remove(local_out)

            except Exception as e:
                print(f"Error processing {task.feature_code}: {e}")
                task.status = 'failed'
                #把错误信息发送给用户
                send_notification(task.email, task.feature_code, success=False, note=str(e))


                session.commit()
                send_notification(task.email, task.feature_code, success=False)

        session.close()
        time.sleep(10)  # 轮询间隔


if __name__ == "__main__":
    run_worker()