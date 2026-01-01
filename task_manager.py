import os
import boto3
import smtplib
from email.mime.text import MIMEText
from sqlalchemy import create_engine, Column, String, DateTime, Enum
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from datetime import datetime

# --- 配置初始化 ---
Base = declarative_base()

# 数据库连接：处理 Railway 的 mysql:// 格式mysql://root:jdakpxILhrVFxpLzuSxCjXlZFluTLGAR@mysql.railway.internal:3306/railway

db_url = os.getenv("MYSQL_URL", "mysql://root:PvklhrFKKovSXiFgwRqvaqifIrwnhpdS@mysql.railway.internal:3306/railway")
if db_url.startswith("mysql://"):
    db_url = db_url.replace("mysql://", "mysql+pymysql://", 1)

engine = create_engine(db_url, pool_recycle=3600)
Session = sessionmaker(bind=engine)

# --- 数据库模型 ---
class ProteinTask(Base):
    __tablename__ = 'protein_tasks'
    feature_code = Column(String(36), primary_key=True)
    email = Column(String(255))
    status = Column(Enum('pending', 'running', 'completed', 'failed'), default='pending')
    input_path = Column(String(500))
    output_path = Column(String(500))
    created_at = Column(DateTime, default=datetime.utcnow)

# 自动创建表
Base.metadata.create_all(engine)

# --- 工具函数 ---
def get_s3_client():
    return boto3.client(
        's3',
        endpoint_url=os.getenv("S3_ENDPOINT"),
        aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
        aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        region_name="auto"
    )

def send_notification(email, code, success=True, note = ''):
    subject = "DGAT 推理任务完成通知" if success else "DGAT 推理任务失败"
    body = f"您的任务已处理完成！\n取件码为：{code}\n请前往网站输入此码查询结果。"
    if not success:
        body = f"抱歉，任务 {code} 处理过程中出现错误，请检查文件格式。"
        if note:
            body += f"\n错误信息：{note}"

    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = os.getenv("SMTP_USER")
    msg['To'] = email

    try:
        with smtplib.SMTP_SSL(os.getenv("SMTP_SERVER"), int(os.getenv("SMTP_PORT", 465))) as server:
            server.login(os.getenv("SMTP_USER"), os.getenv("SMTP_PASS"))
            server.send_message(msg)
    except Exception as e:
        print(f"Mail Error: {e}")