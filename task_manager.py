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
    subject = "DGAT Imputation Finished" if success else "DGAT Failed"
    body = f"\nYour feature code is：{code}\nPlease go to the website -> Your Results to check the result."
    if not success:
        body = f"Your task failed during processing, please check the file format."
        if note:
            body += f"\nError info: {note}"

    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = os.getenv("SMTP_USER")
    msg['To'] = email

    try:
        # 1. 使用 587 端口建立普通 SMTP 连接
        server = smtplib.SMTP('smtp.gmail.com', 587, timeout=30)
        # 2. 将连接升级为加密连接 (必须)
        server.starttls()
        # 3. 登录并发送
        server.login('carlwanghy@gmail.com', 'oilucthfcbwrykjf')
        server.send_message(msg)
        server.quit()
    except Exception as e:
        print(f"Mail Error: {e}")