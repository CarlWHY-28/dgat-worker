import os
import boto3
import smtplib
from email.mime.text import MIMEText
from sqlalchemy import create_engine, Column, String, DateTime, Enum
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import socket
import requests

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

# def send_notification(email, code, success=True, note = ''):
#     subject = "DGAT Imputation Finished" if success else "DGAT Failed"
#     body = f"\nYour feature code is：{code}\nPlease go to the website -> Your Results to check the result."
#     if not success:
#         body = f"Your task failed during processing, please check the file format."
#         if note:
#             body += f"\nError info: {note}"
#
#     msg = MIMEText(body)
#     msg['Subject'] = subject
#     msg['From'] = os.getenv("SMTP_USER")
#     msg['To'] = email
#
#     try:
#         # 1. 使用 587 端口建立普通 SMTP 连接
#         server = smtplib.SMTP('smtp.gmail.com', 587, timeout=30)
#         # 2. 将连接升级为加密连接 (必须)
#         server.starttls()
#         # 3. 登录并发送
#         server.login('carlwanghy@gmail.com', 'oilucthfcbwrykjf')
#         server.send_message(msg)
#         server.quit()
#     except Exception as e:
#         print(f"Mail Error: {e}")


# def send_notification(email, code, success=True, note=''):
#     subject = "DGAT Imputation Finished" if success else "DGAT Failed"
#     body = f"Your feature code is: {code}\n..."
#
#     msg = MIMEText(body)
#     msg['Subject'] = subject
#     msg['From'] = os.getenv("SMTP_USER")
#     msg['To'] = email
#
#     smtp_server = os.getenv("SMTP_SERVER")
#     smtp_port = int(os.getenv("SMTP_PORT", 587))
#     smtp_user = os.getenv("SMTP_USER")
#     smtp_pass = os.getenv("SMTP_PASS")
#
#     try:
#         # 如果是 587 端口，必须使用以下逻辑
#         if smtp_port == 587:
#             server = smtplib.SMTP(smtp_server, smtp_port, timeout=60)
#             server.starttls()  # 关键：升级为安全连接
#         else:
#             # 如果是 465 端口，直接使用 SSL
#             server = smtplib.SMTP_SSL(smtp_server, smtp_port)
#
#         server.login(smtp_user, smtp_pass)
#         server.send_message(msg)
#         server.quit()
#         print(f"Email sent to {email}")
#     except Exception as e:
#         print(f"Mail Error: {e}")
#
# def send_notification(email, code, success=True, note=''):
#     subject = "DGAT Imputation Finished" if success else "DGAT Failed"
#     body = f"Your feature code is: {code}\n{note}"
#
#     msg = MIMEText(body)
#     msg['Subject'] = subject
#     msg['From'] = os.getenv("SMTP_USER")
#     msg['To'] = email
#
#     # Gmail 配置
#     smtp_server = "smtp.gmail.com"
#     # 强烈建议使用 465 (SSL)，比 587 更稳定
#     smtp_port = 465
#
#     smtp_user = os.getenv("SMTP_USER")
#     smtp_pass = os.getenv("SMTP_PASS")
#
#     try:
#         # 使用 SMTP_SSL 直接连接 465 端口
#         # 注意：这里不需要再调用 starttls()，因为连接建立时已经是加密的了
#         server = smtplib.SMTP_SSL(smtp_server, smtp_port, timeout=60)
#
#         server.login(smtp_user, smtp_pass)
#         server.send_message(msg)
#         server.quit()
#         print(f"Email sent successfully to {email}")
#
#     except Exception as e:
#         print(f"Mail Error: {e}")
#
# #
# import os
# import resend
def send_notification(email, code, success=True, note=''):
    # 从环境变量获取配置
    api_key = os.getenv("BREVO_API_KEY")
    sender_email = os.getenv("SENDER_EMAIL")  # 必须是你在 Brevo 验证过的邮箱

    if not api_key or not sender_email:
        print("Brevo config missing")
        return

    # 准备邮件内容
    subject = "DGAT Imputation Finished" if success else "DGAT Failed"
    status_color = "green" if success else "red"
    status_text = "Success" if success else "Failed"

    # HTML 内容（可以写得很漂亮）
    html_content = f"""
    <html>
    <body>
        <h2>DGAT Task Notification</h2>
        <p>Status: <strong style="color: {status_color};">{status_text}</strong></p>
        <p>Your feature code is: <strong>{code}</strong></p>
        <p>Note: {note}</p>
        <hr>
        <p style="font-size: 12px; color: gray;">Sent from DGAT System</p>
    </body>
    </html>
    """

    # Brevo API 地址
    url = "https://api.brevo.com/v3/smtp/email"

    # 请求头
    headers = {
        "accept": "application/json",
        "api-key": api_key,
        "content-type": "application/json"
    }

    # 请求体
    payload = {
        "sender": {"name": "DGAT System", "email": sender_email},
        "to": [{"email": email}],  # 这里的 email 可以是任何真实用户的邮箱
        "subject": subject,
        "htmlContent": html_content
    }

    try:
        response = requests.post(url, json=payload, headers=headers, timeout=10)

        if response.status_code == 201:
            print(f"Email sent successfully to {email}")
        else:
            # 打印错误详情以便调试
            print(f"Brevo Error: {response.status_code} - {response.text}")

    except Exception as e:
        print(f"Network Error: {e}")
#
# def send_notification(email, code, success=True, note=''):
#     # 从环境变量获取配置
#     api_key = os.getenv("BREVO_API_KEY")
#     sender_email = os.getenv("SENDER_EMAIL")  # 必须是你在 Brevo 验证过的邮箱
#
#     if not api_key or not sender_email:
#         print("Brevo config missing")
#         return
#
#     # 准备邮件内容
#     subject = "DGAT Imputation Finished" if success else "DGAT Failed"
#     status_color = "green" if success else "red"
#     status_text = "Success" if success else "Failed"
#
#     # HTML 内容（可以写得很漂亮）
#     html_content = f"""
#     <html>
#     <body>
#         <h2>DGAT Task Notification</h2>
#         <p>Status: <strong style="color: {status_color};">{status_text}</strong></p>
#         <p>Your feature code is: <strong>{code}</strong></p>
#         <p>Note: {note}</p>
#         <hr>
#         <p style="font-size: 12px; color: gray;">Sent from DGAT System</p>
#     </body>
#     </html>
#     """
#
#     # Brevo API 地址
#     url = "https://api.brevo.com/v3/smtp/email"
#
#     # 请求头
#     headers = {
#         "accept": "application/json",
#         "api-key": api_key,
#         "content-type": "application/json"
#     }
#
#     # 请求体
#     payload = {
#         "sender": {"name": "DGAT System", "email": sender_email},
#         "to": [{"email": email}],  # 这里的 email 可以是任何真实用户的邮箱
#         "subject": subject,
#         "htmlContent": html_content
#     }
#
#     try:
#         response = requests.post(url, json=payload, headers=headers, timeout=10)
#
#         if response.status_code == 201:
#             print(f"Email sent successfully to {email}")
#         else:
#             # 打印错误详情以便调试
#             print(f"Brevo Error: {response.status_code} - {response.text}")
#
#     except Exception as e:
#         print(f"Network Error: {e}")