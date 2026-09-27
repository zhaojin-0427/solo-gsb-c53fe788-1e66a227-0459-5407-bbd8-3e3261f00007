"""运行配置（通过环境变量覆盖）。"""
from __future__ import annotations

import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.environ.get("DATA_DIR", os.path.join(BASE_DIR, "data"))
UPLOAD_DIR = os.environ.get("UPLOAD_DIR", os.path.join(DATA_DIR, "uploads"))
DATABASE_URL = os.environ.get(
    "DATABASE_URL", "sqlite:///" + os.path.join(DATA_DIR, "audit.db")
)

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(UPLOAD_DIR, exist_ok=True)
