import os
import tempfile

# 须在导入 app 之前设置，确保测试使用独立的临时数据目录与 SQLite 库
_TMP = tempfile.mkdtemp(prefix="pedigree_audit_test_")
os.environ["DATA_DIR"] = _TMP
os.environ["UPLOAD_DIR"] = os.path.join(_TMP, "uploads")
os.environ["DATABASE_URL"] = "sqlite:///" + os.path.join(_TMP, "test.db")

import pytest  # noqa: E402


@pytest.fixture()
def client():
    from app.database import engine
    from app.models import Base

    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)

    from fastapi.testclient import TestClient

    from app.main import app

    with TestClient(app) as test_client:
        yield test_client
