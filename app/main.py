"""FastAPI 入口与路由。"""
from __future__ import annotations

import json
import os
import shutil
import time
from contextlib import asynccontextmanager

from fastapi import Depends, FastAPI, File, Form, HTTPException, Query, UploadFile
from sqlalchemy.exc import IntegrityError, OperationalError
from sqlalchemy.orm import Session

from . import config, schemas
from .database import engine, get_db
from .models import Audit, Base, SiteResult, ThresholdVersion
from .ploidy import PloidyRuleSet
from .service import run_audit

DISCLAIMER = (
    "本服务仅基于客户端提交的数据做孟德尔传递一致性审计；"
    "不进行任何联网注释，不输出疾病诊断或临床建议。"
)


@asynccontextmanager
async def lifespan(app: FastAPI):
    # 等待数据库就绪（Docker Compose 下 PostgreSQL 启动需要时间）
    for _ in range(60):
        try:
            Base.metadata.create_all(engine)
            break
        except OperationalError:
            time.sleep(1)
    else:
        Base.metadata.create_all(engine)
    yield


app = FastAPI(
    title="家系 VCF 遗传一致性审计 API",
    version="1.0.0",
    description=DISCLAIMER,
    lifespan=lifespan,
)


@app.get("/")
def root():
    return {
        "service": "家系 VCF 遗传一致性审计 API",
        "docs": "/docs",
        "disclaimer": DISCLAIMER,
    }


@app.get("/health")
def health():
    return {"status": "ok"}


# ---------------------------------------------------------------- 阈值版本


@app.post("/threshold-versions", response_model=schemas.ThresholdVersionOut, status_code=201)
def create_threshold_version(payload: schemas.ThresholdVersionCreate,
                             db: Session = Depends(get_db)):
    """登记不可变阈值版本（最小 DP/GQ、最大 PL 差）。创建后不可修改。"""
    tv = ThresholdVersion(**payload.model_dump())
    db.add(tv)
    try:
        db.commit()
    except IntegrityError:
        db.rollback()
        raise HTTPException(409, "相同名称与参数组合的阈值版本已存在")
    db.refresh(tv)
    return tv


@app.get("/threshold-versions", response_model=list[schemas.ThresholdVersionOut])
def list_threshold_versions(db: Session = Depends(get_db)):
    return db.query(ThresholdVersion).order_by(ThresholdVersion.id).all()


@app.get("/threshold-versions/{tv_id}", response_model=schemas.ThresholdVersionOut)
def get_threshold_version(tv_id: int, db: Session = Depends(get_db)):
    tv = db.get(ThresholdVersion, tv_id)
    if tv is None:
        raise HTTPException(404, "阈值版本不存在")
    return tv


# ---------------------------------------------------------------- 审计任务


def _save_upload(dest_dir: str, key: str, upload: UploadFile) -> str:
    filename = os.path.basename(upload.filename or key)
    path = os.path.join(dest_dir, filename)
    with open(path, "wb") as out:
        shutil.copyfileobj(upload.file, out)
    return path


@app.post("/audits", response_model=schemas.AuditOut, status_code=201)
def create_audit(
    threshold_version_id: int = Form(...),
    ploidy_rules: str = Form(..., description="性染色体倍性规则 JSON 字符串"),
    reference_fasta: UploadFile = File(...),
    pedigree: UploadFile = File(...),
    vcf: UploadFile = File(...),
    db: Session = Depends(get_db),
):
    """提交审计：本地参考 FASTA + PED 家系 + 倍性规则 + 多样本 VCF。"""
    tv = db.get(ThresholdVersion, threshold_version_id)
    if tv is None:
        raise HTTPException(404, "阈值版本不存在")
    try:
        rules_obj = json.loads(ploidy_rules)
        PloidyRuleSet.from_dict(rules_obj)  # 提前校验
    except (json.JSONDecodeError, ValueError) as exc:
        raise HTTPException(422, f"倍性规则无效：{exc}")

    audit = Audit(
        threshold_version_id=tv.id,
        ploidy_rules=rules_obj,
        status="running",
    )
    db.add(audit)
    db.flush()  # 取得 audit.id

    dest = os.path.join(config.UPLOAD_DIR, f"audit_{audit.id}")
    os.makedirs(dest, exist_ok=True)
    audit.files = {
        "reference_fasta": _save_upload(dest, "reference_fasta", reference_fasta),
        "pedigree": _save_upload(dest, "pedigree", pedigree),
        "vcf": _save_upload(dest, "vcf", vcf),
    }
    run_audit(db, audit)
    db.refresh(audit)
    return audit


@app.get("/audits", response_model=list[schemas.AuditOut])
def list_audits(db: Session = Depends(get_db)):
    return db.query(Audit).order_by(Audit.id).all()


@app.get("/audits/{audit_id}", response_model=schemas.AuditOut)
def get_audit(audit_id: int, db: Session = Depends(get_db)):
    audit = db.get(Audit, audit_id)
    if audit is None:
        raise HTTPException(404, "审计不存在")
    return audit


@app.get("/audits/{audit_id}/results", response_model=schemas.ResultPage)
def get_audit_results(
    audit_id: int,
    verdict: str | None = Query(None, description="按判定结论过滤"),
    child_id: str | None = Query(None, description="按子代样本过滤"),
    family_id: str | None = Query(None, description="按家系过滤"),
    contig: str | None = Query(None, description="按 contig 过滤"),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
):
    """按确定性排序（rank 升序）分页返回位点级审计结果。"""
    if db.get(Audit, audit_id) is None:
        raise HTTPException(404, "审计不存在")
    q = db.query(SiteResult).filter(SiteResult.audit_id == audit_id)
    if verdict:
        q = q.filter(SiteResult.verdict == verdict)
    if child_id:
        q = q.filter(SiteResult.child_id == child_id)
    if family_id:
        q = q.filter(SiteResult.family_id == family_id)
    if contig:
        q = q.filter(SiteResult.contig == contig)
    total = q.count()
    items = q.order_by(SiteResult.rank).offset(offset).limit(limit).all()
    return {"total": total, "items": items}


@app.post("/audits/{audit_id}/replay", response_model=schemas.AuditOut, status_code=201)
def replay_audit(audit_id: int, payload: schemas.ReplayRequest,
                 db: Session = Depends(get_db)):
    """按另一个阈值版本重放已有审计（复用其上传文件与倍性规则）。"""
    original = db.get(Audit, audit_id)
    if original is None:
        raise HTTPException(404, "审计不存在")
    if db.get(ThresholdVersion, payload.threshold_version_id) is None:
        raise HTTPException(404, "阈值版本不存在")
    audit = Audit(
        threshold_version_id=payload.threshold_version_id,
        replay_of_audit_id=original.id,
        ploidy_rules=original.ploidy_rules,
        files=original.files,
        status="running",
    )
    db.add(audit)
    db.flush()
    run_audit(db, audit)
    db.refresh(audit)
    return audit
