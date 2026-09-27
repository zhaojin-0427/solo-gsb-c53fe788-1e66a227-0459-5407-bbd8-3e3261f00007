"""SQLAlchemy 数据模型。"""
from __future__ import annotations

from datetime import datetime

from sqlalchemy import JSON, DateTime, ForeignKey, String, Text, UniqueConstraint, func
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column


class Base(DeclarativeBase):
    pass


class ThresholdVersion(Base):
    """不可变阈值版本：仅可创建与读取，不提供任何修改/删除接口。"""

    __tablename__ = "threshold_versions"
    __table_args__ = (
        UniqueConstraint("name", "min_dp", "min_gq", "max_pl_diff", name="uq_threshold_version"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(100))
    min_dp: Mapped[int]
    min_gq: Mapped[int]
    max_pl_diff: Mapped[int]
    description: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime, server_default=func.now())


class Audit(Base):
    __tablename__ = "audits"

    id: Mapped[int] = mapped_column(primary_key=True)
    threshold_version_id: Mapped[int] = mapped_column(ForeignKey("threshold_versions.id"))
    status: Mapped[str] = mapped_column(String(20), default="running")  # running/completed/failed
    error: Mapped[str | None] = mapped_column(Text, nullable=True)
    replay_of_audit_id: Mapped[int | None] = mapped_column(nullable=True)
    ploidy_rules: Mapped[dict] = mapped_column(JSON, default=dict)
    files: Mapped[dict] = mapped_column(JSON, default=dict)
    sample_names: Mapped[list] = mapped_column(JSON, default=list)
    contigs: Mapped[list] = mapped_column(JSON, default=list)
    trio_count: Mapped[int] = mapped_column(default=0)
    counts: Mapped[dict] = mapped_column(JSON, default=dict)
    exclusions: Mapped[list] = mapped_column(JSON, default=list)
    created_at: Mapped[datetime] = mapped_column(DateTime, server_default=func.now())


class SiteResult(Base):
    __tablename__ = "site_results"

    id: Mapped[int] = mapped_column(primary_key=True)
    audit_id: Mapped[int] = mapped_column(ForeignKey("audits.id"), index=True)
    rank: Mapped[int]  # 确定性排序序号（1 起）
    contig: Mapped[str] = mapped_column(String(64), index=True)
    pos: Mapped[int]
    ref: Mapped[str] = mapped_column(Text)
    alt: Mapped[str] = mapped_column(Text)
    family_id: Mapped[str] = mapped_column(String(128))
    child_id: Mapped[str] = mapped_column(String(128), index=True)
    region_class: Mapped[str] = mapped_column(String(20))
    verdict: Mapped[str] = mapped_column(String(32), index=True)
    detail: Mapped[dict] = mapped_column(JSON, default=dict)
