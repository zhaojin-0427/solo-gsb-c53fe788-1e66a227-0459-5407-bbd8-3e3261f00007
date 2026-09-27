"""Pydantic 请求/响应模式。"""
from __future__ import annotations

from datetime import datetime

from pydantic import BaseModel, ConfigDict, Field


class ThresholdVersionCreate(BaseModel):
    name: str = Field(min_length=1, max_length=100)
    min_dp: int = Field(ge=0, description="最小深度：样本 DP 须达到该值")
    min_gq: int = Field(ge=0, description="最小基因型质量：样本 GQ 须达到该值")
    max_pl_diff: int = Field(ge=0, description="最大 PL 差：亲本等位支持罚分上限")
    description: str = ""


class ThresholdVersionOut(ThresholdVersionCreate):
    model_config = ConfigDict(from_attributes=True)

    id: int
    created_at: datetime


class ReplayRequest(BaseModel):
    threshold_version_id: int


class AuditOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    status: str
    error: str | None
    threshold_version_id: int
    replay_of_audit_id: int | None
    trio_count: int
    sample_names: list
    counts: dict
    exclusions: list
    created_at: datetime


class SiteResultOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    rank: int
    contig: str
    pos: int
    ref: str
    alt: str
    family_id: str
    child_id: str
    region_class: str
    verdict: str
    detail: dict


class ResultPage(BaseModel):
    total: int
    items: list[SiteResultOut]
