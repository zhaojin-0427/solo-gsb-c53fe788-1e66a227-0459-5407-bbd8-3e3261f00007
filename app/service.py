"""审计任务编排：解析输入、运行引擎、按确定性排序写入结果。"""
from __future__ import annotations

import os

import pysam

from .engine import Thresholds, classify_site
from .models import Audit, SiteResult, ThresholdVersion
from .normalize import normalize_record
from .pedigree import find_trios, parse_ped
from .ploidy import PloidyRuleSet


def _ensure_fai(fasta_path: str) -> None:
    if not os.path.exists(fasta_path + ".fai"):
        pysam.faidx(fasta_path)


def run_audit(db, audit: Audit) -> None:
    """同步执行审计；任何异常都落入 audit.status/audit.error。"""
    try:
        _run(db, audit)
    except Exception as exc:  # noqa: BLE001 - 需要把任意失败持久化到审计记录
        audit.status = "failed"
        audit.error = str(exc)
    db.commit()


def _run(db, audit: Audit) -> None:
    tv = db.get(ThresholdVersion, audit.threshold_version_id)
    if tv is None:
        raise ValueError(f"阈值版本 {audit.threshold_version_id} 不存在")
    th = Thresholds(tv.min_dp, tv.min_gq, tv.max_pl_diff)
    rules = PloidyRuleSet.from_dict(audit.ploidy_rules)

    files = audit.files
    with open(files["pedigree"], encoding="utf-8") as fh:
        persons = parse_ped(fh.read())
    trios, notes = find_trios(persons)
    if not trios:
        raise ValueError("PED 中不存在完整的三人核心家系（子代双亲均须在 PED 中登记）")

    _ensure_fai(files["reference_fasta"])
    fasta = pysam.FastaFile(files["reference_fasta"])
    vcf = pysam.VariantFile(files["vcf"])

    vcf_samples = [str(s) for s in vcf.header.samples]
    contigs = [str(c) for c in vcf.header.contigs]
    needed = sorted({
        sid
        for trio in trios
        for sid in (trio.child.individual_id, trio.father.individual_id, trio.mother.individual_id)
    })
    needed_in_vcf = [s for s in needed if s in vcf_samples]

    rows = []
    for rec in vcf:
        sites, rec_notes = normalize_record(rec, fasta, needed_in_vcf)
        notes.extend(rec_notes)
        for site in sites:
            for trio in trios:
                rows.append(classify_site(site, trio, rules, th))
    vcf.close()
    fasta.close()

    # 确定性排序：contig 按 VCF 头顺序，其后依次为位置、等位、家系、子代
    contig_rank = {c: i for i, c in enumerate(contigs)}
    rows.sort(key=lambda r: (
        contig_rank.get(r["contig"], len(contigs)),
        r["pos"],
        r["ref"],
        r["alt"],
        r["family_id"],
        r["child_id"],
    ))

    counts: dict = {}
    objects = []
    for rank, row in enumerate(rows, start=1):
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
        objects.append(SiteResult(audit_id=audit.id, rank=rank, **row))
    db.add_all(objects)

    audit.status = "completed"
    audit.sample_names = vcf_samples
    audit.contigs = contigs
    audit.trio_count = len(trios)
    audit.counts = counts
    audit.exclusions = notes
