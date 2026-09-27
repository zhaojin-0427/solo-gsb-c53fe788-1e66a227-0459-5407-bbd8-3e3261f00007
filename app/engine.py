"""核心审计引擎：对规范化后的二等位位点执行家系传递一致性判定。

阈值版本（不可变）包含三个参数：

* ``min_dp``       最小深度：样本 DP 须 >= min_dp；
* ``min_gq``       最小基因型质量：样本 GQ 须 >= min_gq（GQ 缺失时由 PL 推导，
                   即 PL 次小值与最小值之差）；
* ``max_pl_diff``  最大 PL 差：判定亲本是否"支持"某等位时，含该等位基因型的
                   最小 PL 与亲本最小 PL 之差须 <= max_pl_diff 方视为有支持。

判定结论（verdict）：

* ``consistent``          一致：存在符合区域传递规则的等位来源组合；
* ``candidate_de_novo``   候选新生变异：子代高置信等位在双亲均无支持，且其余等位可传递；
* ``mendelian_conflict``  孟德尔冲突：证据充分但无可行传递，且不满足候选新生变异条件；
* ``insufficient``        信息不足：必要样本无合格基因型；
* ``inapplicable``        不适用：区域/性别组合不适用（如女性子代的 Y 非 PAR 区）。
"""
from __future__ import annotations

from dataclasses import dataclass

from .normalize import NormalizedSite, SampleCall
from .pedigree import Trio
from .ploidy import PloidyRuleSet

CONSISTENT = "consistent"
CANDIDATE_DE_NOVO = "candidate_de_novo"
MENDELIAN_CONFLICT = "mendelian_conflict"
INSUFFICIENT = "insufficient"
INAPPLICABLE = "inapplicable"

VERDICTS = (CONSISTENT, CANDIDATE_DE_NOVO, MENDELIAN_CONFLICT, INSUFFICIENT, INAPPLICABLE)

ROLE_ZH = {"father": "父亲", "mother": "母亲", "child": "子代"}


@dataclass(frozen=True)
class Thresholds:
    min_dp: int
    min_gq: int
    max_pl_diff: int


def fmt_gt(gt) -> str:
    return "/".join(str(a) for a in gt)


def pl_derived_gq(pl):
    """由 PL 推导 GQ：次小 PL 与最小 PL 之差。"""
    if not pl or len(pl) < 2:
        return None
    ordered = sorted(pl)
    return ordered[1] - ordered[0]


def effective_gq(call: SampleCall):
    if call.gq is not None:
        return call.gq
    return pl_derived_gq(call.pl)


def eligible_genotypes(call, th: Thresholds, ploidy: int):
    """返回 (合格基因型列表, 排除理由列表)。

    合格基因型：通过 DP/GQ 阈值且倍性符合区域要求的检出基因型（按原样采用，
尽可能的确定性；杂合不做向纯合的扩张，亲本不确定性由 max_pl_diff 支持度检验刻画）。
    """
    if call is None:
        return [], ["样本不在 VCF 中"]
    if call.gt is None:
        return [], ["基因型缺失，或含有其他交替等位（多等位拆分后无法映射）"]
    if len(call.gt) != ploidy:
        return [], [f"倍性不符：该区域期望倍性 {ploidy}，实际基因型倍性 {len(call.gt)}"]
    reasons = []
    if call.dp is None:
        reasons.append("DP 缺失，无法验证最小深度")
    elif call.dp < th.min_dp:
        reasons.append(f"DP={call.dp} 低于最小深度 {th.min_dp}")
    gq = effective_gq(call)
    if gq is None:
        reasons.append("GQ 缺失且无法由 PL 推导")
    elif gq < th.min_gq:
        reasons.append(f"GQ={gq} 低于最小基因型质量 {th.min_gq}")
    if reasons:
        return [], reasons
    return [tuple(sorted(call.gt))], []


def allele_support_penalty(pl, allele: int):
    """亲本对某等位的 PL 支持罚分：含该等位基因型的最小 PL 与整体最小 PL 之差。

    返回 None 表示无 PL 数据（按"无支持"处理）。
    """
    if pl is None:
        return None
    if len(pl) == 3:  # 二倍体二等位：基因型顺序 0/0, 0/1, 1/1
        idx = (0, 1) if allele == 0 else (1, 2)
    elif len(pl) == 2:  # 单倍体二等位
        idx = (0,) if allele == 0 else (1,)
    else:
        return None
    return min(pl[i] for i in idx) - min(pl)


def _result(site: NormalizedSite, base: dict, verdict: str, participants: list,
            sources: list, evidence, reasons: list) -> dict:
    detail = {
        "verdict": verdict,
        "region_class": base["region_class"],
        "trio": base["trio"],
        "participants": participants,
        "allelic_sources": sources,
        "de_novo_evidence": evidence,
        "exclusion_reasons": reasons,
    }
    return {
        "contig": site.contig,
        "pos": site.pos0 + 1,
        "ref": site.ref,
        "alt": site.alt,
        "family_id": base["trio"]["family_id"],
        "child_id": base["trio"]["child"],
        "region_class": base["region_class"],
        "verdict": verdict,
        "detail": detail,
    }


def classify_site(site: NormalizedSite, trio: Trio, rules: PloidyRuleSet,
                  th: Thresholds) -> dict:
    """对单个规范化位点与单个三人核心家系给出确定性判定。"""
    kind = rules.region_kind(site.contig, site.pos0)
    child_sex = rules.sex_of(trio.child.sex)
    child_id = trio.child.individual_id
    father_id = trio.father.individual_id
    mother_id = trio.mother.individual_id
    base = {
        "region_class": kind,
        "trio": {
            "family_id": trio.family_id,
            "child": child_id,
            "father": father_id,
            "mother": mother_id,
        },
    }

    # ---- 适用性 ----
    if kind == "y_nonpar" and child_sex != "male":
        who = "女性" if child_sex == "female" else "性别未知的"
        return _result(site, base, INAPPLICABLE, [], [], None,
                       [f"Y 非 PAR 区不适用于{who}子代"])
    if kind == "x_nonpar" and child_sex == "unknown":
        return _result(site, base, INAPPLICABLE, [], [], None,
                       ["子代性别未知，无法确定 X 非 PAR 区倍性"])

    # ---- 参与样本与期望倍性（顺序：亲本在前，子代最后）----
    if kind == "mt":
        parts = [("mother", mother_id, 1), ("child", child_id, 1)]
    elif kind == "y_nonpar":
        parts = [("father", father_id, 1), ("child", child_id, 1)]
    elif kind == "x_nonpar" and child_sex == "male":
        parts = [("mother", mother_id, 2), ("child", child_id, 1)]
    elif kind == "x_nonpar":  # 女性子代
        parts = [("father", father_id, 1), ("mother", mother_id, 2), ("child", child_id, 2)]
    else:  # autosome / par
        parts = [("father", father_id, 2), ("mother", mother_id, 2), ("child", child_id, 2)]

    # ---- 合格基因型 ----
    participants = []
    eligible = {}
    insufficient = []
    for role, sid, ploidy in parts:
        call = site.calls.get(sid)
        elig, reasons = eligible_genotypes(call, th, ploidy)
        eligible[role] = elig
        participants.append({
            "sample": sid,
            "role": role,
            "expected_ploidy": ploidy,
            "genotype": fmt_gt(call.gt) if call is not None and call.gt else None,
            "dp": call.dp if call is not None else None,
            "gq": effective_gq(call) if call is not None else None,
            "eligible_genotypes": [fmt_gt(g) for g in elig],
            "exclusion_reasons": reasons,
        })
        if reasons:
            insufficient.append(f"样本 {sid}（{ROLE_ZH[role]}）：" + "；".join(reasons))
    if insufficient:
        return _result(site, base, INSUFFICIENT, participants, [], None, insufficient)

    child_ploidy = parts[-1][2]
    called_gt = eligible["child"][0]

    # ---- 传递可行性枚举（等位来源）----
    sources = []
    if child_ploidy == 2:
        father_alleles = sorted({a for gt in eligible["father"] for a in gt})
        mother_alleles = sorted({a for gt in eligible["mother"] for a in gt})
        for p_allele in father_alleles:
            for m_allele in mother_alleles:
                if tuple(sorted((p_allele, m_allele))) == called_gt:
                    sources.append({
                        "child_genotype": fmt_gt(called_gt),
                        "copies": [
                            {"allele": p_allele, "from_role": "father", "from_sample": father_id},
                            {"allele": m_allele, "from_role": "mother", "from_sample": mother_id},
                        ],
                    })
    else:
        if kind == "y_nonpar":
            parent_role, parent_id = "father", father_id
        else:  # mt 与 X 非 PAR 男性子代均为母系传递
            parent_role, parent_id = "mother", mother_id
        parent_alleles = sorted({a for gt in eligible[parent_role] for a in gt})
        if called_gt[0] in parent_alleles:
            sources.append({
                "child_genotype": fmt_gt(called_gt),
                "copies": [{"allele": called_gt[0], "from_role": parent_role,
                            "from_sample": parent_id}],
            })

    if sources:
        return _result(site, base, CONSISTENT, participants, sources, None, [])

    # ---- 候选新生变异 ----
    def parent_support(allele):
        support = []
        for role, sid, _ploidy in parts:
            if role == "child":
                continue
            call = site.calls.get(sid)
            penalty = allele_support_penalty(call.pl if call else None, allele)
            support.append({
                "sample": sid,
                "role": role,
                "pl_penalty": penalty,
                "supported": penalty is not None and penalty <= th.max_pl_diff,
            })
        return support

    evidence = None
    if child_ploidy == 2 and len(set(called_gt)) == 2:
        c1, c2 = sorted(set(called_gt))
        all_parent = set(father_alleles) | set(mother_alleles)
        t1, t2 = c1 in all_parent, c2 in all_parent
        if t1 != t2:  # 恰好一个等位不可传递，其余等位可传递
            dn_allele = c2 if t1 else c1
            other = c1 if t1 else c2
            support = parent_support(dn_allele)
            if not any(s["supported"] for s in support):
                transmitters = [
                    role for role, alleles in
                    (("father", father_alleles), ("mother", mother_alleles))
                    if other in alleles
                ]
                evidence = {
                    "allele": dn_allele,
                    "child_genotype": fmt_gt(called_gt),
                    "transmitted_allele": {"allele": other, "from_role_options": transmitters},
                    "parent_support": support,
                }
    elif child_ploidy == 1:
        dn_allele = called_gt[0]
        support = parent_support(dn_allele)
        if not any(s["supported"] for s in support):
            evidence = {
                "allele": dn_allele,
                "child_genotype": fmt_gt(called_gt),
                "transmitted_allele": None,
                "parent_support": support,
            }

    if evidence is not None:
        return _result(site, base, CANDIDATE_DE_NOVO, participants, [], evidence, [])

    return _result(site, base, MENDELIAN_CONFLICT, participants, [], None,
                   ["证据充分但无可行传递组合，且不满足候选新生变异条件"])
