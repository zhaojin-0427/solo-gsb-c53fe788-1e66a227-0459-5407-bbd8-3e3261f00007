"""VCF 多等位变异规范化。

对每条 VCF 记录：

1. 按 ALT 拆分为独立的二等位位点（多等位拆分）；
2. 借助客户端提供的参考 FASTA 对 indel 做左对齐与最小化修剪
   （与 ``bcftools norm -f`` 等价的确定性算法）；
3. 将各样本 GT / PL 重映射到拆分后的二等位坐标系：
   - GT 中等位 0 与当前 ALT 分别映射为 0/1；含有其他 ALT 等位的样本
     在该位点上视为无法映射（GT 置缺失，由审计引擎判为信息不足）；
   - PL 按 VCF 基因型排序规则抽取对应的二等位三元组/二元组。

符号等位（``<DEL>`` 等）、``*`` 等位、REF 与参考基因组不一致、contig 不在
参考基因组中的记录不入审计，而是以"排除理由"的形式记入审计摘要。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import comb

BASES = set("ACGTN")


@dataclass
class SampleCall:
    """重映射到二等位坐标系后的单样本基因型数据。"""

    gt: tuple | None  # 如 (0, 1)；缺失或无法映射为 None
    dp: int | None
    gq: object | None
    pl: tuple | None  # 二倍体 (PL00, PL01, PL11)；单倍体 (PL0, PL1)


@dataclass
class NormalizedSite:
    contig: str
    pos0: int  # 0-based
    ref: str
    alt: str
    calls: dict = field(default_factory=dict)  # sample_id -> SampleCall

    @property
    def pos(self) -> int:  # 1-based，用于输出
        return self.pos0 + 1


def left_align_trim(fetch_base, pos0: int, ref: str, alt: str):
    """对单个 REF/ALT 做左对齐与最小化修剪。

    ``fetch_base(i)`` 返回参考基因组 0-based 位置 ``i`` 处的碱基（大写）。
    返回 ``(pos0, ref, alt)``。算法（与 vt/bcftools 等价）：

    1. 修剪共同后缀（两等位均至少保留 1 个碱基）；
    2. 修剪共同前缀并右移位置；
    3. 若某一等位被修剪为空，则向左侧参考基因组扩展一个碱基；
    4. 若两等位末尾碱基相同（无法继续修剪），则整体向左旋转一个碱基；
    5. 重复直至无法再移动（到达 contig 起点时停止）。
    """
    if len(ref) == 1 and len(alt) == 1:
        return pos0, ref, alt
    while True:
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos0 += 1
        if len(ref) == 0 or len(alt) == 0:
            if pos0 == 0:
                break
            prev = fetch_base(pos0 - 1)
            ref, alt = prev + ref, prev + alt
            pos0 -= 1
            continue
        if ref[-1] == alt[-1]:
            if pos0 == 0:
                break
            prev = fetch_base(pos0 - 1)
            ref, alt = prev + ref[:-1], prev + alt[:-1]
            pos0 -= 1
            continue
        break
    return pos0, ref, alt


def pl_expect(n_alleles: int, ploidy: int) -> int:
    """VCF 中 PL 向量的预期长度：C(n_alleles + ploidy - 1, ploidy)。"""
    return comb(n_alleles + ploidy - 1, ploidy)


def remap_pl(pl: tuple, alt_idx: int, ploidy: int):
    """将多等位 PL 向量重映射为仅含 REF 与第 alt_idx 个 ALT 的二等位向量。

    VCF 二倍体基因型排序为 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, ...，
    基因型 (a, b)（a<=b）的索引为 b*(b+1)//2 + a；单倍体索引即等位编号。
    """
    if ploidy == 1:
        return (pl[0], pl[alt_idx])
    if ploidy == 2:
        i01 = alt_idx * (alt_idx + 1) // 2  # 基因型 (0, alt_idx)
        i11 = alt_idx * (alt_idx + 1) // 2 + alt_idx  # 基因型 (alt_idx, alt_idx)
        return (pl[0], pl[i01], pl[i11])
    return None  # 更高倍性不支持（审计时会因倍性不符被判为信息不足）


def _build_call(sample, alt_idx: int, n_alleles: int) -> SampleCall:
    keys = set(sample.keys())
    gt_raw = sample["GT"] if "GT" in keys else None
    gt = None
    ploidy = 0
    if gt_raw is not None and all(a is not None for a in gt_raw):
        mapped = []
        for allele in gt_raw:
            if allele == 0:
                mapped.append(0)
            elif allele == alt_idx:
                mapped.append(1)
            else:
                mapped = None
                break
        if mapped is not None:
            gt = tuple(mapped)
            ploidy = len(gt)
    pl = sample["PL"] if "PL" in keys else None
    pl_new = None
    if (
        pl is not None
        and ploidy in (1, 2)
        and len(pl) == pl_expect(n_alleles, ploidy)
    ):
        pl_new = remap_pl(tuple(pl), alt_idx, ploidy)
    dp = sample["DP"] if "DP" in keys else None
    gq = sample["GQ"] if "GQ" in keys else None
    return SampleCall(gt=gt, dp=dp, gq=gq, pl=pl_new)


def normalize_record(rec, fasta, sample_ids):
    """规范化一条 VCF 记录。

    返回 ``(sites, notes)``：sites 为 NormalizedSite 列表（每个可审计 ALT
    一个），notes 为排除理由字典列表（记入审计摘要）。
    """
    sites, notes = [], []
    contig = rec.contig
    pos0 = rec.pos - 1
    ref = (rec.ref or "").upper()
    alts = rec.alts or ()

    def note(reason, alt=None):
        notes.append(
            {
                "contig": contig,
                "pos": rec.pos,
                "ref": ref,
                "alt": alt,
                "reason": reason,
            }
        )

    if contig not in fasta.references:
        note("参考基因组中不存在该 contig，整条记录已排除")
        return sites, notes
    try:
        ref_ok = fasta.fetch(contig, pos0, pos0 + len(ref)).upper() == ref
    except ValueError:
        ref_ok = False
    if not ref_ok:
        note("REF 与客户端提供的参考基因组不一致，整条记录已排除")
        return sites, notes

    fetch_base = lambda i: fasta.fetch(contig, i, i + 1).upper()  # noqa: E731

    for alt_idx, alt in enumerate(alts, start=1):
        alt_u = (alt or "").upper()
        if alt_u == "*" or alt_u.startswith("<") or alt_u.endswith(">"):
            note(f"符号/特殊等位 {alt_u} 不参与审计，已跳过", alt=alt_u)
            continue
        if not set(alt_u) <= BASES or not set(ref) <= BASES:
            note("等位含非常规碱基符号，已跳过", alt=alt_u)
            continue
        npos0, nref, nalt = left_align_trim(fetch_base, pos0, ref, alt_u)
        calls = {}
        for sid in sample_ids:
            if sid in rec.samples:
                calls[sid] = _build_call(rec.samples[sid], alt_idx, len(alts) + 1)
        sites.append(NormalizedSite(contig, npos0, nref, nalt, calls))
    return sites, notes
