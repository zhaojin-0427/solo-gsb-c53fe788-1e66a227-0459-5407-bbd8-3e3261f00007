"""规范化模块单元测试。"""
import pysam
import pytest

from app.normalize import (
    left_align_trim,
    normalize_record,
    pl_expect,
    remap_pl,
)


def fetcher(seq):
    return lambda i: seq[i]


class TestLeftAlignTrim:
    def test_snv_unchanged(self):
        assert left_align_trim(fetcher("ACGT"), 1, "C", "T") == (1, "C", "T")

    def test_insertion_left_aligns_to_tract_start(self):
        seq = "AAATTTTTAAA"  # idx 3-7 为 T 串联
        pos0, ref, alt = left_align_trim(fetcher(seq), 5, "T", "TT")
        assert (pos0, ref, alt) == (2, "A", "AT")  # 锚定在串联重复之前的 A

    def test_deletion_left_aligns_to_tract_start(self):
        seq = "AAATTTTTAAA"
        pos0, ref, alt = left_align_trim(fetcher(seq), 5, "TT", "T")
        assert (pos0, ref, alt) == (2, "AT", "A")

    def test_common_suffix_trimmed(self):
        assert left_align_trim(fetcher("ACGT"), 1, "CTG", "CG") == (1, "CT", "C")

    def test_common_prefix_trimmed(self):
        assert left_align_trim(fetcher("ACGT"), 1, "CAA", "CA") == (1, "CA", "C")

    def test_unique_position_not_shifted(self):
        seq = "AATAAA"  # T 唯一，缺失无法左移
        pos0, ref, alt = left_align_trim(fetcher(seq), 1, "AT", "A")
        assert (pos0, ref, alt) == (1, "AT", "A")

    def test_contig_start_stops_shift(self):
        seq = "TTTT"  # 串联重复直达 contig 起点
        pos0, ref, alt = left_align_trim(fetcher(seq), 2, "T", "TT")
        assert pos0 == 0


class TestPlRemap:
    def test_pl_expect(self):
        assert pl_expect(3, 2) == 6
        assert pl_expect(3, 1) == 3
        assert pl_expect(2, 2) == 3

    def test_remap_diploid_third_allele(self):
        # 3 等位：基因型顺序 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
        pl = (0, 30, 300, 60, 120, 600)
        assert remap_pl(pl, alt_idx=2, ploidy=2) == (0, 60, 600)

    def test_remap_diploid_second_allele(self):
        pl = (0, 30, 300, 60, 120, 600)
        assert remap_pl(pl, alt_idx=1, ploidy=2) == (0, 30, 300)

    def test_remap_haploid(self):
        assert remap_pl((0, 50, 90), alt_idx=2, ploidy=1) == (0, 90)


VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=30>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def write_inputs(tmp_path, ref_seq, vcf_body):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + ref_seq + "\n")
    pysam.faidx(str(fa))
    vcf = tmp_path / "in.vcf"
    vcf.write_text(VCF_HEADER + vcf_body)
    return str(fa), str(vcf)


class TestNormalizeRecord:
    def test_multiallelic_split_and_remap(self, tmp_path):
        fa, vcf = write_inputs(
            tmp_path,
            "A" * 30,
            "chr1\t10\t.\tA\tG,T\t.\t.\t.\tGT:DP:GQ:PL\t"
            "0/1:30:99:30,0,300,90,90,900\t0/2:30:99:90,90,900,0,90,600\n",
        )
        with pysam.VariantFile(vcf) as vf, pysam.FastaFile(fa) as fasta:
            rec = next(iter(vf))
            sites, notes = normalize_record(rec, fasta, ["S1", "S2"])
        assert notes == []
        assert len(sites) == 2
        g_site, t_site = sites
        assert (g_site.ref, g_site.alt) == ("A", "G")
        # S1 0/1 -> 二等位 0/1；PL 取 (00, 01, 11)
        assert g_site.calls["S1"].gt == (0, 1)
        assert g_site.calls["S1"].pl == (30, 0, 300)
        # S2 0/2 在 A/G 位点含等位 2 -> 无法映射
        assert g_site.calls["S2"].gt is None
        # A/T 位点：S1 含等位 1 -> 无法映射；S2 0/2 -> (0,1)，PL 取 (00, 02, 22)
        assert t_site.calls["S1"].gt is None
        assert t_site.calls["S2"].gt == (0, 1)
        assert t_site.calls["S2"].pl == (90, 0, 600)

    def test_symbolic_alt_noted(self, tmp_path):
        fa, vcf = write_inputs(
            tmp_path,
            "A" * 30,
            "chr1\t10\t.\tA\t<DEL>\t.\t.\t.\tGT:DP:GQ:PL\t0/0:30:99:0,30,300\t0/0:30:99:0,30,300\n",
        )
        with pysam.VariantFile(vcf) as vf, pysam.FastaFile(fa) as fasta:
            rec = next(iter(vf))
            sites, notes = normalize_record(rec, fasta, ["S1"])
        assert sites == []
        assert len(notes) == 1 and "符号" in notes[0]["reason"]

    def test_ref_mismatch_excluded(self, tmp_path):
        fa, vcf = write_inputs(
            tmp_path,
            "A" * 30,
            "chr1\t10\t.\tC\tT\t.\t.\t.\tGT:DP:GQ:PL\t0/1:30:99:30,0,300\t0/0:30:99:0,30,300\n",
        )
        with pysam.VariantFile(vcf) as vf, pysam.FastaFile(fa) as fasta:
            rec = next(iter(vf))
            sites, notes = normalize_record(rec, fasta, ["S1"])
        assert sites == []
        assert len(notes) == 1 and "不一致" in notes[0]["reason"]

    def test_indel_normalized_against_reference(self, tmp_path):
        seq = "AA" + "T" * 8 + "AA" + "A" * 18  # idx 2-9 poly-T
        fa, vcf = write_inputs(
            tmp_path,
            seq,
            "chr1\t6\t.\tT\tTT\t.\t.\t.\tGT:DP:GQ:PL\t0/1:30:99:30,0,300\t0/0:30:99:0,30,300\n",
        )
        with pysam.VariantFile(vcf) as vf, pysam.FastaFile(fa) as fasta:
            rec = next(iter(vf))
            sites, notes = normalize_record(rec, fasta, ["S1"])
        assert notes == []
        assert (sites[0].pos, sites[0].ref, sites[0].alt) == (2, "A", "AT")
