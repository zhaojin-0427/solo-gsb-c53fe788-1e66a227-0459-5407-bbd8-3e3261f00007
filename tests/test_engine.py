"""审计引擎单元测试。"""
from app.engine import (
    CANDIDATE_DE_NOVO,
    CONSISTENT,
    INAPPLICABLE,
    INSUFFICIENT,
    MENDELIAN_CONFLICT,
    Thresholds,
    allele_support_penalty,
    classify_site,
    eligible_genotypes,
    pl_derived_gq,
)
from app.normalize import NormalizedSite, SampleCall
from app.pedigree import Person, Trio
from app.ploidy import PloidyRuleSet

TH = Thresholds(min_dp=10, min_gq=20, max_pl_diff=20)
RULES = PloidyRuleSet.from_dict({
    "x_contigs": ["chrX"],
    "y_contigs": ["chrY"],
    "mt_contigs": ["chrM"],
    "par_regions": [
        {"contig": "chrX", "start": 0, "end": 50},
        {"contig": "chrY", "start": 0, "end": 50},
    ],
})


def make_trio(child_sex="2"):
    father = Person("FAM1", "S1", None, None, "1", "0")
    mother = Person("FAM1", "S2", None, None, "2", "0")
    child = Person("FAM1", "S3", "S1", "S2", child_sex, "0")
    return Trio("FAM1", child, father, mother)


def call(gt, dp=30, gq=99, pl=None):
    return SampleCall(gt=gt, dp=dp, gq=gq, pl=pl)


def site(calls, contig="chr1", pos0=99, ref="A", alt="G"):
    return NormalizedSite(contig, pos0, ref, alt, calls)


DIP0 = (0, 30, 300)      # 0/0 高置信
DIP0_NS = (0, 40, 400)   # 0/0 高置信，且对等位 1 罚分 40（> max_pl_diff）
DIP0_SUP = (0, 10, 300)  # 0/0 高置信，但对等位 1 罚分 10（<= max_pl_diff）
DIPHET = (30, 0, 300)    # 0/1 高置信
DIP1 = (300, 30, 0)      # 1/1 高置信
HAP0 = (0, 300)
HAP1 = (300, 0)


class TestEligibility:
    def test_low_dp_excluded(self):
        elig, reasons = eligible_genotypes(call((0, 1), dp=5), TH, 2)
        assert elig == [] and any("DP" in r for r in reasons)

    def test_low_gq_excluded(self):
        elig, reasons = eligible_genotypes(call((0, 1), gq=10), TH, 2)
        assert elig == [] and any("GQ" in r for r in reasons)

    def test_gq_derived_from_pl(self):
        assert pl_derived_gq((0, 45, 300)) == 45
        elig, _ = eligible_genotypes(call((0, 0), gq=None, pl=(0, 45, 300)), TH, 2)
        assert elig == [(0, 0)]

    def test_missing_gt(self):
        elig, reasons = eligible_genotypes(call(None), TH, 2)
        assert elig == [] and reasons

    def test_absent_sample(self):
        elig, reasons = eligible_genotypes(None, TH, 2)
        assert elig == [] and any("不在 VCF" in r for r in reasons)

    def test_ploidy_mismatch(self):
        elig, reasons = eligible_genotypes(call((0, 0)), TH, 1)
        assert elig == [] and any("倍性" in r for r in reasons)


class TestSupportPenalty:
    def test_diploid(self):
        assert allele_support_penalty((0, 40, 400), 1) == 40
        assert allele_support_penalty((0, 40, 400), 0) == 0
        assert allele_support_penalty((300, 30, 0), 0) == 30

    def test_haploid(self):
        assert allele_support_penalty((0, 300), 1) == 300

    def test_missing_pl(self):
        assert allele_support_penalty(None, 1) is None


class TestAutosome:
    def test_consistent_with_allelic_sources(self):
        s = site({"S1": call((0, 0), pl=DIP0), "S2": call((0, 1), pl=DIPHET),
                  "S3": call((0, 1), pl=DIPHET)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == CONSISTENT
        sources = r["detail"]["allelic_sources"]
        assert sources == [{
            "child_genotype": "0/1",
            "copies": [
                {"allele": 0, "from_role": "father", "from_sample": "S1"},
                {"allele": 1, "from_role": "mother", "from_sample": "S2"},
            ],
        }]

    def test_de_novo_candidate(self):
        s = site({"S1": call((0, 0), pl=DIP0_NS), "S2": call((0, 0), pl=DIP0_NS),
                  "S3": call((0, 1), pl=DIPHET)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == CANDIDATE_DE_NOVO
        ev = r["detail"]["de_novo_evidence"]
        assert ev["allele"] == 1
        assert ev["transmitted_allele"]["allele"] == 0
        assert all(not p["supported"] for p in ev["parent_support"])

    def test_parental_pl_support_blocks_de_novo(self):
        # 父亲 PL 显示对等位 1 有支持（罚分 10 <= 20）-> 孟德尔冲突而非新生变异
        s = site({"S1": call((0, 0), pl=DIP0_SUP), "S2": call((0, 0), pl=DIP0_NS),
                  "S3": call((0, 1), pl=DIPHET)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == MENDELIAN_CONFLICT

    def test_hom_alt_child_conflict(self):
        s = site({"S1": call((0, 0), pl=DIP0_NS), "S2": call((0, 0), pl=DIP0_NS),
                  "S3": call((1, 1), pl=DIP1)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == MENDELIAN_CONFLICT

    def test_het_parents_hom_child_consistent(self):
        s = site({"S1": call((0, 1), pl=DIPHET), "S2": call((0, 1), pl=DIPHET),
                  "S3": call((1, 1), pl=DIP1)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == CONSISTENT

    def test_insufficient_when_child_low_dp(self):
        s = site({"S1": call((0, 0), pl=DIP0), "S2": call((0, 1), pl=DIPHET),
                  "S3": call((0, 1), dp=2, pl=DIPHET)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == INSUFFICIENT
        assert any("S3" in reason for reason in r["detail"]["exclusion_reasons"])

    def test_insufficient_when_sample_missing(self):
        s = site({"S2": call((0, 1), pl=DIPHET), "S3": call((0, 1), pl=DIPHET)})
        r = classify_site(s, make_trio(), RULES, TH)
        assert r["verdict"] == INSUFFICIENT
        assert any("S1" in reason for reason in r["detail"]["exclusion_reasons"])


class TestSexChromosomesAndMt:
    def test_x_nonpar_daughter_consistent(self):
        s = site({"S1": call((1,), pl=HAP1), "S2": call((0, 1), pl=DIPHET),
                  "S3": call((0, 1), pl=DIPHET)}, contig="chrX", pos0=99)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["verdict"] == CONSISTENT
        copies = r["detail"]["allelic_sources"][0]["copies"]
        assert {c["from_role"]: c["allele"] for c in copies} == {"father": 1, "mother": 0}

    def test_x_nonpar_son_consistent(self):
        s = site({"S2": call((0, 1), pl=DIPHET), "S3": call((1,), pl=HAP1)},
                 contig="chrX", pos0=99)
        r = classify_site(s, make_trio("1"), RULES, TH)
        assert r["verdict"] == CONSISTENT
        copies = r["detail"]["allelic_sources"][0]["copies"]
        assert copies == [{"allele": 1, "from_role": "mother", "from_sample": "S2"}]

    def test_x_nonpar_son_de_novo(self):
        s = site({"S2": call((0, 0), pl=DIP0_NS), "S3": call((1,), pl=HAP1)},
                 contig="chrX", pos0=99)
        r = classify_site(s, make_trio("1"), RULES, TH)
        assert r["verdict"] == CANDIDATE_DE_NOVO

    def test_x_nonpar_father_diploid_is_ploidy_mismatch(self):
        s = site({"S1": call((1, 1), pl=DIP1), "S2": call((0, 1), pl=DIPHET),
                  "S3": call((0, 1), pl=DIPHET)}, contig="chrX", pos0=99)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["verdict"] == INSUFFICIENT
        assert any("倍性" in reason for reason in r["detail"]["exclusion_reasons"])

    def test_par_behaves_like_autosome(self):
        s = site({"S1": call((0, 1), pl=DIPHET), "S2": call((0, 0), pl=DIP0),
                  "S3": call((0, 1), pl=DIPHET)}, contig="chrX", pos0=9)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["region_class"] == "par"
        assert r["verdict"] == CONSISTENT

    def test_y_nonpar_son_consistent(self):
        s = site({"S1": call((1,), pl=HAP1), "S3": call((1,), pl=HAP1)},
                 contig="chrY", pos0=59)
        r = classify_site(s, make_trio("1"), RULES, TH)
        assert r["verdict"] == CONSISTENT

    def test_y_nonpar_daughter_inapplicable(self):
        s = site({"S1": call((1,), pl=HAP1), "S3": call(None)},
                 contig="chrY", pos0=59)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["verdict"] == INAPPLICABLE

    def test_mt_maternal_consistent(self):
        s = site({"S2": call((1,), pl=HAP1), "S3": call((1,), pl=HAP1)},
                 contig="chrM", pos0=29)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["verdict"] == CONSISTENT

    def test_mt_de_novo(self):
        s = site({"S2": call((0,), pl=HAP0), "S3": call((1,), pl=HAP1)},
                 contig="chrM", pos0=59)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["verdict"] == CANDIDATE_DE_NOVO

    def test_mt_conflict_when_mother_supports(self):
        s = site({"S2": call((0,), pl=(0, 5)), "S3": call((1,), pl=HAP1)},
                 contig="chrM", pos0=59)
        r = classify_site(s, make_trio("2"), RULES, TH)
        assert r["verdict"] == MENDELIAN_CONFLICT
