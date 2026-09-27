"""API 集成测试：使用 examples/ 示例数据走完整审计流程。"""
import json
from pathlib import Path

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

EXPECTED_COUNTS = {
    "consistent": 7,
    "candidate_de_novo": 2,
    "mendelian_conflict": 1,
    "insufficient": 2,
    "inapplicable": 1,
}


def create_threshold(client, name="std", min_dp=10, min_gq=20, max_pl_diff=20):
    resp = client.post("/threshold-versions", json={
        "name": name,
        "min_dp": min_dp,
        "min_gq": min_gq,
        "max_pl_diff": max_pl_diff,
        "description": "测试阈值",
    })
    assert resp.status_code == 201, resp.text
    return resp.json()["id"]


def run_example_audit(client, tv_id):
    rules = (EXAMPLES / "ploidy_rules.json").read_text()
    with open(EXAMPLES / "reference.fa", "rb") as fa, \
            open(EXAMPLES / "family.ped", "rb") as ped, \
            open(EXAMPLES / "family.vcf", "rb") as vcf:
        resp = client.post(
            "/audits",
            data={"threshold_version_id": str(tv_id), "ploidy_rules": rules},
            files={
                "reference_fasta": ("reference.fa", fa, "text/plain"),
                "pedigree": ("family.ped", ped, "text/plain"),
                "vcf": ("family.vcf", vcf, "text/plain"),
            },
        )
    assert resp.status_code == 201, resp.text
    body = resp.json()
    assert body["status"] == "completed", body["error"]
    return body


def get_all_results(client, audit_id):
    resp = client.get(f"/audits/{audit_id}/results", params={"limit": 1000})
    assert resp.status_code == 200
    return resp.json()


def test_threshold_version_immutable_and_unique(client):
    tv_id = create_threshold(client)
    # 相同名称+参数组合 -> 409
    dup = client.post("/threshold-versions", json={
        "name": "std", "min_dp": 10, "min_gq": 20, "max_pl_diff": 20,
    })
    assert dup.status_code == 409
    # 不提供修改/删除接口
    assert client.put(f"/threshold-versions/{tv_id}", json={}).status_code == 405
    assert client.delete(f"/threshold-versions/{tv_id}").status_code == 405
    got = client.get(f"/threshold-versions/{tv_id}")
    assert got.status_code == 200 and got.json()["min_dp"] == 10


def test_audit_example_data(client):
    tv_id = create_threshold(client)
    audit = run_example_audit(client, tv_id)
    assert audit["counts"] == EXPECTED_COUNTS
    assert audit["trio_count"] == 1
    assert audit["sample_names"] == ["S1", "S2", "S3"]

    page = get_all_results(client, audit["id"])
    assert page["total"] == sum(EXPECTED_COUNTS.values())
    items = page["items"]
    # 确定性排序：rank 连续递增
    assert [it["rank"] for it in items] == list(range(1, len(items) + 1))

    by_pos = {(it["contig"], it["pos"]): it for it in items}

    # 一致：等位来源指向母亲
    it = by_pos[("chr1", 100)]
    assert it["verdict"] == "consistent"
    copies = it["detail"]["allelic_sources"][0]["copies"]
    assert {c["from_role"]: c["allele"] for c in copies} == {"father": 0, "mother": 1}
    assert [p["sample"] for p in it["detail"]["participants"]] == ["S1", "S2", "S3"]

    # 候选新生变异：chr1:120 与 chrM:60
    it = by_pos[("chr1", 120)]
    assert it["verdict"] == "candidate_de_novo"
    ev = it["detail"]["de_novo_evidence"]
    assert ev["allele"] == 1 and ev["transmitted_allele"]["allele"] == 0
    assert len(ev["parent_support"]) == 2
    it = by_pos[("chrM", 60)]
    assert it["verdict"] == "candidate_de_novo"
    assert it["detail"]["region_class"] == "mt"

    # 孟德尔冲突：chr1:140（子代纯合 ALT）
    assert by_pos[("chr1", 140)]["verdict"] == "mendelian_conflict"

    # 信息不足：chr1:160（子代 DP=5）
    it = by_pos[("chr1", 160)]
    assert it["verdict"] == "insufficient"
    assert any("S3" in r for r in it["detail"]["exclusion_reasons"])

    # 多等位拆分：chr1:180 产生两个位点
    split = [it for it in items if it["contig"] == "chr1" and it["pos"] == 180]
    assert len(split) == 2
    verdicts = {it["alt"]: it["verdict"] for it in split}
    assert verdicts == {"G": "consistent", "T": "insufficient"}

    # indel 左对齐：chr1:200 T>TT 归一化到 chr1:194 A>AT
    it = by_pos[("chr1", 194)]
    assert (it["ref"], it["alt"], it["verdict"]) == ("A", "AT", "consistent")

    # X 非 PAR 区（女儿）：父传唯一 X 等位
    it = by_pos[("chrX", 100)]
    assert it["verdict"] == "consistent"
    assert it["detail"]["region_class"] == "x_nonpar"
    copies = it["detail"]["allelic_sources"][0]["copies"]
    assert {c["from_role"]: c["allele"] for c in copies} == {"father": 1, "mother": 0}

    # PAR 区按常染色体处理
    assert by_pos[("chrX", 10)]["detail"]["region_class"] == "par"

    # Y 非 PAR 区对女性子代不适用
    it = by_pos[("chrY", 60)]
    assert it["verdict"] == "inapplicable"
    assert it["detail"]["exclusion_reasons"]

    # 过滤查询
    resp = client.get(f"/audits/{audit['id']}/results",
                      params={"verdict": "candidate_de_novo"})
    assert resp.json()["total"] == 2
    resp = client.get(f"/audits/{audit['id']}/results", params={"contig": "chrM"})
    assert resp.json()["total"] == 2


def test_replay_deterministic_and_threshold_sensitive(client):
    tv_id = create_threshold(client)
    audit = run_example_audit(client, tv_id)

    # 相同阈值重放 -> 结果完全一致
    same = client.post(f"/audits/{audit['id']}/replay",
                       json={"threshold_version_id": tv_id})
    assert same.status_code == 201
    same = same.json()
    assert same["replay_of_audit_id"] == audit["id"]
    assert same["counts"] == audit["counts"]
    original_items = get_all_results(client, audit["id"])["items"]
    replay_items = get_all_results(client, same["id"])["items"]
    strip = lambda items: [
        {k: it[k] for k in ("rank", "contig", "pos", "ref", "alt", "verdict", "detail")}
        for it in items
    ]
    assert strip(original_items) == strip(replay_items)

    # 更严阈值（min_dp=50）重放 -> 全部位点信息不足/不适用
    strict_id = create_threshold(client, name="strict", min_dp=50)
    strict = client.post(f"/audits/{audit['id']}/replay",
                         json={"threshold_version_id": strict_id}).json()
    assert strict["counts"] == {"insufficient": 12, "inapplicable": 1}


def test_invalid_inputs(client):
    tv_id = create_threshold(client)
    # 倍性规则非法
    with open(EXAMPLES / "reference.fa", "rb") as fa, \
            open(EXAMPLES / "family.ped", "rb") as ped, \
            open(EXAMPLES / "family.vcf", "rb") as vcf:
        resp = client.post(
            "/audits",
            data={"threshold_version_id": str(tv_id), "ploidy_rules": "{not json"},
            files={
                "reference_fasta": ("reference.fa", fa, "text/plain"),
                "pedigree": ("family.ped", ped, "text/plain"),
                "vcf": ("family.vcf", vcf, "text/plain"),
            },
        )
    assert resp.status_code == 422
    # 阈值版本不存在
    with open(EXAMPLES / "reference.fa", "rb") as fa, \
            open(EXAMPLES / "family.ped", "rb") as ped, \
            open(EXAMPLES / "family.vcf", "rb") as vcf:
        resp = client.post(
            "/audits",
            data={"threshold_version_id": "999", "ploidy_rules": "{}"},
            files={
                "reference_fasta": ("reference.fa", fa, "text/plain"),
                "pedigree": ("family.ped", ped, "text/plain"),
                "vcf": ("family.vcf", vcf, "text/plain"),
            },
        )
    assert resp.status_code == 404
    # 重放到不存在的阈值版本
    assert client.post("/audits/1/replay",
                       json={"threshold_version_id": 999}).status_code == 404
