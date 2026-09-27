#!/usr/bin/env bash
# 示例请求：登记阈值版本 -> 提交审计 -> 查询结果 -> 按新阈值版本重放
set -euo pipefail

BASE="${BASE:-http://localhost:8000}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "== 1. 登记不可变阈值版本（最小 DP/GQ、最大 PL 差） =="
TV=$(curl -sS -X POST "$BASE/threshold-versions" \
  -H "Content-Type: application/json" \
  -d '{"name": "standard-v1", "min_dp": 10, "min_gq": 20, "max_pl_diff": 20,
       "description": "标准审计阈值"}')
echo "$TV"
TV_ID=$(echo "$TV" | python3 -c "import sys, json; print(json.load(sys.stdin)['id'])")

echo "== 2. 提交审计（参考 FASTA + PED + 倍性规则 + 多样本 VCF） =="
AUDIT=$(curl -sS -X POST "$BASE/audits" \
  -F "threshold_version_id=$TV_ID" \
  -F "ploidy_rules=<$HERE/ploidy_rules.json" \
  -F "reference_fasta=@$HERE/reference.fa" \
  -F "pedigree=@$HERE/family.ped" \
  -F "vcf=@$HERE/family.vcf")
echo "$AUDIT"
AUDIT_ID=$(echo "$AUDIT" | python3 -c "import sys, json; print(json.load(sys.stdin)['id'])")

echo "== 3. 查询全部结果（按确定性排序） =="
curl -sS "$BASE/audits/$AUDIT_ID/results?limit=100" | python3 -m json.tool

echo "== 4. 仅看候选新生变异 =="
curl -sS "$BASE/audits/$AUDIT_ID/results?verdict=candidate_de_novo" | python3 -m json.tool

echo "== 5. 登记更严格的阈值版本并重放 =="
TV2=$(curl -sS -X POST "$BASE/threshold-versions" \
  -H "Content-Type: application/json" \
  -d '{"name": "strict-v1", "min_dp": 50, "min_gq": 30, "max_pl_diff": 10}')
TV2_ID=$(echo "$TV2" | python3 -c "import sys, json; print(json.load(sys.stdin)['id'])")
curl -sS -X POST "$BASE/audits/$AUDIT_ID/replay" \
  -H "Content-Type: application/json" \
  -d "{\"threshold_version_id\": $TV2_ID}" | python3 -m json.tool
