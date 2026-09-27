# 家系 VCF 遗传一致性审计 API

基于 Python + FastAPI + PostgreSQL + pysam 的家系孟德尔遗传一致性审计服务。
客户端先登记**不可变阈值版本**（最小 DP/GQ、最大 PL 差），再提交本地参考
FASTA、PED 家系、性染色体倍性规则与多样本 VCF；服务对规范化后的二等位位点
枚举亲子等位来源组合，给出逐位点的一致性判定。

> **声明**：本服务仅基于客户端提交的数据做孟德尔传递一致性审计。
> **不进行任何联网注释，不输出疾病诊断或任何临床建议。** PED 中的 phenotype
> 列仅被解析、绝不参与判定或输出解释。

---

## 1. 快速启动（Docker Compose 一键启动）

```bash
docker compose up --build
```

启动后：

| 入口 | 地址 |
| --- | --- |
| API 文档（Swagger UI） | http://localhost:8000/docs |
| OpenAPI Schema | http://localhost:8000/openapi.json |
| 健康检查 | http://localhost:8000/health |
| 服务根信息 | http://localhost:8000/ |

PostgreSQL 数据与上传文件保存在 Docker 卷（`pgdata`、`uploads`）中。

### 本地直接运行（可选，使用 SQLite）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

连接 PostgreSQL 请设置环境变量：

```bash
DATABASE_URL=postgresql+psycopg2://audit:audit@db:5432/pedigree_audit
UPLOAD_DIR=/data/uploads
```

---

## 2. 审计规则

### 2.1 阈值版本（不可变）

阈值版本一经登记即不可修改、删除（API 不提供变更接口），同一
`name + min_dp + min_gq + max_pl_diff` 组合去重。重放审计时引用其 id，
保证结果可按阈值版本精确复现。

| 参数 | 含义 |
| --- | --- |
| `min_dp` | 最小深度：样本 `DP >= min_dp` 才合格 |
| `min_gq` | 最小基因型质量：`GQ >= min_gq` 才合格；GQ 缺失时由 PL 推导（次小 PL − 最小 PL） |
| `max_pl_diff` | 最大 PL 差：亲本"支持"某等位的罚分须 `<= max_pl_diff` |

### 2.2 变异规范化

* 多等位记录按 ALT 拆分为独立**二等位位点**；
* indel 借助参考 FASTA 做**左对齐与最小化修剪**（与 `bcftools norm -f` 等价）；
* 各样本 GT/PL 重映射到二等位坐标系：GT 含其他交替等位的样本在该位点视为
  无法映射（信息不足）；PL 按 VCF 基因型顺序抽取；
* 符号等位（`<DEL>` 等）、`*`、REF 与参考不一致、contig 不在参考中的记录
  不参与审计，作为"排除理由"记入审计摘要 `exclusions`。

### 2.3 性染色体倍性规则（客户端提交 JSON）

`par_regions` 为 **0-based 半开区间** `[start, end)`；`male_sex_value` /
`female_sex_value` 对应 PED 第 5 列取值（默认 1/2）。

```json
{
  "x_contigs": ["chrX"],
  "y_contigs": ["chrY"],
  "mt_contigs": ["chrM"],
  "par_regions": [
    {"contig": "chrX", "start": 0, "end": 50},
    {"contig": "chrY", "start": 0, "end": 50}
  ],
  "male_sex_value": "1",
  "female_sex_value": "2"
}
```

### 2.4 传递规则

| 区域 | 规则 |
| --- | --- |
| `autosome` 常染色体 | 子代二倍体，两个等位分别可来自父亲、母亲 |
| `par` 伪常染色体区 | 按常染色体处理（父亲二倍体） |
| `x_nonpar` X 非 PAR | 男性单倍体、女性二倍体；儿子仅由母亲传递；女儿从父亲处获得其唯一 X 等位 |
| `y_nonpar` Y 非 PAR | 仅男性子代适用，父系单倍体传递；女性/未知性别子代判"不适用" |
| `mt` 线粒体 | 单倍体，仅母系传递，父亲不参与判定 |

亲本性别按家系角色推定（父亲=男、母亲=女）；子代性别取自 PED。

### 2.5 判定结论（verdict）

1. **`insufficient`（信息不足）**：任一必要样本不在 VCF、GT 缺失/无法映射、
   倍性不符、或 DP/GQ 不达标 —— 结果中列出具体样本与排除理由；
2. **`consistent`（一致）**：子代高置信基因型存在符合上述传递规则的等位来源
   组合 —— 结果中列出每一种等位来源（`allelic_sources`）；
3. **`candidate_de_novo`（候选新生变异）**：
   * 二倍体子代：杂合子的一个高置信等位**在双亲均无 PL 支持**（含该等位基因型
     的最小 PL − 亲本最小 PL > `max_pl_diff`；无 PL 时按检出基因型判定），
     且**另一个等位可由任一亲本传递**；
   * 单倍体子代（X 非 PAR 男性 / Y / mt）：子代唯一等位在参与传递的亲本处
     无 PL 支持；
4. **`mendelian_conflict`（孟德尔冲突）**：证据充分但无可行传递，且不满足
   候选新生变异条件（如双亲纯合 0/0、子代纯合 1/1；或亲本 PL 对该等位有支持）；
5. **`inapplicable`（不适用）**：区域/性别组合不适用。

### 2.6 输出结构与确定性排序

每个位点结果包含：参与样本（角色、期望倍性、DP/GQ、合格基因型、排除理由）、
等位来源（子代等位 ← 父/母样本）、新生变异证据（PL 罚分与可传递等位）、
排除理由。结果按

```
(contig 按 VCF 头顺序, pos, ref, alt, family_id, child_id)
```

**确定性排序**并赋予连续 `rank`；所有内部列表同样按固定顺序输出。同一阈值版本
对同一批输入重放，结果逐字节一致。

---

## 3. API 一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/threshold-versions` | 登记不可变阈值版本（重复组合 409） |
| GET | `/threshold-versions` | 列出阈值版本 |
| GET | `/threshold-versions/{id}` | 查询阈值版本 |
| POST | `/audits` | 提交审计（multipart：阈值版本 id、倍性规则 JSON、FASTA、PED、VCF） |
| GET | `/audits` | 列出审计任务 |
| GET | `/audits/{id}` | 审计摘要（各 verdict 计数、排除项、错误） |
| GET | `/audits/{id}/results` | 位点结果，支持 `verdict/child_id/family_id/contig` 过滤与分页 |
| POST | `/audits/{id}/replay` | 按另一个阈值版本重放（复用原文件与倍性规则） |

结果过滤参数：`limit`（默认 100，上限 1000）、`offset`，按 `rank` 升序返回。

---

## 4. 示例请求

示例数据位于 `examples/`：`reference.fa`（chr1/chrX/chrY/chrM 迷你参考）、
`family.ped`（父亲 S1、母亲 S2、女儿 S3）、`ploidy_rules.json`、`family.vcf`
（12 条记录，覆盖一致、候选新生变异、冲突、信息不足、多等位拆分、indel
左对齐、PAR、X/Y 非 PAR、线粒体等场景）。

一键跑完整示例（先 `docker compose up`，另开终端）：

```bash
bash examples/requests.sh
```

或手动执行：

```bash
# 1) 登记阈值版本
curl -X POST http://localhost:8000/threshold-versions \
  -H "Content-Type: application/json" \
  -d '{"name":"standard-v1","min_dp":10,"min_gq":20,"max_pl_diff":20}'

# 2) 提交审计
curl -X POST http://localhost:8000/audits \
  -F "threshold_version_id=1" \
  -F "ploidy_rules=<examples/ploidy_rules.json" \
  -F "reference_fasta=@examples/reference.fa" \
  -F "pedigree=@examples/family.ped" \
  -F "vcf=@examples/family.vcf"

# 3) 查看结果（按确定性排序）
curl "http://localhost:8000/audits/1/results?limit=100"

# 4) 仅看候选新生变异
curl "http://localhost:8000/audits/1/results?verdict=candidate_de_novo"

# 5) 按更严格阈值版本重放
curl -X POST http://localhost:8000/audits/1/replay \
  -H "Content-Type: application/json" -d '{"threshold_version_id":2}'
```

示例审计的预期计数：

| verdict | 数量 | 位点 |
| --- | --- | --- |
| consistent | 7 | chr1:100、chr1:180 G、chr1:194（由 chr1:200 左对齐）、chr1:220、chrX:10、chrX:100、chrM:30 |
| candidate_de_novo | 2 | chr1:120、chrM:60 |
| mendelian_conflict | 1 | chr1:140 |
| insufficient | 2 | chr1:160（子代低 DP）、chr1:180 T（多等位无法映射） |
| inapplicable | 1 | chrY:60（女性子代） |

---

## 5. 测试

```bash
pip install pytest httpx
python -m pytest -q
```

测试覆盖：indel 左对齐算法、多等位 GT/PL 重映射、符号等位与 REF 不一致排除、
各区域传递规则、四种核心判定、API 全流程、重放确定性与阈值敏感性。

---

## 6. 项目结构

```
app/
├── main.py        # FastAPI 路由（阈值版本、审计提交、结果查询、重放）
├── service.py     # 审计编排：解析输入、运行引擎、确定性排序入库
├── engine.py      # 核心判定：合格基因型、传递枚举、新生变异 PL 支持检验
├── normalize.py   # 多等位拆分、indel 左对齐修剪、GT/PL 重映射
├── ploidy.py      # 性染色体倍性规则与区域判定
├── pedigree.py    # PED 解析与核心家系枚举
├── models.py      # SQLAlchemy 模型（阈值版本/审计/位点结果）
├── schemas.py     # Pydantic 模式
├── database.py    # 数据库连接（PostgreSQL，测试用 SQLite）
└── config.py      # 环境变量配置
examples/          # 示例 FASTA/PED/VCF/倍性规则与请求脚本
tests/             # pytest 测试套件
```

## 7. 边界与限制

* 仅审计 PED 中双亲均已登记的三人核心家系；支持同一 VCF 中的多个家系
  （逐家系输出结果）；
* 支持单倍体/二倍体位点（常染色体、PAR、X/Y 非 PAR、线粒体）；更高倍性
  （如三倍体）判为信息不足；
* 变异位点的生物学注释、频率、临床意义不在服务范围内，也永远不会输出。
