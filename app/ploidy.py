"""性染色体倍性规则与区域判定。

规则由客户端在提交审计时以 JSON 提供：

.. code-block:: json

    {
      "x_contigs": ["chrX", "X"],
      "y_contigs": ["chrY", "Y"],
      "mt_contigs": ["chrM", "MT"],
      "par_regions": [
        {"contig": "chrX", "start": 10000, "end": 2781479},
        {"contig": "chrY", "start": 10000, "end": 2781479}
      ],
      "male_sex_value": "1",
      "female_sex_value": "2"
    }

* ``par_regions`` 为伪常染色体区域（PAR），坐标为 **0-based 半开区间** ``[start, end)``。
* ``male_sex_value`` / ``female_sex_value`` 为 PED 第 5 列中表示性别的取值，默认 ``1``/``2``。

区域类别（region_class）：

* ``autosome``   常染色体：二倍体，双亲各传一个等位；
* ``par``        伪常染色体区域：按常染色体处理；
* ``x_nonpar``   X 非 PAR 区：男性单倍体、女性二倍体；儿子仅母系传递，女儿父传其唯一 X、母传其一；
* ``y_nonpar``   Y 非 PAR 区：仅男性子代适用，父系单倍体传递；
* ``mt``         线粒体：单倍体，仅母系传递（父亲不参与判定）。
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ParRegion:
    contig: str
    start: int  # 0-based，含
    end: int  # 0-based，不含


@dataclass(frozen=True)
class PloidyRuleSet:
    x_contigs: frozenset = field(default_factory=frozenset)
    y_contigs: frozenset = field(default_factory=frozenset)
    mt_contigs: frozenset = field(default_factory=frozenset)
    par_regions: tuple = ()
    male_sex_value: str = "1"
    female_sex_value: str = "2"

    @classmethod
    def from_dict(cls, data) -> "PloidyRuleSet":
        if not isinstance(data, dict):
            raise ValueError("倍性规则须为 JSON 对象")

        def names(key):
            value = data.get(key, [])
            if not isinstance(value, list):
                raise ValueError(f"倍性规则字段 {key} 须为数组")
            return frozenset(str(x) for x in value)

        pars = []
        for item in data.get("par_regions", []):
            if not isinstance(item, dict):
                raise ValueError("par_regions 元素须为对象")
            try:
                contig = str(item["contig"])
                start = int(item["start"])
                end = int(item["end"])
            except (KeyError, TypeError, ValueError):
                raise ValueError("par_regions 元素须含 contig/start/end 且坐标为整数")
            if start < 0 or end <= start:
                raise ValueError(f"PAR 区间非法：[{start}, {end})")
            pars.append(ParRegion(contig, start, end))
        return cls(
            x_contigs=names("x_contigs"),
            y_contigs=names("y_contigs"),
            mt_contigs=names("mt_contigs"),
            par_regions=tuple(pars),
            male_sex_value=str(data.get("male_sex_value", "1")),
            female_sex_value=str(data.get("female_sex_value", "2")),
        )

    def sex_of(self, raw: str) -> str:
        """PED 性别原始值 -> male / female / unknown。"""
        if raw == self.male_sex_value:
            return "male"
        if raw == self.female_sex_value:
            return "female"
        return "unknown"

    def region_kind(self, contig: str, pos0: int) -> str:
        """按 contig 与 0-based 位置判定区域类别。"""
        if contig in self.mt_contigs:
            return "mt"
        for par in self.par_regions:
            if par.contig == contig and par.start <= pos0 < par.end:
                return "par"
        if contig in self.x_contigs:
            return "x_nonpar"
        if contig in self.y_contigs:
            return "y_nonpar"
        return "autosome"
