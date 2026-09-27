"""PED 家系文件解析。

PED 为空白分隔的文本，每行至少 6 列：

    family_id individual_id father_id mother_id sex phenotype

* ``father_id`` / ``mother_id`` 为 ``0`` 表示未知；
* ``sex`` 保留原始取值，由倍性规则中的 male/female_sex_value 解释；
* ``phenotype`` 列仅被解析，绝不用于任何判定或输出解释（本服务不做疾病诊断）。

仅当子代的双亲均在 PED 中登记时，才构成可审计的三人核心家系（trio）。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Person:
    family_id: str
    individual_id: str
    father_id: str | None
    mother_id: str | None
    sex: str
    phenotype: str


@dataclass(frozen=True)
class Trio:
    family_id: str
    child: Person
    father: Person
    mother: Person


def parse_ped(text: str) -> dict:
    """解析 PED 文本，返回 {(family_id, individual_id): Person}。"""
    persons: dict = {}
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 6:
            raise ValueError(f"PED 第 {lineno} 行列数不足 6：{line!r}")
        fam, iid, fid, mid, sex, pheno = cols[:6]
        key = (fam, iid)
        if key in persons:
            raise ValueError(f"PED 第 {lineno} 行个体重复：{fam}/{iid}")
        persons[key] = Person(
            family_id=fam,
            individual_id=iid,
            father_id=None if fid in ("0", "") else fid,
            mother_id=None if mid in ("0", "") else mid,
            sex=sex,
            phenotype=pheno,
        )
    if not persons:
        raise ValueError("PED 文件为空")
    return persons


def find_trios(persons: dict) -> tuple:
    """枚举完整三人核心家系，按 (family_id, child_id) 确定性排序。

    返回 (trios, notes)；notes 记录被跳过的不完整家系。
    """
    trios = []
    notes = []
    for person in persons.values():
        if not (person.father_id and person.mother_id):
            continue
        father = persons.get((person.family_id, person.father_id))
        mother = persons.get((person.family_id, person.mother_id))
        if father is None or mother is None:
            notes.append(
                f"个体 {person.family_id}/{person.individual_id} 的亲本未在 PED 中登记，跳过该家系"
            )
            continue
        trios.append(Trio(person.family_id, person, father, mother))
    trios.sort(key=lambda t: (t.family_id, t.child.individual_id))
    return trios, notes
