from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "data" / "manifest.json"


def _manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def test_canonical_workbook_matches_manifest_checksum() -> None:
    manifest = _manifest()
    source = ROOT / manifest["canonical_input"]["path"]

    assert source.is_file()
    assert source.stat().st_size == manifest["canonical_input"]["size_bytes"]
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    assert digest == manifest["canonical_input"]["sha256"]


def test_workbook_schema_and_declared_keys_match_manifest() -> None:
    manifest = _manifest()
    source = ROOT / manifest["canonical_input"]["path"]
    expected_names = [sheet["name"] for sheet in manifest["sheets"]]
    excel = pd.ExcelFile(source)

    assert excel.sheet_names == expected_names
    for contract in manifest["sheets"]:
        frame = excel.parse(contract["name"])
        assert len(frame) == contract["rows"]
        assert frame.columns.tolist() == contract["columns"]
        assert not frame.isna().any().any()
        assert not frame.duplicated().any()
        assert not frame.duplicated(subset=contract["primary_key"]).any()


def test_cohorts_and_qpcr_pool_membership_are_internally_consistent() -> None:
    manifest = _manifest()
    source = ROOT / manifest["canonical_input"]["path"]
    outcomes = pd.read_excel(source, sheet_name="outcomes_6dpf")
    cfos = pd.read_excel(source, sheet_name="cfos_cohort_features")
    pools = pd.read_excel(source, sheet_name="cfos_pools")
    ptz = pd.read_excel(source, sheet_name="ptz_challenge")

    followed_ids = set(outcomes["fish_id"].astype(str))
    cfos_ids = set(cfos["fish_id"].astype(str))
    ptz_ids = set(ptz["fish_id"].astype(str))
    assert len(followed_ids) == 133
    assert len(cfos_ids) == 86
    assert len(ptz_ids) == 34
    assert not (followed_ids & cfos_ids or followed_ids & ptz_ids or cfos_ids & ptz_ids)
    assert len(followed_ids | cfos_ids | ptz_ids) == 253

    pooled_ids: list[str] = []
    for row in pools.itertuples(index=False):
        members = str(row.pooled_fish_ids).split(";")
        assert len(members) == int(row.n_larvae_in_pool) == 4
        assert len(set(members)) == 4
        pooled_ids.extend(members)

    assert len(pooled_ids) == 72
    assert len(set(pooled_ids)) == 72
    assert set(pooled_ids) <= cfos_ids
    assert len(cfos_ids - set(pooled_ids)) == 14
