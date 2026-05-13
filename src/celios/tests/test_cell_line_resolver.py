"""Tests for online-first cell-line identifier resolution.

Sanger API is primary with Cellosaurus fallback, validated with mocked HTTP
responses to keep tests deterministic.
"""

import pandas as pd

from celios.utils.cell_line_resolver import (
    CellLineResolver,
    detect_identifier_type,
    normalize_identifier,
    resolve_identifiers_to_sidm,
    resolve_model_ids_to_sidm,
    resolve_sidm_from_dataframe,
)


class FakeHttpClient:
    def __init__(self, response_map):
        self.response_map = response_map

    def get_json(self, url, params=None, timeout=10):
        key = (url, tuple(sorted((params or {}).items())))
        return self.response_map.get(key)


def test_normalize_identifier_core_patterns():
    assert normalize_identifier(" sidm 872 ") == "SIDM00872"
    assert normalize_identifier("ach_681") == "ACH-000681"
    assert normalize_identifier("RRID CVCL-0023") == "RRID:CVCL_0023"
    assert normalize_identifier("cvcl0023") == "CVCL_0023"


def test_detect_identifier_type():
    assert detect_identifier_type("SIDM00872") == "sidm"
    assert detect_identifier_type("ACH-000681") == "model_id"
    assert detect_identifier_type("RRID:CVCL_0023") == "rrid"
    assert detect_identifier_type("CVCL_0023") == "cvcl"
    assert detect_identifier_type("A549") == "name"


def test_sanger_primary_resolution_by_model_id():
    response_map = {
        (
            "https://api.cellmodelpassports.sanger.ac.uk/search",
            (("q", "ACH-000681"),),
        ): {"results": [{"SangerModelID": "SIDM00872"}]},
    }
    resolver = CellLineResolver(http_client=FakeHttpClient(response_map), use_cellosaurus_fallback=False)
    result = resolver.resolve_one("ACH-000681")
    assert result.status == "resolved"
    assert result.sidm == "SIDM00872"
    assert result.source == "sanger"


def test_cellosaurus_fallback_resolution():
    response_map = {
        # Sanger unresolved
        (
            "https://api.cellmodelpassports.sanger.ac.uk/search",
            (("q", "CVCL_0023"),),
        ): None,
        # Cellosaurus direct SIDM xref
        (
            "https://api.cellosaurus.org/cell-line/CVCL_0023",
            (("fields", "id,ac,sy,xref"), ("format", "json")),
        ): {"xref": ["Cell Model Passports: SIDM00872"]},
        (
            "https://api.cellosaurus.org/search/cell-line",
            (("fields", "id,ac,sy,xref"), ("format", "json"), ("q", "CVCL_0023")),
        ): None,
    }
    resolver = CellLineResolver(http_client=FakeHttpClient(response_map), use_cellosaurus_fallback=True)
    result = resolver.resolve_one("RRID:CVCL_0023")
    assert result.status == "resolved"
    assert result.sidm == "SIDM00872"
    assert result.source == "cellosaurus"


def test_cellosaurus_bridge_to_sanger_resolution():
    response_map = {
        # Sanger unresolved on name
        (
            "https://api.cellmodelpassports.sanger.ac.uk/search",
            (("q", "A549"),),
        ): None,
        # Cellosaurus returns ACH bridge identifier
        (
            "https://api.cellosaurus.org/search/cell-line",
            (("fields", "id,ac,sy,xref"), ("format", "json"), ("q", "A549")),
        ): {"hits": [{"xref": ["DepMap: ACH-000681"]}]},
        # Sanger resolves bridged ACH
        (
            "https://api.cellmodelpassports.sanger.ac.uk/search",
            (("q", "ACH-000681"),),
        ): {"results": [{"SangerModelID": "SIDM00872"}]},
    }
    resolver = CellLineResolver(http_client=FakeHttpClient(response_map), use_cellosaurus_fallback=True)
    result = resolver.resolve_one("A549")
    assert result.status == "resolved"
    assert result.sidm == "SIDM00872"
    assert result.source == "cellosaurus->sanger"


def test_resolve_identifiers_to_sidm_smoke(monkeypatch):
    def fake_resolve(self, identifier):
        if str(identifier).upper() == "A549":
            from celios.utils.cell_line_resolver import CellLineResolutionResult
            return CellLineResolutionResult(
                input_raw="A549",
                input_normalized="A549",
                detected_type="name",
                status="resolved",
                sidm="SIDM00872",
                matched_on="name",
                source="sanger",
            )
        from celios.utils.cell_line_resolver import CellLineResolutionResult
        return CellLineResolutionResult(
            input_raw=str(identifier),
            input_normalized=str(identifier),
            detected_type="name",
            status="unresolved",
            source="sanger",
        )

    monkeypatch.setattr(CellLineResolver, "resolve_one", fake_resolve)
    sidm_dict, not_found, _ = resolve_identifiers_to_sidm(["A549", "UNKNOWN"])
    assert sidm_dict["SIDM00872"] == "A549"
    assert "UNKNOWN" in not_found


def test_resolve_model_ids_to_sidm_smoke(monkeypatch):
    def fake_resolve(self, identifier):
        from celios.utils.cell_line_resolver import CellLineResolutionResult
        normalized = normalize_identifier(identifier)
        if normalized == "ACH-000681":
            return CellLineResolutionResult(
                input_raw=str(identifier),
                input_normalized=normalized,
                detected_type="model_id",
                status="resolved",
                sidm="SIDM00872",
                matched_on="model_id",
                source="sanger",
            )
        return CellLineResolutionResult(
            input_raw=str(identifier),
            input_normalized=normalized,
            detected_type="model_id",
            status="unresolved",
            source="sanger",
        )

    monkeypatch.setattr(CellLineResolver, "resolve_one", fake_resolve)
    model_to_sidm, not_found = resolve_model_ids_to_sidm(["ACH-000681", "ACH-999999"])
    assert model_to_sidm["ACH-000681"] == "SIDM00872"
    assert "ACH-999999" in not_found


def test_resolve_sidm_from_dataframe_identifier_priority(monkeypatch):
    df = pd.DataFrame(
        {
            "RRID": ["RRID:CVCL_0023", None],
            "ModelID": [None, "ACH-000001"],
            "cell_line_name": ["A549", "NIH:OVCAR-3"],
        }
    )

    def fake_resolve(self, identifier):
        from celios.utils.cell_line_resolver import CellLineResolutionResult
        normalized = normalize_identifier(identifier)
        lookup = {
            "RRID:CVCL_0023": "SIDM00872",
            "ACH-000001": "SIDM00105",
            "A549": "SIDM00872",
            "NIH:OVCAR-3": "SIDM00105",
        }
        sidm = lookup.get(normalized)
        if sidm:
            return CellLineResolutionResult(
                input_raw=str(identifier),
                input_normalized=normalized,
                detected_type=detect_identifier_type(normalized),
                status="resolved",
                sidm=sidm,
                matched_on="test",
                source="sanger",
            )
        return CellLineResolutionResult(
            input_raw=str(identifier),
            input_normalized=normalized,
            detected_type=detect_identifier_type(normalized),
            status="unresolved",
            source="sanger",
        )

    monkeypatch.setattr(CellLineResolver, "resolve_one", fake_resolve)
    sidm_dict, not_found, resolution_report = resolve_sidm_from_dataframe(df)

    assert sidm_dict["SIDM00872"] == "A549"
    assert sidm_dict["SIDM00105"] == "NIH:OVCAR-3"
    assert not not_found
    assert resolution_report["resolved"] == 2
