"""Cell line identifier normalization and SIDM resolution.

Online-first implementation:
- Normalize common identifier formats (SIDM, ACH, RRID/CVCL)
- Detect identifier type
- Resolve identifiers primarily via Sanger Cell Model Passports API
- Optional fallback via Cellosaurus API

This resolver is intentionally schema-tolerant to handle minor API response
shape differences across versions.
"""

import re
import unicodedata
from dataclasses import dataclass
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import json
from typing import Any, Dict, List, Optional, Tuple
import hashlib
import os
from pathlib import Path

import pandas as pd


SIDM_PATTERN = re.compile(r"^SIDM\s*[-_]?\s*(\d+)$", re.IGNORECASE)
ACH_PATTERN = re.compile(r"^ACH\s*[-_]?\s*(\d+)$", re.IGNORECASE)
CVCL_PATTERN = re.compile(r"^CVCL\s*[-_]?\s*([A-Z0-9]+)$", re.IGNORECASE)
SIDM_SCAN_PATTERN = re.compile(r"SIDM\d{5,}", re.IGNORECASE)
ACH_SCAN_PATTERN = re.compile(r"ACH-\d{6}", re.IGNORECASE)
CVCL_SCAN_PATTERN = re.compile(r"CVCL_[A-Z0-9]+", re.IGNORECASE)

_RECORD_LIST_KEYS = (
    "results",
    "hits",
    "items",
    "records",
    "data",
    "models",
    "cell_lines",
    "cellLines",
    "entries",
)

_SIDM_FIELD_KEYS = {"sidm", "sangermodelid", "sanger_model_id", "sangerid"}
_IDENTIFIER_FIELD_KEYS = {"identifier", "identifiers"}
_TRUSTED_SIDM_PATHS = {
    ("sidm",),
    ("sangermodelid",),
    ("model", "sidm"),
    ("model", "sangermodelid"),
}
_ACH_FIELD_KEYS = {"ach", "modelid", "model_id", "sangermodelid", "depmapid", "depmap_id", "id", "ac", "accession"}
_RRID_FIELD_KEYS = {"rrid", "id", "ac", "accession"}
_CVCL_FIELD_KEYS = {"cvcl", "rrid", "id", "ac", "accession"}
_NAME_FIELD_KEYS = {
    "name",
    "names",
    "synonym",
    "synonyms",
    "alias",
    "aliases",
    "cell_line_name",
    "celllinename",
    "display_name",
    "displayname",
    "label",
    "title",
    "symbol",
    "primary_name",
    "preferred_name",
}


def _collapse_whitespace(value: str) -> str:
    return re.sub(r"\s+", " ", value).strip()


def normalize_identifier(value: object) -> str:
    """Normalize an arbitrary cell-line identifier string.

    Rules:
    - NFKC normalize + trim + collapse whitespace
    - Uppercase
    - Canonical prefixes for SIDM/ACH/CVCL/RRID
    """
    if value is None:
        return ""

    text = unicodedata.normalize("NFKC", str(value))
    text = _collapse_whitespace(text)
    if not text:
        return ""

    upper = text.upper()

    sidm_match = SIDM_PATTERN.match(upper)
    if sidm_match:
        digits = sidm_match.group(1)
        if len(digits) < 5:
            digits = digits.zfill(5)
        return f"SIDM{digits}"

    ach_match = ACH_PATTERN.match(upper)
    if ach_match:
        digits = ach_match.group(1).zfill(6)
        return f"ACH-{digits}"

    cvcl = extract_cvcl(upper)
    if cvcl:
        if upper.startswith("RRID"):
            return f"RRID:{cvcl}"
        if upper.startswith("CVCL"):
            return cvcl
        if "RRID" in upper:
            return f"RRID:{cvcl}"

    # Default: uppercase + collapsed whitespace
    return upper


def extract_cvcl(value: object) -> Optional[str]:
    """Extract and normalize CVCL token from arbitrary text."""
    if value is None:
        return None

    text = unicodedata.normalize("NFKC", str(value)).upper()
    match = re.search(r"CVCL\s*[-_]?\s*([A-Z0-9]+)", text)
    if not match:
        return None
    return f"CVCL_{match.group(1)}"


def detect_identifier_type(value: object) -> str:
    """Detect identifier type after normalization.

    Returns one of: sidm, model_id, rrid, cvcl, name
    """
    normalized = normalize_identifier(value)
    if not normalized:
        return "name"
    if normalized.startswith("SIDM"):
        return "sidm"
    if normalized.startswith("ACH-"):
        return "model_id"
    if normalized.startswith("RRID:CVCL_"):
        return "rrid"
    if normalized.startswith("CVCL_"):
        return "cvcl"
    raw_text = unicodedata.normalize("NFKC", str(value)).strip().upper() if value is not None else ""
    if raw_text.startswith("RRID") and extract_cvcl(raw_text):
        return "rrid"
    if raw_text.startswith("CVCL") and extract_cvcl(raw_text):
        return "cvcl"
    return "name"


def normalize_name_key(value: object) -> str:
    """Build a comparison key for name/synonym matching."""
    if value is None:
        return ""
    text = unicodedata.normalize("NFKC", str(value)).upper().strip()
    return re.sub(r"[^A-Z0-9]", "", text)


def _normalize_path_tokens(path: Tuple[str, ...]) -> List[str]:
    tokens: List[str] = []
    for part in path:
        cleaned = re.sub(r"[^a-z0-9]+", "", str(part).lower())
        if cleaned:
            tokens.append(cleaned)
    return tokens


def _path_matches(path: Tuple[str, ...], allowed_keys: set) -> bool:
    tokens = _normalize_path_tokens(path)
    return any(
        token in allowed_keys or any(token.startswith(allowed) for allowed in allowed_keys)
        for token in tokens
    )


def _walk_record_fields(value: Any, path: Tuple[str, ...] = ()):
    if isinstance(value, dict):
        for key, sub_value in value.items():
            yield from _walk_record_fields(sub_value, path + (str(key),))
    elif isinstance(value, list):
        for item in value:
            yield from _walk_record_fields(item, path)
    else:
        yield path, "" if value is None else str(value)


def _payload_to_records(payload: Any) -> List[Dict[str, Any]]:
    if payload is None:
        return []
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        for key in _RECORD_LIST_KEYS:
            value = payload.get(key)
            if isinstance(value, list):
                records = [item for item in value if isinstance(item, dict)]
                if records:
                    return records
            if isinstance(value, dict):
                nested_records = _payload_to_records(value)
                if nested_records:
                    return nested_records
        return [payload]
    return []


def _jsonapi_payload_records(payload: Any) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    data_records: List[Dict[str, Any]] = []
    included_records: List[Dict[str, Any]] = []
    if isinstance(payload, dict):
        data = payload.get("data")
        if isinstance(data, list):
            data_records = [item for item in data if isinstance(item, dict)]
        elif isinstance(data, dict):
            data_records = [data]
        included = payload.get("included")
        if isinstance(included, list):
            included_records = [item for item in included if isinstance(item, dict)]
    return data_records, included_records


def _record_relationship_identifier_ids(record: Dict[str, Any]) -> List[str]:
    relationships = record.get("relationships")
    if not isinstance(relationships, dict):
        return []
    identifiers = relationships.get("identifiers")
    if not isinstance(identifiers, dict):
        return []
    data = identifiers.get("data")
    if not isinstance(data, list):
        return []
    ids: List[str] = []
    for item in data:
        if isinstance(item, dict) and item.get("id") is not None:
            item_id = str(item.get("id"))
            if item_id not in ids:
                ids.append(item_id)
    return ids


def _record_identifier_value(record: Dict[str, Any]) -> Optional[str]:
    attributes = record.get("attributes")
    if isinstance(attributes, dict):
        for key in ("identifier", "accession", "value"):
            value = attributes.get(key)
            if value is not None and str(value).strip():
                return str(value).strip()
    for key in ("identifier", "accession", "value"):
        value = record.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _unique_values(values: List[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def _dedupe_records(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen = set()
    ordered: List[Dict[str, Any]] = []
    for record in records:
        key = json.dumps(record, sort_keys=True, default=str)
        if key not in seen:
            seen.add(key)
            ordered.append(record)
    return ordered


def _resolve_records_from_payloads(payloads: List[Any], normalized: str, identifier_type: str):
    raw_hit_count = 0
    matched_records: List[Dict[str, Any]] = []
    for payload in payloads:
        records = _payload_to_records(payload)
        raw_hit_count += len(records)
        for record in records:
            if _record_matches_identifier(record, normalized, identifier_type):
                matched_records.append(record)
    matched_records = _dedupe_records(matched_records)
    sidm_candidates = _unique_values(
        [sidm for record in matched_records for sidm in _extract_sidms_from_record(record)]
    )
    return raw_hit_count, matched_records, sidm_candidates


def _record_matches_identifier(record: Dict[str, Any], normalized: str, identifier_type: str) -> bool:
    if not record:
        return False

    name_key = normalize_name_key(normalized)
    cvcl_token = extract_cvcl(normalized)

    for path, value in _walk_record_fields(record):
        if not value:
            continue
        if identifier_type == "sidm":
            if _path_matches(path, _SIDM_FIELD_KEYS):
                if normalize_identifier(value) == normalized:
                    return True
                sidm_match = SIDM_SCAN_PATTERN.search(value)
                if sidm_match and normalize_identifier(sidm_match.group(0)) == normalized:
                    return True
        elif identifier_type == "model_id":
            if _path_matches(path, _ACH_FIELD_KEYS) and normalize_identifier(value) == normalized:
                return True
        elif identifier_type in {"cvcl", "rrid"}:
            if _path_matches(path, _CVCL_FIELD_KEYS | _RRID_FIELD_KEYS | _IDENTIFIER_FIELD_KEYS) or tuple(_normalize_path_tokens(path)) == ("xref",):
                normalized_value = normalize_identifier(value)
                if normalized_value == normalized:
                    return True
                if cvcl_token and extract_cvcl(value) == cvcl_token:
                    return True
        else:
            if _path_matches(path, _NAME_FIELD_KEYS) and normalize_name_key(value) == name_key:
                return True

    return False


def _extract_sidms_from_record(record: Dict[str, Any], allow_xref_sidm: bool = False) -> List[str]:
    sidms: List[str] = []
    record_id = record.get("id")
    record_type = record.get("type")
    if record_type == "model" and record_id is not None:
        normalized_id = normalize_identifier(record_id)
        if detect_identifier_type(normalized_id) == "sidm" and normalized_id not in sidms:
            sidms.append(normalized_id)
    for path, value in _walk_record_fields(record):
        if not value:
            continue
        path_norm = tuple(_normalize_path_tokens(path))
        if allow_xref_sidm and path_norm == ("xref",):
            sidm_match = SIDM_SCAN_PATTERN.search(value)
            if sidm_match:
                normalized = normalize_identifier(sidm_match.group(0))
                if detect_identifier_type(normalized) == "sidm" and normalized not in sidms:
                    sidms.append(normalized)
        elif path_norm in _TRUSTED_SIDM_PATHS:
            normalized = normalize_identifier(value)
            if detect_identifier_type(normalized) == "sidm" and normalized not in sidms:
                sidms.append(normalized)
    return sidms


def _extract_trusted_sidm_sources(record: Dict[str, Any], allow_xref_sidm: bool = False) -> List[Dict[str, str]]:
    sources: List[Dict[str, str]] = []
    seen = set()
    record_id = record.get("id")
    record_type = record.get("type")
    if record_type == "model" and record_id is not None:
        normalized_id = normalize_identifier(record_id)
        if detect_identifier_type(normalized_id) == "sidm":
            key = (normalized_id, "id")
            sources.append({"sidm": normalized_id, "source_field": "id", "raw_value": str(record_id)})
            seen.add(key)
    for path, value in _walk_record_fields(record):
        if not value:
            continue
        path_norm = tuple(_normalize_path_tokens(path))
        source_field = ".".join(path)
        if allow_xref_sidm and path_norm == ("xref",):
            sidm_match = SIDM_SCAN_PATTERN.search(value)
            if not sidm_match:
                continue
            normalized = normalize_identifier(sidm_match.group(0))
        elif path_norm in _TRUSTED_SIDM_PATHS:
            normalized = normalize_identifier(value)
            if detect_identifier_type(normalized) != "sidm":
                continue
        else:
            continue
        key = (normalized, source_field)
        if key in seen:
            continue
        seen.add(key)
        sources.append({"sidm": normalized, "source_field": source_field, "raw_value": str(value)})
    return sources


def _extract_trusted_ach_aliases(record: Dict[str, Any]) -> List[str]:
    aliases: List[str] = []
    record_id = record.get("id")
    record_type = record.get("type")
    if record_type == "model" and record_id is not None:
        normalized_id = normalize_identifier(record_id)
        if detect_identifier_type(normalized_id) == "model_id" and normalized_id not in aliases:
            aliases.append(normalized_id)
    for path, value in _walk_record_fields(record):
        if not value or not _path_matches(path, _ACH_FIELD_KEYS):
            continue
        normalized = normalize_identifier(value)
        if detect_identifier_type(normalized) == "model_id" and normalized not in aliases:
            aliases.append(normalized)
    return aliases


def _extract_bridge_ids_from_record(record: Dict[str, Any]) -> List[str]:
    bridge_ids: List[str] = []
    for path, value in _walk_record_fields(record):
        if not value:
            continue
        if _path_matches(path, _ACH_FIELD_KEYS):
            normalized = normalize_identifier(value)
            if normalized.startswith("ACH-") and normalized not in bridge_ids:
                bridge_ids.append(normalized)
        if _path_matches(path, _CVCL_FIELD_KEYS | _RRID_FIELD_KEYS) or tuple(_normalize_path_tokens(path)) == ("xref",):
            cvcl = extract_cvcl(value)
            if cvcl and cvcl not in bridge_ids:
                bridge_ids.append(cvcl)
            normalized = normalize_identifier(value)
            if normalized.startswith("RRID:CVCL_") and normalized not in bridge_ids:
                bridge_ids.append(normalized)
        if tuple(_normalize_path_tokens(path)) == ("xref",):
            for token in ACH_SCAN_PATTERN.findall(value):
                normalized = normalize_identifier(token)
                if normalized.startswith("ACH-") and normalized not in bridge_ids:
                    bridge_ids.append(normalized)
            for token in CVCL_SCAN_PATTERN.findall(value):
                normalized = normalize_identifier(token)
                if normalized.startswith("CVCL_") and normalized not in bridge_ids:
                    bridge_ids.append(normalized)
    return bridge_ids


def _extract_name_keys_from_record(record: Dict[str, Any]) -> List[str]:
    names: List[str] = []
    for path, value in _walk_record_fields(record):
        if not value or not _path_matches(path, _NAME_FIELD_KEYS):
            continue
        key = normalize_name_key(value)
        if key and key not in names:
            names.append(key)
    return names


def _extract_identifier_tokens_from_record(record: Dict[str, Any], identifier_type: str) -> List[str]:
    tokens: List[str] = []
    for path, value in _walk_record_fields(record):
        if not value:
            continue
        path_norm = tuple(_normalize_path_tokens(path))
        if identifier_type == "model_id" and _path_matches(path, _ACH_FIELD_KEYS):
            normalized = normalize_identifier(value)
            if normalized.startswith("ACH-") and normalized not in tokens:
                tokens.append(normalized)
        elif identifier_type in {"cvcl", "rrid"} and path_norm in {("cvcl",), ("rrid",), ("id",), ("ac",), ("accession",), ("model", "cvcl"), ("model", "rrid"), ("model", "id"), ("model", "ac"), ("model", "accession") }:
            normalized = normalize_identifier(value)
            cvcl = extract_cvcl(value)
            for candidate in (normalized, cvcl):
                if candidate and candidate not in tokens:
                    tokens.append(candidate)
        elif identifier_type == "name" and _path_matches(path, _NAME_FIELD_KEYS):
            key = normalize_name_key(value)
            if key and key not in tokens:
                tokens.append(key)
        elif identifier_type == "sidm" and path_norm in _TRUSTED_SIDM_PATHS:
            normalized = normalize_identifier(value)
            if normalized.startswith("SIDM") and normalized not in tokens:
                tokens.append(normalized)
    return tokens


@dataclass
class CellLineResolutionResult:
    input_raw: str
    input_normalized: str
    detected_type: str
    status: str
    sidm: Optional[str] = None
    matched_on: Optional[str] = None
    candidate_sidms: Optional[List[str]] = None
    source: Optional[str] = None
    raw_hit_count: int = 0
    exact_matched_record_count: int = 0
    matched_records: Optional[List[Dict[str, Any]]] = None
    sidm_source_fields: Optional[List[Dict[str, str]]] = None
    linked_identifier_records: Optional[List[Dict[str, Any]]] = None


class JsonHttpClient:
    """Small JSON HTTP client for resolver API calls."""

    def __init__(self):
        self._cache: Dict[Tuple[str, Tuple[Tuple[str, str], ...], int], Any] = {}

    def get_json(self, url: str, params: Optional[Dict[str, str]] = None, timeout: int = 10):
        params_key = tuple(sorted((params or {}).items()))
        cache_key = (url, params_key, timeout)
        if cache_key in self._cache:
            return self._cache[cache_key]

        if params:
            query = urlencode(params)
            separator = "&" if "?" in url else "?"
            url = f"{url}{separator}{query}"

        request = Request(url=url, headers={"Accept": "application/json"}, method="GET")
        try:
            with urlopen(request, timeout=timeout) as response:
                raw = response.read().decode("utf-8")
                if not raw:
                    self._cache[cache_key] = None
                    return None
                payload = json.loads(raw)
                self._cache[cache_key] = payload
                return payload
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
            self._cache[cache_key] = None
            return None

    def get_cache_stats(self) -> Dict[str, int]:
        return {"request_cache_enabled": True, "cached_requests": len(self._cache)}


def _walk_leaf_strings(value):
    if isinstance(value, dict):
        for key, sub_value in value.items():
            yield str(key)
            yield from _walk_leaf_strings(sub_value)
    elif isinstance(value, list):
        for item in value:
            yield from _walk_leaf_strings(item)
    else:
        yield str(value)


def _extract_matches(payload, pattern: re.Pattern) -> List[str]:
    matches = []
    for text in _walk_leaf_strings(payload):
        matches.extend(pattern.findall(text))
    seen = set()
    out = []
    for match in matches:
        normalized = normalize_identifier(match)
        if normalized and normalized not in seen:
            seen.add(normalized)
            out.append(normalized)
    return out


class SangerApiResolver:
    """Resolve identifiers against Sanger Cell Model Passports API."""

    DEFAULT_BASE_URLS = (
        "https://api.cellmodelpassports.sanger.ac.uk",
        "https://cellmodelpassports.sanger.ac.uk/api",
    )

    def __init__(
        self,
        http_client: Optional[JsonHttpClient] = None,
        base_urls: Optional[List[str]] = None,
        timeout: int = 10,
        deep_debug: bool = False,
    ):
        self.http = http_client or JsonHttpClient()
        self.base_urls = tuple(base_urls or self.DEFAULT_BASE_URLS)
        self.timeout = timeout
        self.deep_debug = deep_debug

    def _candidate_requests(self, normalized: str, identifier_type: str):
        """Generate candidate requests in priority order.
        
        For CVCL/RRID: high-confidence endpoints first to reduce API calls.
        - First try: /models?q=<CVCL>&include=identifiers (best match)
        - Then: Detail lookups by direct path
        - Finally: Fallback search endpoints
        """
        requests = []
        include_param = {"include": "identifiers"} if identifier_type in {"cvcl", "rrid"} else None
        
        search_candidates = [normalized]
        if identifier_type == "rrid":
            cvcl = extract_cvcl(normalized)
            if cvcl:
                search_candidates.append(cvcl)

        # For CVCL/RRID: try high-confidence endpoints first
        if identifier_type in {"cvcl", "rrid"}:
            for base_url in self.base_urls:
                # High priority: /models?q=<CVCL>&include=identifiers
                for token in search_candidates:
                    requests.append((f"{base_url}/models", {"q": token, "include": "identifiers"}))
            
            # Medium priority: Direct path lookups
            detail_paths = [
                f"/models/{normalized}",
                f"/model/{normalized}",
            ]
            for base_url in self.base_urls:
                for path in detail_paths:
                    requests.append((f"{base_url}{path}", include_param))
            
            # Lower priority: Other endpoints
            for base_url in self.base_urls:
                for token in search_candidates:
                    requests.append((f"{base_url}/search", {"q": token, "include": "identifiers"}))
                    requests.append((f"{base_url}/models", {"search": token, "include": "identifiers"}))
                    requests.append((f"{base_url}/passports", {"search": token, "include": "identifiers"}))
            
            # Fallback: Less likely to match
            fallback_paths = [
                f"/passports/{normalized}",
                f"/cell-lines/{normalized}",
                f"/cell_lines/{normalized}",
            ]
            for base_url in self.base_urls:
                for path in fallback_paths:
                    requests.append((f"{base_url}{path}", include_param))
                for token in search_candidates:
                    requests.append((f"{base_url}/model_list", {"search": token, "include": "identifiers"}))
        else:
            # For other identifier types: use original order (no optimization)
            detail_paths = [
                f"/models/{normalized}",
                f"/model/{normalized}",
                f"/passports/{normalized}",
                f"/cell-lines/{normalized}",
                f"/cell_lines/{normalized}",
            ]
            for base_url in self.base_urls:
                for path in detail_paths:
                    requests.append((f"{base_url}{path}", None))
                for token in search_candidates:
                    requests.append((f"{base_url}/search", {"q": token}))
                    requests.append((f"{base_url}/models", {"search": token}))
                    requests.append((f"{base_url}/models", {"q": token}))
                    requests.append((f"{base_url}/passports", {"search": token}))
                    requests.append((f"{base_url}/model_list", {"search": token}))
        
        return requests

    def resolve_one(self, identifier: object) -> CellLineResolutionResult:
        raw = "" if identifier is None else str(identifier)
        normalized = normalize_identifier(raw)
        identifier_type = detect_identifier_type(normalized)

        if not normalized:
            return CellLineResolutionResult(raw, normalized, identifier_type, "unresolved", source="sanger")

        # Canonical SIDM provided by user is accepted directly.
        if identifier_type == "sidm":
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=normalized,
                matched_on="sidm",
                source="sanger",
            )

        if identifier_type in {"cvcl", "rrid"}:
            raw_hit_count = 0
            matched_records: List[Dict[str, Any]] = []
            sidm_sources: List[Dict[str, str]] = []
            matched_identifier_ids: List[str] = []
            linked_identifier_records: List[Dict[str, Any]] = []
            printed_first_exact_model = False
            cvcl_token = extract_cvcl(normalized)

            for url, params in self._candidate_requests(normalized, identifier_type):
                payload = self.http.get_json(url, params=params, timeout=self.timeout)
                if payload is None or not isinstance(payload, dict):
                    continue
                data_records, included_records = _jsonapi_payload_records(payload)
                raw_hit_count += len(data_records) + len(included_records)

                for included in included_records:
                    if included.get("type") != "model_identifier":
                        continue
                    included_value = _record_identifier_value(included)
                    if not included_value:
                        continue
                    normalized_value = normalize_identifier(included_value)
                    if normalized_value == normalized or (cvcl_token and extract_cvcl(normalized_value) == cvcl_token):
                        included_id = str(included.get("id"))
                        if included_id not in matched_identifier_ids:
                            matched_identifier_ids.append(included_id)

                for record in data_records:
                    relationship_ids = _record_relationship_identifier_ids(record)
                    if relationship_ids and any(identifier_id in matched_identifier_ids for identifier_id in relationship_ids):
                        matched_records.append(record)
                        if not printed_first_exact_model:
                            linked_identifier_records = [
                                included
                                for included in included_records
                                if included.get("type") == "model_identifier"
                                and str(included.get("id")) in relationship_ids
                            ]
                            model_sidms = _extract_sidms_from_record(record)
                            model_sidm = model_sidms[0] if model_sidms else None
                            if self.deep_debug:
                                print(f"[SIDM][debug][exact_model] model_sidm={model_sidm}")
                                print(f"[SIDM][debug][exact_model] model_record_type={record.get('type')} model_record_id={record.get('id')}")
                                print(f"[SIDM][debug][exact_model] linked_identifier_records={linked_identifier_records}")
                                for linked_record in linked_identifier_records:
                                    identifier_value = _record_identifier_value(linked_record)
                                    identifier_type_name = linked_record.get("type")
                                    identifier_name = None
                                    attributes = linked_record.get("attributes")
                                    if isinstance(attributes, dict):
                                        identifier_name = attributes.get("name") or attributes.get("label") or attributes.get("identifier")
                                    if identifier_name is None:
                                        identifier_name = linked_record.get("name")
                                    print(
                                        f"[SIDM][debug][exact_model] linked_identifier type={identifier_type_name} name={identifier_name} value={identifier_value}"
                                    )
                            printed_first_exact_model = True

            matched_records = _dedupe_records(matched_records)
            if matched_records:
                sidm_candidates = _unique_values([
                    sidm for record in matched_records for sidm in _extract_sidms_from_record(record)
                ])
                for record in matched_records:
                    sidm_sources.extend(_extract_trusted_sidm_sources(record))

                if len(sidm_candidates) == 1:
                    return CellLineResolutionResult(
                        input_raw=raw,
                        input_normalized=normalized,
                        detected_type=identifier_type,
                        status="resolved",
                        sidm=sidm_candidates[0],
                        matched_on=identifier_type,
                        source="sanger",
                        raw_hit_count=raw_hit_count,
                        exact_matched_record_count=len(matched_records),
                        candidate_sidms=sidm_candidates,
                        matched_records=matched_records,
                        sidm_source_fields=sidm_sources,
                        linked_identifier_records=linked_identifier_records,
                    )
                if len(sidm_candidates) > 1:
                    return CellLineResolutionResult(
                        input_raw=raw,
                        input_normalized=normalized,
                        detected_type=identifier_type,
                        status="ambiguous",
                        matched_on=identifier_type,
                        candidate_sidms=sidm_candidates,
                        source="sanger",
                        raw_hit_count=raw_hit_count,
                        exact_matched_record_count=len(matched_records),
                        matched_records=matched_records,
                        sidm_source_fields=sidm_sources,
                        linked_identifier_records=linked_identifier_records,
                    )

        payloads = [self.http.get_json(url, params=params, timeout=self.timeout) for url, params in self._candidate_requests(normalized, identifier_type)]
        payloads = [payload for payload in payloads if payload is not None]
        raw_hit_count, matched_records, all_sidm_matches = _resolve_records_from_payloads(payloads, normalized, identifier_type)
        sidm_sources: List[Dict[str, str]] = []
        for record in matched_records:
            sidm_sources.extend(_extract_trusted_sidm_sources(record))

        if len(all_sidm_matches) == 1:
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=all_sidm_matches[0],
                matched_on=identifier_type,
                source="sanger",
                raw_hit_count=raw_hit_count,
                exact_matched_record_count=len(matched_records),
                candidate_sidms=all_sidm_matches,
                matched_records=matched_records,
                sidm_source_fields=sidm_sources,
            )
        if len(all_sidm_matches) > 1:
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="ambiguous",
                matched_on=identifier_type,
                candidate_sidms=all_sidm_matches,
                source="sanger",
                raw_hit_count=raw_hit_count,
                exact_matched_record_count=len(matched_records),
                matched_records=matched_records,
                sidm_source_fields=sidm_sources,
            )
        return CellLineResolutionResult(
            input_raw=raw,
            input_normalized=normalized,
            detected_type=identifier_type,
            status="unresolved",
            source="sanger",
            raw_hit_count=raw_hit_count,
            exact_matched_record_count=len(matched_records),
            matched_records=matched_records,
            sidm_source_fields=sidm_sources,
        )


class CellosaurusApiResolver:
    """Resolve identifiers using Cellosaurus API.

    This resolver can sometimes provide direct SIDM cross-references and can also
    provide bridge identifiers (ACH/CVCL) for a second Sanger lookup.
    """

    def __init__(self, http_client: Optional[JsonHttpClient] = None, base_url: str = "https://api.cellosaurus.org", timeout: int = 10):
        self.http = http_client or JsonHttpClient()
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _query_payload_entries(self, normalized: str, identifier_type: str):
        payloads = []
        cvcl = extract_cvcl(normalized)
        if cvcl:
            detail_url = f"{self.base_url}/cell-line/{cvcl}"
            payload = self.http.get_json(
                detail_url,
                params={"fields": "id,ac,sy,xref", "format": "json"},
                timeout=self.timeout,
            )
            if payload is not None:
                payloads.append((detail_url, payload))

        search_url = f"{self.base_url}/search/cell-line"
        payload = self.http.get_json(
            search_url,
            params={"q": normalized, "fields": "id,ac,sy,xref", "format": "json"},
            timeout=self.timeout,
        )
        if payload is not None:
            payloads.append((search_url, payload))
        return payloads

    def _query_payloads(self, normalized: str, identifier_type: str):
        return [payload for _, payload in self._query_payload_entries(normalized, identifier_type)]

    def _query_exact_records(self, normalized: str, identifier_type: str):
        if identifier_type in {"cvcl", "rrid"}:
            raw_hit_count = 0
            matched_records: List[Dict[str, Any]] = []
            sidm_sources: List[Dict[str, str]] = []
            matched_identifier_ids: List[str] = []

            for _, payload in self._query_payload_entries(normalized, identifier_type):
                data_records, included_records = _jsonapi_payload_records(payload)
                raw_hit_count += len(data_records) + len(included_records)

                for included in included_records:
                    if included.get("type") != "model_identifier":
                        continue
                    included_value = _record_identifier_value(included)
                    if not included_value:
                        continue
                    normalized_value = normalize_identifier(included_value)
                    if normalized_value == normalized or (extract_cvcl(normalized_value) and extract_cvcl(normalized_value) == extract_cvcl(normalized)):
                        included_id = str(included.get("id"))
                        if included_id not in matched_identifier_ids:
                            matched_identifier_ids.append(included_id)

                for record in data_records:
                    if _record_relationship_identifier_ids(record) and any(identifier_id in matched_identifier_ids for identifier_id in _record_relationship_identifier_ids(record)):
                        matched_records.append(record)

            matched_records = _dedupe_records(matched_records)
            if not matched_records:
                # Fallback to direct field inspection for non-JSONAPI payload shapes.
                for _, payload in self._query_payload_entries(normalized, identifier_type):
                    records = _payload_to_records(payload)
                    raw_hit_count += len(records)
                    for record in records:
                        if _record_matches_identifier(record, normalized, identifier_type):
                            matched_records.append(record)
                matched_records = _dedupe_records(matched_records)

            sidm_candidates = _unique_values([
                sidm for record in matched_records for sidm in _extract_sidms_from_record(record)
            ])
            for record in matched_records:
                sidm_sources.extend(_extract_trusted_sidm_sources(record))
            return raw_hit_count, matched_records, sidm_candidates, sidm_sources

        raw_hit_count = 0
        matched_records: List[Dict[str, Any]] = []
        for url, payload in self._query_payload_entries(normalized, identifier_type):
            records = _payload_to_records(payload)
            raw_hit_count += len(records)
            if "/cell-line/" in url:
                matched_records.extend(records)
                continue
            for record in records:
                if _record_matches_identifier(record, normalized, identifier_type):
                    matched_records.append(record)

        matched_records = _dedupe_records(matched_records)
        sidm_candidates = _unique_values([
            sidm for record in matched_records for sidm in _extract_sidms_from_record(record)
        ])
        sidm_sources: List[Dict[str, str]] = []
        for record in matched_records:
            sidm_sources.extend(_extract_trusted_sidm_sources(record))
        return raw_hit_count, matched_records, sidm_candidates, sidm_sources

    def resolve_one(self, identifier: object) -> CellLineResolutionResult:
        raw = "" if identifier is None else str(identifier)
        normalized = normalize_identifier(raw)
        identifier_type = detect_identifier_type(normalized)

        if not normalized:
            return CellLineResolutionResult(raw, normalized, identifier_type, "unresolved", source="cellosaurus")

        raw_hit_count, matched_records, all_sidm_matches, sidm_sources = self._query_exact_records(normalized, identifier_type)

        if len(all_sidm_matches) == 1:
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=all_sidm_matches[0],
                matched_on=identifier_type,
                source="cellosaurus",
                raw_hit_count=raw_hit_count,
                exact_matched_record_count=len(matched_records),
                candidate_sidms=all_sidm_matches,
                matched_records=matched_records,
                sidm_source_fields=sidm_sources,
            )
        if len(all_sidm_matches) > 1:
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="ambiguous",
                matched_on=identifier_type,
                candidate_sidms=all_sidm_matches,
                source="cellosaurus",
                raw_hit_count=raw_hit_count,
                exact_matched_record_count=len(matched_records),
                matched_records=matched_records,
                sidm_source_fields=sidm_sources,
            )
        return CellLineResolutionResult(
            input_raw=raw,
            input_normalized=normalized,
            detected_type=identifier_type,
            status="unresolved",
            source="cellosaurus",
            raw_hit_count=raw_hit_count,
            exact_matched_record_count=len(matched_records),
            matched_records=matched_records,
            sidm_source_fields=sidm_sources,
        )

    def extract_bridge_identifiers(self, identifier: object) -> List[str]:
        """Extract possible bridge IDs (ACH/CVCL) from Cellosaurus payloads."""
        normalized = normalize_identifier(identifier)
        identifier_type = detect_identifier_type(normalized)
        _, matched_records, _, _ = self._query_exact_records(normalized, identifier_type)
        candidates: List[str] = []
        for record in matched_records:
            for token in _extract_bridge_ids_from_record(record):
                if token not in candidates:
                    candidates.append(token)
        return candidates


class CellLineCache:
    """Disk cache for resolved cell lines to avoid repeated API calls.
    
    Cache key is based on input file path + content hash.
    Stores: sidm_list, alias_to_sidm, resolution_report, cache_stats.
    """

    def __init__(self, cache_dir: Optional[str] = None):
        """Initialize cache manager.
        
        Args:
            cache_dir: Directory for cache files. If None, uses ~/.celios/cache/resolution/
        """
        if cache_dir is None:
            cache_dir = os.path.join(Path.home(), ".celios", "cache", "resolution")
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _hash_input(self, file_path: str, df: pd.DataFrame) -> str:
        """Generate cache key from file path and content hash."""
        # Combine file path and dataframe content
        path_hash = hashlib.md5(file_path.encode()).hexdigest()
        # Hash rows (without full content, just shape + first/last row)
        df_key = f"{df.shape}:{df.iloc[0].to_string() if len(df) > 0 else 'empty'}"
        df_hash = hashlib.md5(df_key.encode()).hexdigest()
        return f"{path_hash}_{df_hash}"

    def get(self, file_path: str, df: pd.DataFrame) -> Optional[Dict[str, Any]]:
        """Load cache if exists and valid.
        
        Returns:
            Dict with keys: sidm_list, alias_to_sidm, resolution_report, cache_stats
            None if cache miss or error.
        """
        cache_key = self._hash_input(file_path, df)
        cache_file = self.cache_dir / f"{cache_key}.json"
        
        if not cache_file.exists():
            return None
        
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except Exception:
            return None

    def put(self, file_path: str, df: pd.DataFrame, data: Dict[str, Any]) -> bool:
        """Save resolved cell lines to cache.
        
        Args:
            file_path: Path to input file
            df: Input dataframe
            data: Dict with keys: sidm_list, alias_to_sidm, resolution_report, cache_stats
        
        Returns:
            True if saved successfully, False otherwise.
        """
        cache_key = self._hash_input(file_path, df)
        cache_file = self.cache_dir / f"{cache_key}.json"
        
        try:
            # Ensure data is JSON serializable
            clean_data = {
                'sidm_list': list(data.get('sidm_list', [])),
                'alias_to_sidm': dict(data.get('alias_to_sidm', {})),
                'resolution_report': str(data.get('resolution_report', '')),
                'cache_stats': dict(data.get('cache_stats', {})),
            }
            with open(cache_file, 'w') as f:
                json.dump(clean_data, f, indent=2)
            return True
        except Exception:
            return False


class CellLineResolver:
    """Composite resolver: Sanger API primary, Cellosaurus fallback.
    
    Includes in-memory caching of resolved identifiers to avoid repeated
    API calls during batch operations.
    """

    def __init__(
        self,
        use_cellosaurus_fallback: bool = True,
        http_client: Optional[JsonHttpClient] = None,
        sanger_base_urls: Optional[List[str]] = None,
        cellosaurus_base_url: str = "https://api.cellosaurus.org",
        timeout: int = 10,
        enable_cache: bool = True,
        deep_debug: bool = False,
    ):
        self.sanger = SangerApiResolver(
            http_client=http_client,
            base_urls=sanger_base_urls,
            timeout=timeout,
            deep_debug=deep_debug,
        )
        self.cellosaurus = CellosaurusApiResolver(
            http_client=http_client,
            base_url=cellosaurus_base_url,
            timeout=timeout,
        )
        self.use_cellosaurus_fallback = use_cellosaurus_fallback
        self._cache: Dict[str, CellLineResolutionResult] = {} if enable_cache else None
        self.sidm_to_sidm: Dict[str, str] = {}
        self.cvcl_to_sidm: Dict[str, str] = {}
        self.rrid_to_sidm: Dict[str, str] = {}
        self.ach_to_sidm: Dict[str, str] = {}
        self.name_to_sidm: Dict[str, str] = {}

    def _register_alias(self, lookup: Dict[str, str], alias: str, sidm: str) -> None:
        if not alias or not sidm:
            return
        if alias not in lookup:
            lookup[alias] = sidm

    def _register_result(self, result: CellLineResolutionResult) -> None:
        if result.status != "resolved" or not result.sidm:
            return

        sidm = result.sidm
        self._register_alias(self.sidm_to_sidm, sidm, sidm)
        self._register_alias(self.sidm_to_sidm, result.input_normalized, sidm)
        self._register_alias(self.sidm_to_sidm, result.input_raw, sidm)

        if result.detected_type == "model_id":
            self._register_alias(self.ach_to_sidm, result.input_normalized, sidm)
            self._register_alias(self.ach_to_sidm, normalize_identifier(result.input_raw), sidm)
        elif result.detected_type == "rrid":
            self._register_alias(self.rrid_to_sidm, result.input_normalized, sidm)
            cvcl = extract_cvcl(result.input_normalized)
            if cvcl:
                self._register_alias(self.cvcl_to_sidm, cvcl, sidm)
                self._register_alias(self.rrid_to_sidm, f"RRID:{cvcl}", sidm)
        elif result.detected_type == "cvcl":
            cvcl = extract_cvcl(result.input_normalized) or result.input_normalized
            self._register_alias(self.cvcl_to_sidm, cvcl, sidm)
            self._register_alias(self.rrid_to_sidm, f"RRID:{cvcl}", sidm)
        elif result.detected_type == "name":
            self._register_alias(self.name_to_sidm, normalize_name_key(result.input_normalized), sidm)
            self._register_alias(self.name_to_sidm, normalize_name_key(result.input_raw), sidm)

        for record in result.matched_records or []:
            for token in _extract_identifier_tokens_from_record(record, result.detected_type):
                if result.detected_type == "model_id" and token.startswith("ACH-"):
                    self._register_alias(self.ach_to_sidm, token, sidm)
                elif result.detected_type in {"cvcl", "rrid"}:
                    if token.startswith("RRID:CVCL_"):
                        self._register_alias(self.rrid_to_sidm, token, sidm)
                        cvcl = extract_cvcl(token)
                        if cvcl:
                            self._register_alias(self.cvcl_to_sidm, cvcl, sidm)
                    elif token.startswith("CVCL_"):
                        self._register_alias(self.cvcl_to_sidm, token, sidm)
                        self._register_alias(self.rrid_to_sidm, f"RRID:{token}", sidm)
                elif result.detected_type == "name":
                    self._register_alias(self.name_to_sidm, token, sidm)
                elif result.detected_type == "sidm" and token.startswith("SIDM"):
                    self._register_alias(self.sidm_to_sidm, token, sidm)

        for record in result.linked_identifier_records or []:
            linked_value = _record_identifier_value(record)
            if not linked_value:
                continue
            normalized_linked_value = normalize_identifier(linked_value)
            if normalized_linked_value.startswith("ACH-"):
                self._register_alias(self.ach_to_sidm, normalized_linked_value, sidm)

        if result.detected_type in {"cvcl", "rrid"}:
            cvcl = extract_cvcl(result.input_normalized)
            if cvcl:
                self._register_alias(self.cvcl_to_sidm, cvcl, sidm)
                self._register_alias(self.rrid_to_sidm, f"RRID:{cvcl}", sidm)

    def _lookup_local(self, normalized: str, identifier_type: str) -> Optional[str]:
        if identifier_type == "sidm":
            return self.sidm_to_sidm.get(normalized) or normalized
        if identifier_type == "model_id":
            return self.ach_to_sidm.get(normalized)
        if identifier_type == "rrid":
            return self.rrid_to_sidm.get(normalized) or self.cvcl_to_sidm.get(extract_cvcl(normalized) or "")
        if identifier_type == "cvcl":
            cvcl = extract_cvcl(normalized) or normalized
            return self.cvcl_to_sidm.get(cvcl) or self.rrid_to_sidm.get(f"RRID:{cvcl}")
        if identifier_type == "name":
            return self.name_to_sidm.get(normalize_name_key(normalized))
        return None

    def get_alias_to_sidm(self) -> Dict[str, str]:
        alias_to_sidm: Dict[str, str] = {}
        for lookup in (self.sidm_to_sidm, self.cvcl_to_sidm, self.rrid_to_sidm, self.ach_to_sidm, self.name_to_sidm):
            for alias, sidm in lookup.items():
                if alias and alias not in alias_to_sidm:
                    alias_to_sidm[alias] = sidm
        return alias_to_sidm

    def resolve_one(self, identifier: object) -> CellLineResolutionResult:
        raw = "" if identifier is None else str(identifier)
        normalized = normalize_identifier(raw)
        identifier_type = detect_identifier_type(normalized)

        if not normalized:
            result = CellLineResolutionResult(raw, normalized, identifier_type, "unresolved", source="sanger")
            if self._cache is not None:
                self._cache[normalized] = result
            return result

        if self._cache is not None and normalized in self._cache:
            return self._cache[normalized]

        local_sidm = self._lookup_local(normalized, identifier_type)
        if local_sidm:
            result = CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=local_sidm,
                matched_on=identifier_type,
                source="cache",
                candidate_sidms=[local_sidm],
                raw_hit_count=0,
                exact_matched_record_count=0,
                sidm_source_fields=[],
            )
            if self._cache is not None:
                self._cache[normalized] = result
            return result

        primary = self.sanger.resolve_one(identifier)
        if primary.status == "resolved" or primary.status == "ambiguous":
            if primary.status == "resolved":
                self._register_result(primary)
            if self._cache is not None:
                self._cache[normalized] = primary
            return primary

        if not self.use_cellosaurus_fallback:
            if self._cache is not None:
                self._cache[normalized] = primary
            return primary

        fallback = self.cellosaurus.resolve_one(identifier)
        if fallback.status == "resolved" or fallback.status == "ambiguous":
            if fallback.status == "resolved":
                self._register_result(fallback)
            if self._cache is not None:
                self._cache[normalized] = fallback
            return fallback

        bridge_candidates = []
        for token in self.cellosaurus.extract_bridge_identifiers(identifier):
            bridged = self.sanger.resolve_one(token)
            if bridged.status == "resolved" and bridged.sidm:
                if bridged.sidm not in bridge_candidates:
                    bridge_candidates.append(bridged.sidm)
                self._register_result(bridged)
            elif bridged.status == "ambiguous" and bridged.candidate_sidms:
                for sidm in bridged.candidate_sidms:
                    if sidm not in bridge_candidates:
                        bridge_candidates.append(sidm)

        if len(bridge_candidates) == 1:
            result = CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=bridge_candidates[0],
                matched_on=identifier_type,
                source="cellosaurus->sanger",
                candidate_sidms=bridge_candidates,
            )
            if self._cache is not None:
                self._cache[normalized] = result
            self._register_result(result)
            return result
        if len(bridge_candidates) > 1:
            result = CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="ambiguous",
                matched_on=identifier_type,
                candidate_sidms=bridge_candidates,
                source="cellosaurus->sanger",
            )
            if self._cache is not None:
                self._cache[normalized] = result
            return result

        if self._cache is not None:
            self._cache[normalized] = primary
        return primary
    
    def clear_cache(self) -> None:
        """Clear the in-memory resolution cache."""
        if self._cache is not None:
            self._cache.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """Return cache statistics (e.g., for logging)."""
        request_stats = self.sanger.http.get_cache_stats() if hasattr(self.sanger.http, "get_cache_stats") else {"cached_requests": 0}
        if self._cache is None:
            return {
                "cache_enabled": False,
                "cached_identifiers": 0,
                "cached_requests": request_stats.get("cached_requests", 0),
                "sidm_to_sidm": len(self.sidm_to_sidm),
                "cvcl_to_sidm": len(self.cvcl_to_sidm),
                "rrid_to_sidm": len(self.rrid_to_sidm),
                "ach_to_sidm": len(self.ach_to_sidm),
                "name_to_sidm": len(self.name_to_sidm),
            }
        return {
            "cache_enabled": True,
            "cached_identifiers": len(self._cache),
            "cached_requests": request_stats.get("cached_requests", 0),
            "sidm_to_sidm": len(self.sidm_to_sidm),
            "cvcl_to_sidm": len(self.cvcl_to_sidm),
            "rrid_to_sidm": len(self.rrid_to_sidm),
            "ach_to_sidm": len(self.ach_to_sidm),
            "name_to_sidm": len(self.name_to_sidm),
        }


def resolve_identifiers_to_sidm(
    identifiers: List[object], model_registry: Optional[str] = None
) -> Tuple[Dict[str, str], List[str], List[CellLineResolutionResult]]:
    """Resolve arbitrary identifiers to SIDM.

    Returns:
      - sidm_dict: {SIDM -> original_identifier}
      - not_found: unresolved input identifiers
      - results: detailed per-input resolution results
    """
    resolver = CellLineResolver(use_cellosaurus_fallback=True)
    sidm_dict: Dict[str, str] = {}
    not_found: List[str] = []
    results: List[CellLineResolutionResult] = []

    for value in identifiers:
        try:
            result = resolver.resolve_one(value)
        except Exception:
            result = CellLineResolutionResult(
                input_raw="" if value is None else str(value),
                input_normalized=normalize_identifier(value),
                detected_type=detect_identifier_type(value),
                status="unresolved",
                source="sanger",
            )
        results.append(result)
        if result.status == "resolved" and result.sidm is not None:
            sidm_dict[result.sidm] = result.input_raw
        else:
            not_found.append(result.input_raw)

    return sidm_dict, not_found, results


def resolve_model_ids_to_sidm(
    model_ids: List[object], model_registry: Optional[str] = None
) -> Tuple[Dict[str, str], List[str]]:
    """Resolve DepMap ModelIDs to SIDM.

    If model_registry is provided, it is used directly for deterministic local
    mapping (backward-compatible behavior for explicit custom registries).
    Otherwise the online resolver is used.

    Returns:
      - model_to_sidm: {ModelID -> SIDM}
      - not_found: unresolved ModelID values
    """
    model_to_sidm: Dict[str, str] = {}
    not_found: List[str] = []
    seen_missing = set()

    if model_registry:
        registry_df = pd.read_csv(model_registry)
        if "ModelID" not in registry_df.columns or "SangerModelID" not in registry_df.columns:
            raise ValueError("Model registry must contain 'ModelID' and 'SangerModelID' columns")

        lookup = {}
        for _, row in registry_df[["ModelID", "SangerModelID"]].dropna().iterrows():
            model_id = normalize_identifier(row["ModelID"])
            sidm = normalize_identifier(row["SangerModelID"])
            if model_id and sidm:
                lookup[model_id] = sidm

        for value in model_ids:
            normalized_model_id = normalize_identifier(value)
            sidm = lookup.get(normalized_model_id)
            if sidm:
                model_to_sidm[normalized_model_id] = sidm
            elif normalized_model_id not in seen_missing:
                not_found.append(normalized_model_id)
                seen_missing.add(normalized_model_id)

        return model_to_sidm, not_found

    resolver = CellLineResolver(use_cellosaurus_fallback=True)
    for value in model_ids:
        normalized_model_id = normalize_identifier(value)
        try:
            result = resolver.resolve_one(normalized_model_id)
        except Exception:
            result = CellLineResolutionResult(
                input_raw=str(value),
                input_normalized=normalized_model_id,
                detected_type=detect_identifier_type(normalized_model_id),
                status="unresolved",
                source="sanger",
            )
        if result.status == "resolved" and result.sidm is not None:
            model_to_sidm[normalized_model_id] = result.sidm
        elif normalized_model_id not in seen_missing:
            not_found.append(normalized_model_id)
            seen_missing.add(normalized_model_id)

    return model_to_sidm, not_found


def resolve_sidm_from_dataframe(
    df: pd.DataFrame, 
    model_registry: Optional[str] = None, 
    verbose: bool = False,
    cell_line_file: Optional[str] = None,
    enable_cache: bool = True,
    cache_dir: Optional[str] = None,
    deep_debug: bool = False,
) -> Tuple[Dict[str, str], List[str], Dict[str, int]]:
    """Resolve SIDMs row-wise from a user-provided cell-line table.

    Resolution strategy per row:
    1) Try best available identifier column by priority
    2) If unresolved, try display name columns as fallback

    Args:
        df: Input dataframe with cell line identifiers
        model_registry: Optional model registry URL
        verbose: Print resolution details (summary only)
        cell_line_file: Path to input file (used for cache key generation)
        enable_cache: Whether to use disk cache
        cache_dir: Directory for cache files (default: ~/.celios/cache/resolution/)
        deep_debug: Print per-row details (verbose must also be True)

    Returns:
      - sidm_dict: {SIDM -> display_name}
      - not_found: identifiers/names that could not be resolved
      - resolution_report: {
          'total_rows': int,
          'resolved': int,
          'ambiguous': int,
          'unresolved': int,
          'cache_stats': dict
        }
    """
    # Check cache first if enabled and input file provided
    cache = None
    if enable_cache and cell_line_file:
        cache = CellLineCache(cache_dir)
        cached_result = cache.get(cell_line_file, df)
        if cached_result:
            if verbose:
                print(f"[SIDM] Cache hit for {cell_line_file}")
            sidm_dict = {}
            for sidm, name in cached_result['alias_to_sidm'].items():
                if sidm.startswith('SIDM'):
                    sidm_dict[sidm] = name
            return sidm_dict, [], cached_result['cache_stats']
    
    id_columns_priority = [
        "SIDM",
        "sidm",
        "SangerModelID",
        "ModelID",
        "model_id",
        "ACH",
        "RRID",
        "rrid",
        "CVCL",
        "cvcl",
        "CCLE_ID",
        "CCLEName",
        "cell_line_name",
        "CellLineName",
        "NAME",
        "name",
    ]
    display_columns_priority = ["cell_line_name", "CellLineName", "NAME", "name"]

    available_id_cols = [col for col in id_columns_priority if col in df.columns]
    if not available_id_cols:
        raise ValueError(
            "cell_line_file must contain at least one identifier column. "
            "Accepted columns include SIDM, ModelID, RRID, CVCL, CCLE_ID, or cell_line_name."
        )

    # Check if SIDM column exists - if so, use directly without API calls (Todo 3)
    has_sidm_column = any(col for col in id_columns_priority if col in {"SIDM", "sidm"} and col in df.columns)
    
    if has_sidm_column and verbose:
        print("[SIDM] SIDM column detected - using directly without API calls")

    resolver = CellLineResolver(use_cellosaurus_fallback=True)
    sidm_dict: Dict[str, str] = {}
    not_found: List[str] = []
    # Track resolution statistics
    resolution_counts = {
        'total_rows': len(df),
        'resolved': 0,
        'ambiguous': 0,
        'unresolved': 0,
    }

    # Keep detailed per-row results and alias map for resolved rows
    resolution_results = []
    cvcl_debug_emitted = False

    def _print_cvcl_record_debug(query_identifier: str, row_index: int, deep_debug: bool = False) -> None:
        if not deep_debug:
            return
        normalized_query = normalize_identifier(query_identifier)
        identifier_type = detect_identifier_type(normalized_query)
        raw_records: List[Dict[str, Any]] = []
        for url, params in resolver.sanger._candidate_requests(normalized_query, identifier_type):
            payload = resolver.sanger.http.get_json(url, params=params, timeout=resolver.sanger.timeout)
            if payload is None:
                continue
            raw_records.extend(_payload_to_records(payload))
        if not raw_records:
            print(f"[SIDM][debug][row={row_index}] no raw records found for {query_identifier}")
            return

        keywords = ("cvcl", "cellosaurus", "rrid", "xref", "identifier", "accession")
        #print(f"[SIDM][debug][row={row_index}] raw_record_count={len(raw_records)}")
        for record_index, record in enumerate(raw_records[:2], start=1):
            matching_fields: List[Dict[str, str]] = []
            for path, value in _walk_record_fields(record):
                path_text = ".".join(path).lower()
                value_text = str(value).lower()
                if any(keyword in path_text or keyword in value_text for keyword in keywords):
                    matching_fields.append({"field": ".".join(path), "value": str(value)})
            #print(f"[SIDM][debug][row={row_index}] raw_record_{record_index}_top_level_keys={top_level_keys}")
            #print(f"[SIDM][debug][row={row_index}] raw_record_{record_index}_matching_fields={matching_fields}")
    for row_idx, row in df.iterrows():
        row_identifier = None
        for col in available_id_cols:
            value = row[col]
            if pd.notna(value) and str(value).strip():
                row_identifier = str(value).strip()
                break

        if row_identifier is None:
            resolution_counts['unresolved'] += 1
            resolution_results.append(
                CellLineResolutionResult(
                    input_raw="",
                    input_normalized="",
                    detected_type="name",
                    status="unresolved",
                    source="sanger",
                )
            )
            if verbose and deep_debug:
                print(
                    f"[SIDM][row={row_idx}] input_cell_line_name='' input_identifier='' "
                    f"detected_type=name raw_hit_count=0 exact_matched_record_count=0 sidm_candidates=[] status=unresolved"
                )
            continue

        display_name = row_identifier
        for col in display_columns_priority:
            if col not in df.columns:
                continue
            value = row[col]
            if pd.notna(value) and str(value).strip():
                display_name = str(value).strip()
                break

        # If SIDM column exists and current identifier is SIDM, use it directly (Todo 3)
        normalized_id = normalize_identifier(row_identifier)
        id_type = detect_identifier_type(normalized_id)
        
        if has_sidm_column and id_type == "sidm":
            result = CellLineResolutionResult(
                input_raw=row_identifier,
                input_normalized=normalized_id,
                detected_type="sidm",
                status="resolved",
                sidm=normalized_id,
                matched_on="sidm",
                source="sanger",
            )
        else:
            # API resolution for non-SIDM identifiers
            if verbose and not cvcl_debug_emitted and id_type in {"cvcl", "rrid"}:
                _print_cvcl_record_debug(row_identifier, row_idx, deep_debug=False)
                cvcl_debug_emitted = True

            try:
                result = resolver.resolve_one(row_identifier)
            except Exception as exc:
                if verbose and deep_debug:
                    print(f"[SIDM][row={row_idx}] resolution error for '{row_identifier}': {exc}")
                result = CellLineResolutionResult(
                    input_raw=row_identifier,
                    input_normalized=normalized_id,
                    detected_type=id_type,
                    status="unresolved",
                    source="sanger",
                )

            if result.status == "unresolved" and display_name != row_identifier:
                try:
                    fallback_result = resolver.resolve_one(display_name)
                except Exception as exc:
                    if verbose and deep_debug:
                        print(f"[SIDM][row={row_idx}] display-name fallback failed for '{display_name}': {exc}")
                    fallback_result = None
                if fallback_result is not None and fallback_result.status != "unresolved":
                    result = fallback_result

        resolution_results.append(result)

        if result.status == "resolved" and result.sidm is not None:
            sidm_dict[result.sidm] = display_name
            resolution_counts['resolved'] += 1
            if verbose and deep_debug:
                print(
                    f"[SIDM][row={row_idx}] input_cell_line_name='{display_name}' "
                    f"input_identifier='{row_identifier}' detected_type={result.detected_type} "
                    f"raw_hit_count={getattr(result, 'raw_hit_count', 0)} "
                    f"exact_matched_record_count={getattr(result, 'exact_matched_record_count', 0)} "
                    f"sidm_candidates={(result.candidate_sidms or [result.sidm])} status={result.status}"
                )
        elif result.status == "ambiguous":
            resolution_counts['ambiguous'] += 1
            not_found.append(row_identifier)
            if verbose and deep_debug:
                matched_record_keys = [list(record.keys()) for record in (result.matched_records or [])]
                print(
                    f"[SIDM][row={row_idx}] input_cell_line_name='{display_name}' "
                    f"input_identifier='{row_identifier}' detected_type={result.detected_type} "
                    f"raw_hit_count={getattr(result, 'raw_hit_count', 0)} "
                    f"exact_matched_record_count={getattr(result, 'exact_matched_record_count', 0)} "
                    f"sidm_candidates={result.candidate_sidms or []} status={result.status}"
                )
                print(f"[SIDM][row={row_idx}] matched_record_keys={matched_record_keys}")
                print(f"[SIDM][row={row_idx}] sidm_source_fields={result.sidm_source_fields or []}")
        else:
            resolution_counts['unresolved'] += 1
            not_found.append(row_identifier)
            if verbose and deep_debug:
                print(
                    f"[SIDM][row={row_idx}] input_cell_line_name='{display_name}' "
                    f"input_identifier='{row_identifier}' detected_type={result.detected_type} "
                    f"raw_hit_count={getattr(result, 'raw_hit_count', 0)} "
                    f"exact_matched_record_count={getattr(result, 'exact_matched_record_count', 0)} "
                    f"sidm_candidates={result.candidate_sidms or []} status={result.status}"
                )

    alias_to_sidm = resolver.get_alias_to_sidm()
    resolution_counts['cache_stats'] = resolver.get_cache_stats()
    resolution_counts['alias_to_sidm'] = alias_to_sidm
    resolution_counts['detailed_results'] = resolution_results

    # Save to disk cache if enabled
    if enable_cache and cache and cell_line_file:
        cache_data = {
            'sidm_list': list(sidm_dict.keys()),
            'alias_to_sidm': alias_to_sidm,
            'resolution_report': str(resolution_counts),
            'cache_stats': resolution_counts.get('cache_stats', {}),
        }
        if cache.put(cell_line_file, df, cache_data):
            if verbose and deep_debug:
                print(f"[SIDM] Saved resolution cache for {cell_line_file}")

    if verbose:
        # Simplified Step 2 verbose output (move details behind deep_debug)
        ach_aliases = sorted(alias for alias in alias_to_sidm if normalize_identifier(alias).startswith("ACH-"))
        print("[SIDM] Resolution summary:")
        print(f"  total_rows: {resolution_counts['total_rows']}")
        print(f"  resolved: {resolution_counts['resolved']}")
        print(f"  ambiguous: {resolution_counts['ambiguous']}")
        print(f"  unresolved: {resolution_counts['unresolved']}")
        print(f"  alias_count: {len(alias_to_sidm)}")
        print(f"  ach_alias_count: {len(ach_aliases)}")
        
        if deep_debug:
            # Detailed diagnostics only with deep_debug=True
            ach_aliases = sorted(alias for alias in alias_to_sidm if normalize_identifier(alias).startswith("ACH-"))
            sidm_alias_groups: Dict[str, List[str]] = {}
            for alias, sidm in alias_to_sidm.items():
                sidm_alias_groups.setdefault(sidm, []).append(alias)
            print(f"  ach_aliases_preview: {ach_aliases[:10]}")
            for sidm in sorted(sidm_alias_groups):
                print(f"  aliases_for_{sidm}: {sorted(sidm_alias_groups[sidm])}")
            print(f"  cache_stats: {resolution_counts['cache_stats']}")

    return sidm_dict, not_found, resolution_counts


def save_identifier_mapping(
    df: pd.DataFrame,
    resolution_report: Dict[str, Any],
    alias_to_sidm: Dict[str, str],
    output_path: str,
    verbose: bool = False,
) -> bool:
    """Save user-readable identifier mapping to CSV.
    
    Args:
        df: Input dataframe with cell line identifiers
        resolution_report: Resolution report dict containing detailed_results
        alias_to_sidm: Map of all aliases to SIDM
        output_path: Path to save identifiers.csv
        verbose: Print confirmation message
    
    Returns:
        True if saved successfully, False otherwise
    """
    try:
        detailed_results = resolution_report.get('detailed_results', [])
        if not detailed_results:
            return False
        
        rows = []
        
        for idx, result in enumerate(detailed_results):
            # Extract identifiers for this SIDM
            sidm = result.sidm if result else None
            
            # Find all aliases for this SIDM
            all_aliases = []
            if sidm:
                for alias, mapped_sidm in alias_to_sidm.items():
                    if mapped_sidm == sidm:
                        all_aliases.append(alias)
            
            # Extract specific identifier types
            cvcl = None
            rrid = None
            ach = None
            
            if result and result.detected_type == "cvcl":
                cvcl = extract_cvcl(result.input_normalized)
            elif result and result.detected_type == "rrid":
                cvcl = extract_cvcl(result.input_normalized)
                rrid = result.input_normalized
            elif result and result.detected_type == "model_id":
                ach = result.input_normalized
            
            # Build row
            row = {
                'input_cell_line_name': df.iloc[idx].get('cell_line_name') if 'cell_line_name' in df.columns else '',
                'input_identifier': result.input_raw if result else '',
                'detected_type': result.detected_type if result else 'unknown',
                'status': result.status if result else 'unresolved',
                'SIDM': sidm or '',
                'matched_on': result.matched_on if result else '',
                'source': result.source if result else '',
                'CVCL': cvcl or '',
                'RRID': rrid or '',
                'ACH': ach or '',
                'all_aliases': '; '.join(sorted(all_aliases)) if all_aliases else '',
            }
            rows.append(row)
        
        # Create dataframe and save
        mapping_df = pd.DataFrame(rows)
        mapping_df.to_csv(output_path, index=False)
        
        if verbose:
            print(f"[SIDM] Saved identifier mapping to: {output_path}")
        
        return True
    
    except Exception as exc:
        if verbose:
            print(f"[SIDM] Error saving identifier mapping: {exc}")
        return False


__all__ = [
    "CellLineCache",
    "CellLineResolutionResult",
    "CellLineResolver",
    "CellosaurusApiResolver",
    "JsonHttpClient",
    "SangerApiResolver",
    "detect_identifier_type",
    "extract_cvcl",
    "normalize_identifier",
    "normalize_name_key",
    "resolve_identifiers_to_sidm",
    "resolve_model_ids_to_sidm",
    "resolve_sidm_from_dataframe",
    "save_identifier_mapping",
]
