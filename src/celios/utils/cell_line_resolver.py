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
from typing import Dict, List, Optional, Tuple

import pandas as pd


SIDM_PATTERN = re.compile(r"^SIDM\s*[-_]?\s*(\d+)$", re.IGNORECASE)
ACH_PATTERN = re.compile(r"^ACH\s*[-_]?\s*(\d+)$", re.IGNORECASE)
CVCL_PATTERN = re.compile(r"^CVCL\s*[-_]?\s*([A-Z0-9]+)$", re.IGNORECASE)
SIDM_SCAN_PATTERN = re.compile(r"SIDM\d{5,}", re.IGNORECASE)
ACH_SCAN_PATTERN = re.compile(r"ACH-\d{6}", re.IGNORECASE)
CVCL_SCAN_PATTERN = re.compile(r"CVCL_[A-Z0-9]+", re.IGNORECASE)


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

    # Accept RRID variants and normalize to RRID:CVCL_XXXX
    # e.g. RRID CVCL_0023, RRID-CVCL0023, RRID:CVCL-0023
    if "RRID" in upper or "CVCL" in upper:
        cvcl = extract_cvcl(upper)
        if cvcl:
            if upper.startswith("RRID"):
                return f"RRID:{cvcl}"
            # Raw CVCL input goes through canonical CVCL form
            return cvcl

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
    return "name"


def normalize_name_key(value: object) -> str:
    """Build a comparison key for name/synonym matching."""
    if value is None:
        return ""
    text = unicodedata.normalize("NFKC", str(value)).upper().strip()
    return re.sub(r"[^A-Z0-9]", "", text)


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


class JsonHttpClient:
    """Small JSON HTTP client for resolver API calls."""

    def get_json(self, url: str, params: Optional[Dict[str, str]] = None, timeout: int = 10):
        if params:
            query = urlencode(params)
            separator = "&" if "?" in url else "?"
            url = f"{url}{separator}{query}"

        request = Request(url=url, headers={"Accept": "application/json"}, method="GET")
        try:
            with urlopen(request, timeout=timeout) as response:
                raw = response.read().decode("utf-8")
                if not raw:
                    return None
                return json.loads(raw)
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
            return None


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
    ):
        self.http = http_client or JsonHttpClient()
        self.base_urls = tuple(base_urls or self.DEFAULT_BASE_URLS)
        self.timeout = timeout

    def _candidate_requests(self, normalized: str, identifier_type: str):
        requests = []
        detail_paths = [
            f"/models/{normalized}",
            f"/model/{normalized}",
            f"/passports/{normalized}",
            f"/cell-lines/{normalized}",
            f"/cell_lines/{normalized}",
        ]
        search_candidates = [normalized]
        if identifier_type == "rrid":
            cvcl = extract_cvcl(normalized)
            if cvcl:
                search_candidates.append(cvcl)

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

        all_sidm_matches = []
        for url, params in self._candidate_requests(normalized, identifier_type):
            payload = self.http.get_json(url, params=params, timeout=self.timeout)
            if payload is None:
                continue
            sidms = _extract_matches(payload, SIDM_SCAN_PATTERN)
            if sidms:
                for sidm in sidms:
                    if sidm not in all_sidm_matches:
                        all_sidm_matches.append(sidm)

        if len(all_sidm_matches) == 1:
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=all_sidm_matches[0],
                matched_on=identifier_type,
                source="sanger",
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
            )
        return CellLineResolutionResult(
            input_raw=raw,
            input_normalized=normalized,
            detected_type=identifier_type,
            status="unresolved",
            source="sanger",
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

    def _query_payloads(self, normalized: str, identifier_type: str):
        payloads = []
        cvcl = extract_cvcl(normalized)
        if cvcl:
            payload = self.http.get_json(
                f"{self.base_url}/cell-line/{cvcl}",
                params={"fields": "id,ac,sy,xref", "format": "json"},
                timeout=self.timeout,
            )
            if payload is not None:
                payloads.append(payload)

        payload = self.http.get_json(
            f"{self.base_url}/search/cell-line",
            params={"q": normalized, "fields": "id,ac,sy,xref", "format": "json"},
            timeout=self.timeout,
        )
        if payload is not None:
            payloads.append(payload)
        return payloads

    def resolve_one(self, identifier: object) -> CellLineResolutionResult:
        raw = "" if identifier is None else str(identifier)
        normalized = normalize_identifier(raw)
        identifier_type = detect_identifier_type(normalized)

        if not normalized:
            return CellLineResolutionResult(raw, normalized, identifier_type, "unresolved", source="cellosaurus")

        all_sidm_matches = []
        for payload in self._query_payloads(normalized, identifier_type):
            sidms = _extract_matches(payload, SIDM_SCAN_PATTERN)
            for sidm in sidms:
                if sidm not in all_sidm_matches:
                    all_sidm_matches.append(sidm)

        if len(all_sidm_matches) == 1:
            return CellLineResolutionResult(
                input_raw=raw,
                input_normalized=normalized,
                detected_type=identifier_type,
                status="resolved",
                sidm=all_sidm_matches[0],
                matched_on=identifier_type,
                source="cellosaurus",
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
            )
        return CellLineResolutionResult(
            input_raw=raw,
            input_normalized=normalized,
            detected_type=identifier_type,
            status="unresolved",
            source="cellosaurus",
        )

    def extract_bridge_identifiers(self, identifier: object) -> List[str]:
        """Extract possible bridge IDs (ACH/CVCL) from Cellosaurus payloads."""
        normalized = normalize_identifier(identifier)
        identifier_type = detect_identifier_type(normalized)
        candidates = []
        for payload in self._query_payloads(normalized, identifier_type):
            ach_values = _extract_matches(payload, ACH_SCAN_PATTERN)
            cvcl_values = _extract_matches(payload, CVCL_SCAN_PATTERN)
            for token in ach_values + cvcl_values:
                if token not in candidates:
                    candidates.append(token)
        return candidates


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
    ):
        self.sanger = SangerApiResolver(
            http_client=http_client,
            base_urls=sanger_base_urls,
            timeout=timeout,
        )
        self.cellosaurus = CellosaurusApiResolver(
            http_client=http_client,
            base_url=cellosaurus_base_url,
            timeout=timeout,
        )
        self.use_cellosaurus_fallback = use_cellosaurus_fallback
        self._cache: Dict[str, CellLineResolutionResult] = {} if enable_cache else None

    def resolve_one(self, identifier: object) -> CellLineResolutionResult:
        raw = "" if identifier is None else str(identifier)
        normalized = normalize_identifier(raw)
        
        # Check cache first
        if self._cache is not None and normalized in self._cache:
            return self._cache[normalized]
        
        identifier_type = detect_identifier_type(normalized)

        if not normalized:
            result = CellLineResolutionResult(raw, normalized, identifier_type, "unresolved", source="sanger")
            if self._cache is not None:
                self._cache[normalized] = result
            return result

        primary = self.sanger.resolve_one(identifier)
        if primary.status == "resolved" or primary.status == "ambiguous":
            if self._cache is not None:
                self._cache[normalized] = primary
            return primary

        if not self.use_cellosaurus_fallback:
            if self._cache is not None:
                self._cache[normalized] = primary
            return primary

        fallback = self.cellosaurus.resolve_one(identifier)
        if fallback.status == "resolved" or fallback.status == "ambiguous":
            if self._cache is not None:
                self._cache[normalized] = fallback
            return fallback

        # Try bridge identifiers discovered from Cellosaurus (e.g., ACH/CVCL).
        for token in self.cellosaurus.extract_bridge_identifiers(identifier):
            bridged = self.sanger.resolve_one(token)
            if bridged.status == "resolved":
                bridged.source = "cellosaurus->sanger"
                if self._cache is not None:
                    self._cache[normalized] = bridged
                return bridged

        if self._cache is not None:
            self._cache[normalized] = primary
        return primary
    
    def clear_cache(self) -> None:
        """Clear the in-memory resolution cache."""
        if self._cache is not None:
            self._cache.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """Return cache statistics (e.g., for logging)."""
        if self._cache is None:
            return {"cache_enabled": False}
        return {"cache_enabled": True, "cached_identifiers": len(self._cache)}


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
        result = resolver.resolve_one(value)
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
        result = resolver.resolve_one(normalized_model_id)
        if result.status == "resolved" and result.sidm is not None:
            model_to_sidm[normalized_model_id] = result.sidm
        elif normalized_model_id not in seen_missing:
            not_found.append(normalized_model_id)
            seen_missing.add(normalized_model_id)

    return model_to_sidm, not_found


def resolve_sidm_from_dataframe(
    df: pd.DataFrame, model_registry: Optional[str] = None, verbose: bool = False
) -> Tuple[Dict[str, str], List[str], Dict[str, int]]:
    """Resolve SIDMs row-wise from a user-provided cell-line table.

    Resolution strategy per row:
    1) Try best available identifier column by priority
    2) If unresolved, try display name columns as fallback

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
    alias_to_sidm: Dict[str, str] = {}

    for row_idx, row in df.iterrows():
        row_identifier = None
        for col in available_id_cols:
            value = row[col]
            if pd.notna(value) and str(value).strip():
                row_identifier = str(value).strip()
                break

        if row_identifier is None:
            resolution_counts['unresolved'] += 1
            if verbose:
                print(f"[SIDM][row={row_idx}] unresolved: empty identifier row")
            continue

        display_name = row_identifier
        for col in display_columns_priority:
            if col not in df.columns:
                continue
            value = row[col]
            if pd.notna(value) and str(value).strip():
                display_name = str(value).strip()
                break

        result = resolver.resolve_one(row_identifier)
        if result.status != "resolved":
            # Fallback to display name if it differs from identifier.
            if display_name != row_identifier:
                result = resolver.resolve_one(display_name)

        # record detailed result
        resolution_results.append(result)

        if result.status == "resolved" and result.sidm is not None:
            sidm = result.sidm
            sidm_dict[sidm] = display_name
            resolution_counts['resolved'] += 1
            if verbose:
                print(
                    f"[SIDM][row={row_idx}] resolved '{row_identifier}' -> {sidm} "
                    f"(source={result.source}, matched_on={result.matched_on})"
                )
            # add canonicalized aliases for this resolved row
            if result.input_normalized:
                alias_to_sidm[result.input_normalized] = sidm
            # also map the raw input and the chosen display name
            alias_to_sidm[result.input_raw] = sidm
            if display_name and display_name != result.input_raw:
                alias_to_sidm[display_name] = sidm
        elif result.status == "ambiguous":
            resolution_counts['ambiguous'] += 1
            not_found.append(row_identifier)
            if verbose:
                print(
                    f"[SIDM][row={row_idx}] ambiguous '{row_identifier}' "
                    f"candidates={result.candidate_sidms}"
                )
        else:
            resolution_counts['unresolved'] += 1
            not_found.append(row_identifier)
            if verbose:
                print(f"[SIDM][row={row_idx}] unresolved '{row_identifier}'")

    # Add cache statistics
    resolution_counts['cache_stats'] = resolver.get_cache_stats()
    # expose alias mapping and detailed results to callers (training will use alias map)
    resolution_counts['alias_to_sidm'] = alias_to_sidm
    resolution_counts['detailed_results'] = resolution_results

    if verbose:
        print("[SIDM] Resolution summary:")
        print(f"  total_rows: {resolution_counts['total_rows']}")
        print(f"  resolved: {resolution_counts['resolved']}")
        print(f"  ambiguous: {resolution_counts['ambiguous']}")
        print(f"  unresolved: {resolution_counts['unresolved']}")
        print(f"  alias_count: {len(alias_to_sidm)}")

    return sidm_dict, not_found, resolution_counts


__all__ = [
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
]
