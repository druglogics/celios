# CELIOS Improvement Roadmap

This document outlines planned improvements for the CELIOS pipeline, prioritized by impact and feasibility.

## High Priority (Foundation)

### 1. Configuration Schema Validation
**Status:** Planned  
**Feasibility:** Medium (1 week)  
**Benefit:** High  
**Description:** Implement Pydantic-based validation for config dictionaries to catch errors early and provide better user feedback.

**Steps:**
- Phase 1: Schema Design (1-2 weeks) - Analyze config structure, define schemas
- Phase 2: Core Validation Logic (2-3 weeks) - Create validation classes, integrate into run_celios
- Phase 3: Update Configs & Docs (1-2 weeks) - Update defaults and documentation
- Phase 4: Testing & Edge Cases (1 week) - Test invalid configs and ensure compatibility

**Attention Points:** Avoid breaking changes; provide migration guide for existing configs.

## Medium Priority (Performance)

### 2. Optimize Memory Usage
**Status:** Planned  
**Feasibility:** Hard (4-8 weeks)  
**Benefit:** High  
**Description:** Replace in-memory DataFrame processing with chunked/lazy loading to handle larger datasets without OOM errors.

**Steps:**
- Phase 1: Assessment & Setup (1-2 weeks) - Profile memory usage, choose library (Dask/Polars)
- Phase 2: Core Refactoring (2-4 weeks) - Refactor matrix building and I/O operations
- Phase 3: Integration & Testing (1-2 weeks) - Integrate with fallback, test on large datasets
- Phase 4: Documentation & Rollout (1 week) - Update docs and add config options

**Attention Points:** Test on limited RAM hardware; ensure backward compatibility.

## Low Priority (Extensibility)

### 3. Implement Plugin Architecture
**Status:** Planned  
**Feasibility:** Hard (6-10 weeks)  
**Benefit:** Medium  
**Description:** Create a plugin system for data sources to allow custom omics data without modifying core code.

**Steps:**
- Phase 1: Design & Prototyping (2-3 weeks) - Define plugin interface, prototype with existing source
- Phase 2: Core Plugin Framework (3-4 weeks) - Create plugin loader and discovery system
- Phase 3: Migrate Existing Sources (2-3 weeks) - Convert hardcoded sources to plugins
- Phase 4: Testing & Documentation (2 weeks) - Create examples and update docs

**Attention Points:** Handle versioning and security; provide clear development guide.

## Implementation Notes

- **Timeline:** Start with Configuration Schema Validation as it's foundational and quickest.
- **Dependencies:** Each improvement may require new libraries (Pydantic, Dask, etc.).
- **Testing:** All changes must maintain backward compatibility and pass existing tests.
- **Documentation:** Update PROJECT_STRUCTURE.md and README.md after each improvement.

## Completed Improvements

- ✅ Tissue-aware output organization (implemented in tissue.py)
- ✅ Modular tissue logic extraction
- ✅ Updated test suite for both legacy and tissue modes
- ✅ Enhanced documentation and project structure

---

**Last Updated:** February 2026