# CELIOS

**CEll LIne OmicS processor** — extract omics data into integrated activity datasets for Boolean model calibration.

[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Used in the [DrugLogics](https://github.com/druglogics) and [TRAFIKK](https://github.com/druglogics/trafikk) pipelines.

---

## Quick Start

**Install:**
```bash
pip install -e .
```

**Run pipeline:**
```bash
celios run --config config.yaml --verbose
```

**Get help:**
```bash
celios --help
```

---

## 📚 Documentation

Detailed documentation is in the root-level markdown files:

| Document | Purpose |
|----------|---------|
| **[QUICKSTART.md](QUICKSTART.md)** | 5-minute quick reference with common commands |
| **[INSTALL.md](INSTALL.md)** | Installation guide, virtual environments, troubleshooting |
| **[PIPELINE.md](PIPELINE.md)** | Detailed 3-step pipeline overview with configurations |
| **[USAGE.md](USAGE.md)** | Full usage guide, advanced examples, and API reference |
| **[CONFIGURATION.md](CONFIGURATION.md)** | Configuration reference with all options |
| **[OUTPUTS.md](OUTPUTS.md)** | Output file formats and interpretation |
| **[ARCHITECTURE.md](ARCHITECTURE.md)** | Code organization, module descriptions, architecture |
| **[notebooks/](notebooks/)** | Interactive Jupyter notebooks with examples |

---

## Features

- **Step 1:** Node dictionary generation from biological networks (SIF format)
- **Step 2:** Cell-line identifier resolution and tissue-aware organization
- **Step 3:** Multi-source omics integration (mutations, CNV, TF activity, expression)
- **Format support:** Legacy CCLE (genes × SIDM) and 26Q1 (ModelID × genes) formats
- **Configuration-driven:** JSON or YAML configs for easy reproducibility
- **Tissue organization:** Optional per-tissue output structure for DrugLogics compatibility

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) (if available) or submit issues to GitHub.

---

## License

[MIT License](LICENSE)

