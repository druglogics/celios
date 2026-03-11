"""Features package: collected pipeline-facing modules.

This package provides a stable namespace `celios.features` so we can
gradually move implementations here without breaking existing imports.
"""

from ..utils import io
from . import node, sifbase, training, tissue

__all__ = ["io", "node", "sifbase", "training", "tissue"]
