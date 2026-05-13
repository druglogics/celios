"""
Internal database resource management for DREXPA.

Provides unified access to the drug-target interaction database shipped with the package.
"""
from pathlib import Path


def get_internal_database_path() -> str:
    """
    Get absolute path to the internal drug-target interaction database.

    The database is shipped with the package under drexpa/resources/
    and managed by the project team (not user-configurable).

    Returns:
        str: Absolute path to DrugTargetInteractionDB.db

    Raises:
        FileNotFoundError: If the database does not exist in the package resources.
    """
    # Locate the database relative to this module
    db_dir = Path(__file__).parent
    db_path = db_dir / "DrugTargetInteractionDB.db"

    if not db_path.exists():
        raise FileNotFoundError(
            f"Internal drug-target interaction database not found at {db_path}. "
            "This may indicate a broken DREXPA installation. "
            "Please reinstall the package: pip install --upgrade drexpa"
        )

    return str(db_path)
