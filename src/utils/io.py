"""Reusable input/output utilities."""

from pathlib import Path

import pandas as pd

from src.config.constants import (
    RAW_DATA_DIR,
    CLEAN_DATA_DIR,
    PANDAS_DISPLAY_OPTIONS,
)


def read_csv_safe(path: str | Path) -> pd.DataFrame | None:
    try:
        path = Path(path)
        if not path.exists():
            print(f"File not found: {path}")
            return None
        return pd.read_csv(path)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return None
    

def read_parquet_safe(path: str | Path) -> pd.DataFrame | None:
    try:
        path = Path(path)
        if not path.exists():
            print(f"File not found: {path}")
            return None
        return pd.read_parquet(path)
    except Exception as e:
        print(f"Error reading Parquet file: {e}")
        return None


def export_parquet(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    df.to_parquet(path, index=False)


def set_pandas_display_options():
    """Set pandas display options for consistent formatting across the notebook."""
    for option, value in PANDAS_DISPLAY_OPTIONS.items():
        pd.set_option(option, value)


def load_data(
        stage: str = "clean"
) -> pd.DataFrame:
    """
    Load the dataset from the specified data directory.
    """
    if stage.lower().strip() not in ["raw", "clean"]:
        raise ValueError(f"Invalid stage: {stage}. Must be 'raw' or 'clean'.")
    
    if stage.lower().strip() == "clean":
        path = Path(CLEAN_DATA_DIR) / "ntrarogyaseva.parquet"

        return pd.read_parquet(path)
    
    else:
        path = Path(RAW_DATA_DIR) / "ntrarogyaseva.csv"
        
        return pd.read_csv(path)