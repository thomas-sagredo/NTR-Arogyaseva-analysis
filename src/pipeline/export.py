"""Export public-ready datasets."""

from pathlib import Path

from src.utils.io import export_parquet

from src.config.constants import CLEAN_DATA_DIR


def export_dataset(df):
    output_path = Path(CLEAN_DATA_DIR) / "ntrarogyaseva.parquet"
    export_parquet(df, output_path)