"""Load raw dataset."""

from pathlib import Path

from src.utils.io import read_csv_safe
    
from src.config.constants import RAW_DATA_DIR

def load_dataset():
    path = Path(RAW_DATA_DIR) / "ntrarogyaseva.csv"
    return read_csv_safe(path)