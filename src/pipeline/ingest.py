"""Load raw dataset."""

from pathlib import Path

import pandas as pd
    
from src.config.constants import RAW_DATA_DIR

def load_dataset():
    path = Path(RAW_DATA_DIR) / "ntrarogyaseva.csv"
    try:
        if not path.exists():
            print(f"File not found: {path}")
            return None
        return pd.read_csv(path, index_col=[0])
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return None