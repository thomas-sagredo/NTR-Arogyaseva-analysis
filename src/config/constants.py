"""Global constants used across the pipeline."""

# Directory paths for data
RAW_DATA_DIR = "data/raw"
CLEAN_DATA_DIR = "data/clean"

# Display options for pandas DataFrames to ensure consistent formatting across the pipeline.
PANDAS_DISPLAY_OPTIONS = {
    "display.max_columns": 50,
    "display.max_rows": 200,
    "display.width": 140,
    "display.float_format": "{:.4f}".format,
    "display.max_colwidth": 80,
}