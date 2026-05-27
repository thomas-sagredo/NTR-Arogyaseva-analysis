"""End-to-end clean dataset builder."""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.pipeline.clean import (
    convert_to_datetime,
    normalize_column_names,
    normalize_sex_column,
    normalize_text_columns,
)
from src.pipeline.export import export_dataset
from src.pipeline.ingest import load_dataset
from src.pipeline.validate import (
    flag_mortality_without_date,
    validate_one_to_one_relationships,
    validate_row_count,
)


def main():
    print("Loading raw dataset...")

    raw_df = load_dataset()
    assert raw_df is not None, "Failed to load raw dataset"

    print("Applying structural cleaning...")

    clean_df = normalize_column_names(raw_df)
    clean_df = convert_to_datetime(clean_df)
    clean_df = normalize_text_columns(clean_df)
    clean_df = normalize_sex_column(clean_df)

    print("Running validation checks...")

    validate_row_count(raw_df, clean_df)
    clean_df = validate_one_to_one_relationships(clean_df)
    clean_df = flag_mortality_without_date(clean_df, verbose=True)

    print("Exporting clean dataset...")

    export_dataset(clean_df)

    print("Clean dataset successfully built.")


if __name__ == "__main__":
    main()
