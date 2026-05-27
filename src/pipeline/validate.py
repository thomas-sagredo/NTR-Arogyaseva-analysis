"""Validation checks for public datasets."""

import pandas as pd

from src.config.mappings import (
    one_to_one_mappings, 
    unidirectional_mappings)


def validate_row_count(raw_df: pd.DataFrame, public_df: pd.DataFrame) -> None:
    assert len(raw_df) == len(public_df), (
        "Row count mismatch between raw and public datasets."
    )


def flag_one_to_one_inconsistencies(
    df: pd.DataFrame,
    first_col: str,
    second_col: str,
    check_both_directions: bool = True,
    first_to_second_col: str | None = None,
    second_to_first_col: str | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Flag rows that violate a one-to-one relationship between two columns.

    When check_both_directions is False, only first_col -> second_col is checked.
    """
    df = df.copy()

    first_to_second_col = (
        first_to_second_col
        or f"inconsistent_{first_col}_to_{second_col}"
    )

    first_to_second = (
        df.groupby(first_col, dropna=False)[second_col]
        .transform(lambda values: values.nunique(dropna=False))
        > 1
    )
    if first_to_second.any():
        df[first_to_second_col] = first_to_second.astype(bool)

        if verbose:
            print(
                f"Found {first_to_second.sum()} rows where {first_col} maps to multiple {second_col} values."
            )
    else:
        if verbose:
            print(f"No inconsistencies found for {first_col} -> {second_col}.")
            print("The inconsistency flag column will not be added to the DataFrame.")

    if check_both_directions:
        df = flag_one_to_one_inconsistencies(
            df=df,
            first_col=second_col,
            second_col=first_col,
            check_both_directions=False,
            first_to_second_col=second_to_first_col,
            verbose=verbose,
        )

    return df


def validate_one_to_one_relationships(df: pd.DataFrame) -> pd.DataFrame:
    """Validate one-to-one relationships defined in the mappings."""
    for first_col, second_col in one_to_one_mappings:
        df = flag_one_to_one_inconsistencies(
            df=df,
            first_col=first_col,
            second_col=second_col,
            check_both_directions=True,
            verbose=True,
        )
    for first_col, second_col in unidirectional_mappings:
        df = flag_one_to_one_inconsistencies(
            df=df,
            first_col=first_col,
            second_col=second_col,
            check_both_directions=False,
            verbose=True,
        )
    return df


def flag_mortality_without_date(
    df: pd.DataFrame,
    flag_col: str = "mortality_without_date",
    verbose: bool = False,
) -> pd.DataFrame:
    """Flag rows where mortality is Y and mortality_date is null."""
    df = df.copy()

    mortality_yes = (
        df["mortality y / n"]
        .astype("string")
        .str.strip()
        .str.upper()
        .eq("Y")
    )
    mortality_date_is_null = df["mortality_date"].isna()

    inconsistent = mortality_yes & mortality_date_is_null

    if inconsistent.any():
        df[flag_col] = inconsistent

        if verbose:
            print(
                f"Column {flag_col} created with {inconsistent.sum()} flagged rows."
            )
    else:
        if verbose:
            print(
                "No rows found where mortality y / n is Y and mortality_date is null."
            )
            print("The inconsistency flag column will not be added to the DataFrame.")

    return df


