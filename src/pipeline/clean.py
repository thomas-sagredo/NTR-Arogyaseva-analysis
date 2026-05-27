"""Structural cleaning logic."""

import pandas as pd

from src.utils.text import normalize_text


def normalize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names to a consistent format."""
    df = df.copy()
    df = df.rename(columns=lambda x: normalize_text(x))
    return df


def convert_to_datetime(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert specified columns to datetime format.
    """
    df = df.copy()

    date_cols = df.select_dtypes(include=["object"]).columns
    date_cols = date_cols[date_cols.str.endswith("_date")]

    for col in date_cols:
        try:
            df[col] = pd.to_datetime(df[col], format="%d/%m/%Y %H:%M:%S", dayfirst=True, errors="coerce")
            print(f"Column `{col}` successfully converted to datetime.")
        except Exception as e:
            print(f"Error converting column `{col}`: {e}")

    return df


def normalize_text_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize text columns in the dataframe using the provided normalization function."""

    text_cols = df.select_dtypes(include=["object"]).columns

    df = df.copy()

    df[text_cols] = df[text_cols].map(normalize_text)

    print(f"Text columns normalized: {text_cols}")    

    return df


def normalize_sex_column(df: pd.DataFrame) -> pd.DataFrame:
    """Map sex values to male/female and flag child rows by age."""
    df = df.copy()

    sex = df["sex"].astype("string").str.strip().str.lower()
    age = pd.to_numeric(df["age"], errors="coerce") # Convert age to numeric, coercing errors to NaN
    is_female = sex.str.contains("female", na=False)
    is_male = sex.str.contains("male", na=False) & ~is_female

    df["is_child"] = age <= 14
    df.loc[is_female, "sex"] = "female"
    df.loc[is_male, "sex"] = "male"

    print("Sex column normalized and is_child column created.")

    return df
