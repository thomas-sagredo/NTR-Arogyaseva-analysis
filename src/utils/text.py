"""Generic text normalization utilities."""

import re
import unicodedata


def remove_accents(text: str) -> str:
    normalized = unicodedata.normalize("NFKD", text)
    return "".join(char for char in normalized if not unicodedata.combining(char))


def normalize_text(text: str) -> str:
    if text is None:
        return text

    text = str(text).strip().lower()
    text = remove_accents(text)
    text = re.sub(r"\s+", " ", text)

    return text


def to_snake_case(text: str) -> str:
    text = normalize_text(text)
    text = re.sub(r"[^a-z0-9]+", "_", text)

    return text.strip("_")