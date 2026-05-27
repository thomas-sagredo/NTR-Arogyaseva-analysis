"""
Relationship mappings for the NTR Arogyaseva dataset.
"""

one_to_one_mappings = [
    ("surgery", "surgery_code"),
    ("category_name", "category_code"),
]

unidirectional_mappings = [
    ("surgery", "category_name"),
    ("hosp_name", "hosp_type"),
    ("hosp_name", "hosp_district"),
    ("hosp_name", "hosp_location"),
]