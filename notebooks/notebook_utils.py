"""Utilities for notebooks: locate and set the repository root cleanly.

Expose `find_repo_root` and `ensure_repo_root` so notebooks can remain terse.
"""

import sys
import os
from pathlib import Path


def find_repo_root(markers=(".git", "LICENSE", "requirements.txt")) -> Path:
    """Find the repository root by looking for common markers."""
    current_path = Path(__file__).resolve()
    for parent in current_path.parents:
        if any((parent / marker).exists() for marker in markers):
            return parent
    raise FileNotFoundError("Could not find repository root (no marker found).")


def ensure_repo_root(markers=(".git", "LICENSE", "requirements.txt")) -> Path:
    """Find the repository root and change the current working directory to it."""
    root = find_repo_root(markers)
    sys.path.insert(0, str(root))
    os.chdir(str(root))
    return root
