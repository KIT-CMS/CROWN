from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import Any

# File with helper functions for the CROWN code generation


def is_empty(value: Any) -> bool:
    """
    Check if a value is empty.

    Args:
        value: The value that should be checked.

    Returns:
        bool: Whether the input value is considered 'empty'
    """
    empty_values = [None]

    try:
        length = len(value)
    except TypeError:
        length = -1

    return value in empty_values or length == 0
