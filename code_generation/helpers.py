from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import Any, Dict, Generator
from contextlib import contextmanager
import contextvars
import inspect
import re

# File with helper functions for the CROWN code generation

# Context registry for managing default values across code generation
CONTEXT_REGISTRY: Dict[str, contextvars.ContextVar] = {
    k: contextvars.ContextVar(k, default=None)
    for k in [
        "name",
        "scopes",
        "shift_key",
        "shift_map",
        "producers",
        "ignore_producers",
        "samples",
        "exclude_samples",
        "input",
        "subproducers",
        "call",
        "output",
        "vec_configs",
    ]
}


@contextmanager
def defaults(**kwargs: Any) -> Generator[None, None, None]:
    """Context manager for setting default values for producer and systematic shift configuration.
    
    Args:
        **kwargs: Default values to set for various parameters
        
    Example:
        with defaults(scopes=['global'], call='myFunction({input})'):
            producer = Producer(...)
            
        with defaults(shift_key='scale', shift_map={'Up': [1.1], 'Down': [0.9]}):
            add_shift(name='jes', producers=[producer])
    """
    tokens = []
    try:
        for key, value in kwargs.items():
            try:
                tokens.append(CONTEXT_REGISTRY[key].set(value))
            except KeyError:
                raise ValueError(f"Unknown context variable: {key}")
        yield
    finally:
        for token in reversed(tokens):
            token.var.reset(token)


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


def get_variable_name() -> str:
    """Automatically determine the variable name from the calling context.
    
    Returns:
        The variable name being assigned to
        
    Raises:
        RuntimeError: If variable name cannot be determined from context
    """
    frame = inspect.currentframe().f_back.f_back
    code_context = inspect.getframeinfo(frame).code_context
    if code_context:
        call_line = code_context[0].strip()
        match = re.match(r"([\w\d_]+)\s*(=|:=)", call_line)
        if match:
            return match.group(1)
    raise RuntimeError("Could not determine variable name from context")


class MissingValue(Exception):
    """Exception raised when a required value is missing."""
    def __init__(self, variable_name: str):
        super().__init__(
            f"Missing value for variable '{variable_name}'. "
            "Please provide a value either as an argument or "
            "through a 'with defaults(...)' context."
        )


class NameNotDetermined(Exception):
    """Exception raised when the name cannot be automatically determined."""
    def __init__(self):
        super().__init__(
            "Name could not be determined. "
            "This should not happened. Workaround: provide the name explicitly"
        )
