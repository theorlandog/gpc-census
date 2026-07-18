"""gpc-census: exact extremal states for fermionic natural-occupation-number polytopes."""

from importlib.metadata import PackageNotFoundError, version

from gpc_census.core import slater_vertices

__all__ = ["__version__", "slater_vertices"]

try:
    __version__ = version("gpc-census")
except PackageNotFoundError:  # running from a checkout without installation
    __version__ = "0+unknown"
