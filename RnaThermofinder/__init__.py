"""
RNA Thermometer Finder
A tool for identifying and analyzing RNA thermometers in bacterial sequences
"""

__version__ = "1.0.0"
__author__ = "Roy Vaknin"
__email__ = "roycyber13@gmail.com"

# Make key modules easily accessible
from .core import FastaParse
from .core import HairpinAnalysis

__all__ = [
    'FastaParse',
    'HairpinAnalysis',
    '__version__'
]