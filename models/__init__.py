"""
Models for ILS
"""
from .shortList_frag import fragments as short_fragments
from .mediumList_frag import fragments as medium_fragments
from .longList_frag import fragments as long_fragments

__all__ = ['short_fragments', 'medium_fragments', 'long_fragments']