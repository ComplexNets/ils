"""
Fragment List Selector Utility
Allows switching between different fragment lists (short or long)
"""

from typing import List, Dict
from models import shortList_frag, longList_frag

class FragmentSelector:
    """Manages fragment list selection"""
    
    LIST_OPTIONS = {
        'short': shortList_frag.fragments,
        'long': longList_frag.fragments
    }
    
    @classmethod
    def get_available_lists(cls) -> List[str]:
        """Get names of available fragment lists"""
        return list(cls.LIST_OPTIONS.keys())
    
    @classmethod
    def get_fragments(cls, list_name: str = 'short') -> List[Dict]:
        """
        Get fragments from specified list
        Args:
            list_name: Name of list to use ('short' or 'long')
        Returns:
            List of fragment dictionaries
        Raises:
            ValueError if list_name is invalid
        """
        list_name = list_name.lower()
        if list_name not in cls.LIST_OPTIONS:
            raise ValueError(f"Invalid list name: {list_name}. Available options: {', '.join(cls.get_available_lists())}")
        
        return cls.LIST_OPTIONS[list_name]
    
    @classmethod
    def get_fragment_counts(cls, list_name: str = 'short') -> Dict[str, int]:
        """
        Get count of fragments by type in specified list
        Args:
            list_name: Name of list to use ('short' or 'long')
        Returns:
            Dictionary of counts by fragment type
        """
        fragments = cls.get_fragments(list_name)
        counts = {
            'cation': len([f for f in fragments if f.get('fragment_type') == 'cation']),
            'anion': len([f for f in fragments if f.get('fragment_type') == 'anion']),
            'alkyl_chain': len([f for f in fragments if f.get('fragment_type') == 'alkyl_chain']),
            'total': len(fragments)
        }
        return counts
