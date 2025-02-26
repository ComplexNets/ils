"""
Solubility calculation module for ionic liquids using Group Contribution Methods.
Primary method: Modified UNIFAC-IL method for solubility prediction.
Reference: https://doi.org/10.1016/j.fluid.2008.02.020
"""

from typing import Dict, Optional, List, Tuple
import math
from utils.utils import get_fragment_properties

# Constants for UNIFAC-IL parameters
T_REF = 298.15  # K (reference temperature)

# Environment factors for different chemical environments
ENVIRONMENT_FACTORS = {
    'hydroxyl': 1.25,
    'carboxyl': 1.30,  # Separate factor for COOH groups
    'aromatic': 1.15,
    'conjugated': 1.10,
    'ether': 1.20,
}

# Maximum number of environment factors to apply
MAX_ENVIRONMENT_FACTORS = 2  # Limit stacking of environment factors

def calculate_group_contribution(fragment: Dict) -> float:
    """
    Calculate the contribution of a molecular fragment to solubility.
    
    This is a simplified group contribution approach based on broad molecular fragments,
    not a full UNIFAC-IL implementation. It provides quick screening estimates but does
    not include the complete parameter matrix and activity coefficients of standard UNIFAC.
    
    Args:
        fragment: Dictionary containing fragment properties
        
    Returns:
        float: Contribution value for the fragment
    """
    base_contribution = 0.0
    environment_multiplier = 1.0
    environment_count = 0
    
    # Handle special cases for hydroxyl vs carboxyl
    if 'COOH' in fragment.get('name', ''):
        if environment_count < MAX_ENVIRONMENT_FACTORS:
            environment_multiplier *= ENVIRONMENT_FACTORS['carboxyl']
            environment_count += 1
    elif 'OH' in fragment.get('name', ''):
        if environment_count < MAX_ENVIRONMENT_FACTORS:
            environment_multiplier *= ENVIRONMENT_FACTORS['hydroxyl']
            environment_count += 1
            
    # Check other environmental factors
    for env, factor in ENVIRONMENT_FACTORS.items():
        if env not in ['hydroxyl', 'carboxyl'] and env in fragment.get('name', '').lower():
            if environment_count < MAX_ENVIRONMENT_FACTORS:
                environment_multiplier *= factor
                environment_count += 1
                
    # More granular handling of alkyl groups
    if 'methyl' in fragment.get('name', '').lower():
        base_contribution += 0.15  # CH3 contribution
    elif 'ethyl' in fragment.get('name', '').lower():
        base_contribution += 0.15 + 0.12  # CH3 + CH2
    elif 'propyl' in fragment.get('name', '').lower():
        base_contribution += 0.15 + 2 * 0.12  # CH3 + 2*CH2
    elif 'butyl' in fragment.get('name', '').lower():
        base_contribution += 0.15 + 3 * 0.12  # CH3 + 3*CH2
    
    # Base contributions for common fragments
    if 'imidazolium' in fragment.get('name', '').lower():
        base_contribution += 0.5
    elif 'pyridinium' in fragment.get('name', '').lower():
        base_contribution += 0.45
    elif 'ammonium' in fragment.get('name', '').lower():
        base_contribution += 0.4
    elif 'bf4' in fragment.get('name', '').lower() or 'tetrafluoroborate' in fragment.get('name', '').lower():
        base_contribution += 0.3
    elif 'pf6' in fragment.get('name', '').lower() or 'hexafluorophosphate' in fragment.get('name', '').lower():
        base_contribution += 0.25
    
    return base_contribution * environment_multiplier

def calculate_solubility(cation: Dict, anion: Dict, alkyl_chain: Dict, functional_group: Dict = None,
                      temperature: float = 298.15, max_solubility: float = 1000.0) -> Optional[float]:
    """
    Calculate water solubility of an ionic liquid in g/L
    
    Args:
        cation: Dictionary containing cation information
        anion: Dictionary containing anion information
        alkyl_chain: Dictionary containing alkyl chain information
        functional_group: Dictionary containing functional group information (optional)
        temperature: Temperature in K (default: 298.15 K)
        max_solubility: Maximum solubility value to cap results (default: 1000 g/L)
        
    Returns:
        Solubility in g/L or None if calculation fails
    """
    # Calculate individual contributions
    cation_contribution = calculate_group_contribution(cation)
    anion_contribution = calculate_group_contribution(anion)
    alkyl_contribution = calculate_group_contribution(alkyl_chain)
    
    if functional_group:
        functional_group_contribution = calculate_group_contribution(functional_group)
    else:
        functional_group_contribution = 0.0
    
    # Temperature correction (simplified Arrhenius-like relationship)
    temp_factor = math.exp(0.05 * (temperature - 298.15) / 298.15)
    
    # Calculate total solubility
    total_contribution = cation_contribution + anion_contribution + alkyl_contribution + functional_group_contribution
    solubility = math.exp(total_contribution) * temp_factor
    
    # Validate result
    if validate_solubility(solubility, max_solubility=max_solubility):
        return solubility
    return None

def validate_solubility(solubility: Optional[float], 
                       min_solubility: float = 0.001,  # g/L
                       max_solubility: float = 1000.0  # g/L
                       ) -> bool:
    """
    Validate if the calculated solubility is within reasonable bounds for ionic liquids.
    Args:
        solubility: Calculated solubility in g/L
        min_solubility: Minimum acceptable solubility (default: 0.001 g/L)
        max_solubility: Maximum acceptable solubility (default: 1000 g/L)
    Returns:
        True if solubility is valid, False otherwise
    """
    if solubility is None:
        print("Validation failed: solubility is None")
        return False
    
    if not isinstance(solubility, (int, float)):
        print(f"Validation failed: solubility is not a number, got {type(solubility)}")
        return False
    
    if math.isnan(solubility) or math.isinf(solubility):
        print(f"Validation failed: solubility is {solubility}")
        return False
    
    if solubility < min_solubility:
        print(f"Validation failed: solubility {solubility} g/L is below minimum {min_solubility} g/L")
        return False
    
    if solubility > max_solubility:
        print(f"Validation failed: solubility {solubility} g/L is above maximum {max_solubility} g/L")
        return False
    
    print(f"Validation passed: solubility {solubility} g/L is within bounds [{min_solubility}, {max_solubility}]")
    return True

def screen_fragments_by_solubility(fragments_data: Dict[str, List[Dict]], 
                                 target_range: Tuple[float, float]
                                 ) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated solubility
    Args:
        fragments_data: Dictionary of fragments organized by type
        target_range: Tuple of (min_solubility, max_solubility) in g/L
    Returns:
        Dictionary of fragments organized by type that meet the target range
    """
    min_solubility, max_solubility = target_range
    screened_fragments = {
        'cation': [],
        'anion': [],
        'alkyl_chain': []
    }
    
    # For each fragment type
    for frag_type in fragments_data:
        for fragment in fragments_data[frag_type]:
            # Create a test combination with this fragment
            test_combo = {
                'cation': fragment if frag_type == 'cation' else fragments_data['cation'][0],
                'anion': fragment if frag_type == 'anion' else fragments_data['anion'][0],
                'alkyl_chain': fragment if frag_type == 'alkyl_chain' else fragments_data['alkyl_chain'][0]
            }
            
            # Calculate solubility
            solubility = calculate_solubility(
                test_combo['cation'],
                test_combo['anion'],
                test_combo['alkyl_chain']
            )
            
            # Check if within range
            if solubility and min_solubility <= solubility <= max_solubility:
                screened_fragments[frag_type].append(fragment)
    
    return screened_fragments

def test_solubility_calculations(fragments_data: Optional[Dict[str, List[Dict]]] = None,
                               combinations: Optional[List[Dict]] = None,
                               num_test_combinations: int = 3):
    """
    Test function to show detailed solubility calculations
    Args:
        fragments_data: Optional dictionary of fragments to test screening
        combinations: Optional list of ionic liquid combinations to test
        num_test_combinations: Number of combinations to test
    """
    if combinations:
        print("\nTesting specific combinations:")
        for combo in combinations[:num_test_combinations]:
            solubility = calculate_solubility(
                combo['cation'],
                combo['anion'],
                combo['alkyl_chain']
            )
            print(f"\nCombination:")
            print(f"  Cation: {combo['cation']['name']}")
            print(f"  Anion: {combo['anion']['name']}")
            print(f"  Alkyl chain: {combo['alkyl_chain']['name']}")
            print(f"  Calculated solubility: {solubility:.2f} g/L")

if __name__ == "__main__":
    print("\n=== Solubility Calculation Module Test ===")
    test_solubility_calculations()
