# combine_fragHeatCapacity.py - functions to calculate heat capacity of ILs
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from core.combine_fragments import combine_fragments

# UNIFAC group parameters for heat capacity calculation
unifac_params = {
    'IM+': {'A': 150.0, 'B': 0.15, 'C': -0.0002},  # Imidazolium
    'PY+': {'A': 145.0, 'B': 0.12, 'C': -0.0002},  # Pyridinium
    'BF4-': {'A': 130.0, 'B': 0.10, 'C': -0.0001}, # Tetrafluoroborate
    'PF6-': {'A': 140.0, 'B': 0.11, 'C': -0.0001}, # Hexafluorophosphate
    'SCN-': {'A': 120.0, 'B': 0.09, 'C': -0.0001}, # Thiocyanate
    'DCA-': {'A': 125.0, 'B': 0.10, 'C': -0.0001}, # Dicyanamide
    'CH2': {'A': 30.0, 'B': 0.05, 'C': -0.0001},   # Methylene group
    'CH3': {'A': 35.0, 'B': 0.05, 'C': -0.0001}    # Methyl group
}

def get_cation_type(cation_name: str) -> str:
    """Map cation names to their UNIFAC types"""
    if 'imidazolium' in cation_name.lower():
        return 'IM+'
    elif 'pyridinium' in cation_name.lower():
        return 'PY+'
    return 'IM+'  # Default to imidazolium if unknown

def get_anion_type(anion_name: str) -> str:
    """Map anion names to their UNIFAC types"""
    name_lower = anion_name.lower()
    if 'tetrafluoroborate' in name_lower or 'bf4' in name_lower:
        return 'BF4-'
    elif 'hexafluorophosphate' in name_lower or 'pf6' in name_lower:
        return 'PF6-'
    elif 'thiocyanate' in name_lower or 'scn' in name_lower:
        return 'SCN-'
    elif 'dicyanamide' in name_lower or 'dca' in name_lower:
        return 'DCA-'
    return 'BF4-'  # Default to BF4 if unknown

def add_alkyl_contributions(chain_length: int, T: float) -> float:
    """Calculate heat capacity contribution from alkyl chain"""
    ch2_count = max(0, chain_length - 1)  # Number of CH2 groups
    ch3_count = 1  # Always one CH3 terminal group
    
    ch2_contrib = ch2_count * unifac_contribution(unifac_params, 'CH2', T)
    ch3_contrib = ch3_count * unifac_contribution(unifac_params, 'CH3', T)
    
    return ch2_contrib + ch3_contrib

def unifac_contribution(params: Dict, group_type: str, T: float) -> float:
    """Calculate UNIFAC contribution with temperature dependence"""
    if group_type not in params:
        return 0.0
    dT = T - 298.15
    return params[group_type]['A'] + params[group_type]['B'] * dT + params[group_type]['C'] * dT**2

def calculate_heat_capacity(ionic_liquid: Dict, T: float = 298.15) -> Optional[float]:
    """
    Calculate heat capacity for a given ionic liquid combination using UNIFAC method
    Args:
        ionic_liquid: Dictionary containing cation, anion, and alkyl chain information
        T: Temperature in Kelvin (default: 298.15 K)
    Returns:
        Heat capacity in J/mol·K
    """
    try:
        # Get fragment types
        cation_type = get_cation_type(ionic_liquid['cation']['Name'])
        anion_type = get_anion_type(ionic_liquid['anion']['Name'])
        
        # Calculate contributions from each part
        cation_contrib = unifac_contribution(unifac_params, cation_type, T)
        anion_contrib = unifac_contribution(unifac_params, anion_type, T)
        
        # Estimate chain length from molecular weight
        alkyl_mw = ionic_liquid['alkyl']['MolecularWeight']
        chain_length = max(1, round((alkyl_mw - 15) / 14))  # Approximate CH2 count
        alkyl_contrib = add_alkyl_contributions(chain_length, T)
        
        # Sum all contributions
        total_cp = cation_contrib + anion_contrib + alkyl_contrib
        
        # Add corrections based on molecular properties
        total_cp *= (1 + 0.02 * ionic_liquid['cation']['HydrogenBondDonorCount'])
        total_cp *= (1 + 0.01 * ionic_liquid['anion']['HydrogenBondAcceptorCount'])
        
        return max(200, min(400, total_cp))  # Constrain to reasonable range
        
    except Exception as e:
        print(f"Error calculating heat capacity: {str(e)}")
        return None

def calculate_all_heat_capacities(combinations: List[Dict]) -> List[Dict]:
    """Calculate heat capacities for a list of ionic liquid combinations"""
    for combination in combinations:
        cp = calculate_heat_capacity(combination)
        if cp is not None:
            combination['cp_unifac'] = cp
    return combinations

def main():
    """Main function to get combinations and calculate their heat capacities"""
    try:
        # Get valid combinations
        combinations = combine_fragments()
        
        # Calculate heat capacities
        combinations_with_cp = calculate_all_heat_capacities(combinations)
        return combinations_with_cp
        
    except Exception as e:
        print(f"Error in heat capacity calculation main: {str(e)}")
        return []

if __name__ == "__main__":
    print("\n=== Running Ionic Liquid Heat Capacity Calculator ===\n")
    results = main()
    print(f"\nCalculated heat capacities for {len(results)} combinations")
