import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from core.combine_fragments import get_filtered_fragments, combine_fragments
from utils.utils import get_fragment_properties

# UNIFAC group parameters for heat capacity calculation
# Values derived from literature data for common ionic liquids
unifac_params = {
    'IM+': {'A': 142.0, 'B': 0.15, 'C': -0.0002},  # Imidazolium core contribution
    'PY+': {'A': 138.0, 'B': 0.14, 'C': -0.0002},  # Pyridinium (scaled relative to IM+)
    'BF4-': {'A': 120.0, 'B': 0.12, 'C': -0.0001}, # From [EMIM][BF4] and [BMIM][BF4] data
    'PF6-': {'A': 130.0, 'B': 0.13, 'C': -0.0001}, # Literature data for PF6-based ILs
    'NTf2-': {'A': 140.0, 'B': 0.12, 'C': -0.0001}, # Reduced from 180.0 to better match literature
    'SCN-': {'A': 95.0, 'B': 0.10, 'C': -0.0001},  # Literature data for thiocyanate ILs
    'DCA-': {'A': 100.0, 'B': 0.11, 'C': -0.0001}, # Literature data for dicyanamide ILs
    'CH2': {'A': 31.0, 'B': 0.05, 'C': -0.0001},   # Derived from alkyl chain length variations
    'CH3': {'A': 35.0, 'B': 0.05, 'C': -0.0001}    # Terminal methyl contribution
}

def estimate_fragment_heat_capacity(fragment: Dict) -> Optional[float]:
    """
    Estimate heat capacity for a single fragment
    Based on literature values for common ionic liquids:
    - [EMIM][BF4]: ~320 J/mol·K at 298.15K
    - [BMIM][PF6]: ~350 J/mol·K at 298.15K
    """
    try:
        # Get properties from database
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props:
            print(f"  No properties found for {fragment['name']}")
            return None
            
        # Get required properties
        mw = props['molecular_weight']
        heavy_atoms = props['heavy_atom_count']
        rotatable_bonds = props['rotatable_bond_count']
        h_donors = props['hydrogen_bond_donor_count']
        h_acceptors = props['hydrogen_bond_acceptor_count']
        
        if not all([mw, heavy_atoms is not None]):
            print(f"  Missing required properties for {fragment['name']}")
            return None
            
        # Base heat capacity estimation based on literature correlations
        # 1. Molecular weight contribution (from group contribution methods)
        mw_contribution = 0.92 * mw  # Reduced from 0.95 to better match NTf2 ILs
        
        # 2. Heavy atom contribution (based on atomic group additivity)
        atom_contribution = 4.8 * heavy_atoms  # Reduced from 5.2 for better NTf2 prediction
        
        # 3. Rotatable bond contribution (from conformational analysis)
        flexibility_contribution = 2.8 * rotatable_bonds
        
        # 4. Hydrogen bonding contribution (from H-bond network studies)
        h_bond_contribution = 1.8 * h_donors + 0.9 * h_acceptors
        
        # Fragment-specific adjustment factors from ionic liquid literature
        if fragment['fragment_type'] == 'Cation':
            base_adjustment = 1.15
        elif fragment['fragment_type'] == 'Anion':
            # Apply additional scaling for large anions like NTf2
            if fragment['name'] in ['bis(trifluoromethanesulfonyl)imide', 'NTf2']:
                base_adjustment = 0.95  # Reduced adjustment for NTf2
            else:
                base_adjustment = 1.05
        else:  # Alkyl Chain
            base_adjustment = 0.95
            
        # Calculate total heat capacity
        total_cp = (
            mw_contribution +
            atom_contribution +
            flexibility_contribution +
            h_bond_contribution
        ) * base_adjustment
        
        # Add UNIFAC contribution with literature-based scaling
        unifac_type = get_fragment_type(fragment['name'], fragment['fragment_type'])
        if unifac_type:
            unifac_cp = unifac_contribution(unifac_type)
            # Additional scaling for NTf2-containing ILs
            if unifac_type == 'NTf2-' or fragment.get('name') in ['bis(trifluoromethanesulfonyl)imide', 'NTf2']:
                total_cp = (total_cp + unifac_cp) * 0.85
            else:
                total_cp = (total_cp + unifac_cp) * 0.95
        
        # Print detailed contributions for debugging
        print(f"  Heat capacity contributions for {fragment['name']}:")
        print(f"    MW contribution: {mw_contribution:.1f}")
        print(f"    Atom contribution: {atom_contribution:.1f}")
        print(f"    Flexibility contribution: {flexibility_contribution:.1f}")
        print(f"    H-bond contribution: {h_bond_contribution:.1f}")
        print(f"    Base adjustment: {base_adjustment:.2f}")
        print(f"    Total: {total_cp:.1f}")
        
        return total_cp
        
    except Exception as e:
        print(f"Error estimating fragment heat capacity: {e}")
        return None

def get_fragment_type(fragment_name: str, fragment_type: str) -> str:
    """Map fragment names to their UNIFAC types"""
    name_lower = fragment_name.lower()
    
    if fragment_type == 'Cation':
        if 'imidazolium' in name_lower:
            return 'IM+'
        elif 'pyridinium' in name_lower:
            return 'PY+'
        return 'IM+'  # Default to imidazolium if unknown
        
    elif fragment_type == 'Anion':
        if 'tetrafluoroborate' in name_lower or 'bf4' in name_lower:
            return 'BF4-'
        elif 'hexafluorophosphate' in name_lower or 'pf6' in name_lower:
            return 'PF6-'
        elif 'thiocyanate' in name_lower or 'scn' in name_lower:
            return 'SCN-'
        elif 'dicyanamide' in name_lower or 'dca' in name_lower:
            return 'DCA-'
        return 'BF4-'  # Default to BF4 if unknown
    
    return None

def unifac_contribution(group_type: str, T: float = 298.15) -> float:
    """Calculate UNIFAC contribution with temperature dependence"""
    if group_type not in unifac_params:
        return 0.0
    dT = T - 298.15
    return unifac_params[group_type]['A'] + unifac_params[group_type]['B'] * dT + unifac_params[group_type]['C'] * dT**2

def add_alkyl_contributions(chain_length: int, T: float = 298.15) -> float:
    """Calculate heat capacity contribution from alkyl chain"""
    ch2_count = max(0, chain_length - 1)  # Number of CH2 groups
    ch3_count = 1  # Always one CH3 terminal group
    
    ch2_contrib = ch2_count * unifac_contribution('CH2', T)
    ch3_contrib = ch3_count * unifac_contribution('CH3', T)
    
    return ch2_contrib + ch3_contrib

def calculate_ionic_liquid_heat_capacity(ionic_liquid: Dict) -> Optional[float]:
    """
    Calculate heat capacity for a complete ionic liquid
    Args:
        ionic_liquid: Dictionary containing cation, anion, and alkyl chain information
    Returns:
        Heat capacity in J/mol·K
    """
    try:
        print(f"\nCalculating heat capacity for {ionic_liquid['cation'].get('name', 'Unknown')} {ionic_liquid['anion'].get('name', 'Unknown')} with {ionic_liquid['alkyl_chain'].get('name', 'Unknown')}:")
        
        # Calculate heat capacity for each component
        cation_cp = estimate_fragment_heat_capacity(ionic_liquid['cation'])
        anion_cp = estimate_fragment_heat_capacity(ionic_liquid['anion'])
        alkyl_cp = estimate_fragment_heat_capacity(ionic_liquid['alkyl_chain'])
        
        print(f"\nComponent heat capacities:")
        print(f"  Cation: {f'{cation_cp:.1f}' if cation_cp is not None else 'None'} J/mol·K")
        print(f"  Anion: {f'{anion_cp:.1f}' if anion_cp is not None else 'None'} J/mol·K")
        print(f"  Alkyl: {f'{alkyl_cp:.1f}' if alkyl_cp is not None else 'None'} J/mol·K")
        
        if any(cp is None for cp in [cation_cp, anion_cp, alkyl_cp]):
            print("Error: One or more heat capacities are None")
            return None
            
        # Total heat capacity is sum of fragments plus interaction terms
        total_cp = cation_cp + anion_cp + alkyl_cp
        
        # Add interaction terms based on hydrogen bonding
        cation_h_donors = ionic_liquid['cation'].get('h_bond_donors', 0)
        cation_h_acceptors = ionic_liquid['cation'].get('h_bond_acceptors', 0)
        anion_h_donors = ionic_liquid['anion'].get('h_bond_donors', 0)
        anion_h_acceptors = ionic_liquid['anion'].get('h_bond_acceptors', 0)
        
        # Each H-bond pair contributes ~5 J/mol·K
        h_bond_pairs = min(cation_h_donors, anion_h_acceptors) + min(anion_h_donors, cation_h_acceptors)
        h_bond_contribution = 5.0 * h_bond_pairs
        
        print(f"  H-bond pairs: {h_bond_pairs}")
        print(f"  H-bond contribution: {h_bond_contribution:.1f} J/mol·K")
        
        total_cp += h_bond_contribution
        
        print(f"  Total heat capacity: {total_cp:.1f} J/mol·K")
        return total_cp
        
    except Exception as e:
        print(f"Error calculating ionic liquid heat capacity: {e}")
        return None

def screen_fragments_by_heat_capacity(fragments_data: Dict[str, List[Dict]], target_range: tuple) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated heat capacity
    Args:
        fragments_data: Dictionary of fragments organized by type (cation, anion, alkyl_chain)
        target_range: Tuple of (min_capacity, max_capacity) in J/mol·K
    Returns:
        Dictionary of fragments organized by type that meet the target range
    """
    min_capacity, max_capacity = target_range
    
    # Initialize output dictionary with same structure as input
    screened_fragments = {frag_type: [] for frag_type in fragments_data.keys()}
    
    for frag_type, fragments in fragments_data.items():
        for fragment in fragments:
            estimated_capacity = estimate_fragment_heat_capacity(fragment)
            
            if estimated_capacity is not None:
                # Use configurable scaling factors for component contribution
                scaling_factor = 1/3  # Can be adjusted or passed as parameter
                if min_capacity * scaling_factor <= estimated_capacity <= max_capacity * scaling_factor:
                    screened_fragments[frag_type].append({
                        **fragment,  # Preserve all original properties
                        'estimated_heat_capacity': estimated_capacity
                    })
    
    return screened_fragments

def test_heat_capacity_calculations(fragments_data: Optional[Dict[str, List[Dict]]] = None,
                                  combinations: Optional[List[Dict]] = None,
                                  num_test_combinations: int = 3):
    """
    Test function to show detailed heat capacity calculations
    Args:
        fragments_data: Optional dictionary of fragments to test screening
        combinations: Optional list of ionic liquid combinations to test
        num_test_combinations: Number of combinations to test
    """
    print("\n=== Testing Heat Capacity Calculations ===")
    
    if combinations:
        # Test provided combinations
        for i, il in enumerate(combinations[:num_test_combinations]):
            print(f"\nTesting combination {i+1}:")
            print(f"Cation: {il['cation']['name']}")
            print(f"Anion: {il['anion']['name']}")
            if 'alkyl_chain' in il:
                print(f"Alkyl chain: {il['alkyl_chain']['name']}")
                
            heat_capacity = calculate_ionic_liquid_heat_capacity(il)
            if heat_capacity is not None:
                print(f"Calculated heat capacity: {heat_capacity:.2f} J/mol·K")
    
    if fragments_data:
        # Test fragment screening
        print("\nTesting fragment screening:")
        target_range = (100, 400)  # J/mol·K, can be parameterized
        screened = screen_fragments_by_heat_capacity(fragments_data, target_range)
        
        for frag_type, fragments in screened.items():
            if fragments:  # Only print if there are fragments that passed screening
                print(f"\n{frag_type.capitalize()} fragments within heat capacity range {target_range}:")
                for frag in fragments:
                    print(f"  {frag['name']}: {frag['estimated_heat_capacity']:.2f} J/mol·K")

if __name__ == "__main__":
    print("\n=== Heat Capacity Calculation Module Test ===")
    test_heat_capacity_calculations()
