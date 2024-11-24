import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from models.shortList_frag import fragments
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
        heavy_atoms = props['heavy_atoms']
        rotatable_bonds = props['rotatable_bonds']
        h_donors = props['h_bond_donors']
        h_acceptors = props['h_bond_acceptors']
        
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

def screen_fragments_by_heat_capacity(fragments: List[Dict], target_range: tuple) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated heat capacity
    Returns dictionary of fragments organized by type that meet the target range
    """
    min_capacity, max_capacity = target_range
    
    # Initialize dictionary for organized fragments
    fragments_data = {
        'Cation': [],
        'Anion': [],
        'Alkyl Chain': []
    }
    
    print("\nScreening fragments by heat capacity...")
    
    for fragment in fragments:
        frag_type = fragment['fragment_type']
        if frag_type in fragments_data:
            estimated_capacity = estimate_fragment_heat_capacity(fragment)
            
            if estimated_capacity is not None:
                # Check if estimated capacity contributes to target range
                if min_capacity/3 <= estimated_capacity <= max_capacity/2:  # Rough estimation since this is just one component
                    fragments_data[frag_type].append({
                        'name': fragment['name'],
                        'smiles': fragment['smiles'],
                        'heat_capacity': estimated_capacity,
                        'fragment_type': frag_type
                    })
    
    return fragments_data

def test_heat_capacity_calculations():
    """Test function to show detailed heat capacity calculations"""
    print("\nTesting with fragments from shortList_frag.py...")
    
    test_il = {
        'cation': {
            'name': 'Ethylimidazolium',
            'fragment_type': 'cation',
            'molecular_weight': 111.2,
            'h_bond_donors': 1,
            'h_bond_acceptors': 2
        },
        'anion': {
            'name': 'Tetrafluoroborate',
            'fragment_type': 'anion',
            'molecular_weight': 86.8,
            'h_bond_donors': 0,
            'h_bond_acceptors': 4
        },
        'alkyl_chain': {
            'name': 'Ethyl',
            'fragment_type': 'alkyl_chain',
            'molecular_weight': 29.1,
            'h_bond_donors': 0,
            'h_bond_acceptors': 0
        }
    }
    
    print(f"\nTest Ionic Liquid: {test_il['cation']['name']} {test_il['anion']['name']} with {test_il['alkyl_chain']['name']}")
    
    # Calculate full ionic liquid heat capacity
    il_cp = calculate_ionic_liquid_heat_capacity(test_il)
    
    if il_cp is not None:
        print(f"\nFinal Ionic Liquid Heat Capacity: {il_cp:.1f} J/mol·K")
    else:
        print("\nError: Could not calculate heat capacity")

if __name__ == "__main__":
    print("\n=== Heat Capacity Calculation Module Test ===")
    test_heat_capacity_calculations()
