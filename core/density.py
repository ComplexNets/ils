import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional, Tuple
from models.shortList_frag import fragments
from utils.utils import get_fragment_properties
from core.combine_fragments import get_filtered_fragments, combine_fragments

# Group contribution parameters for density calculation
group_params = {
    'cation': {
        'imidazolium': {'volume': 25.0, 'weight_factor': 1.2},
        'pyridinium': {'volume': 24.0, 'weight_factor': 1.1}
    },
    'anion': {
        'tetrafluoroborate': {'volume': 23.0, 'weight_factor': 1.3},
        'hexafluorophosphate': {'volume': 24.0, 'weight_factor': 1.4},
        'bis(trifluoromethylsulfonyl)imide': {'volume': 28.0, 'weight_factor': 1.5}
    },
    'alkyl_chain': {
        'base': {'volume': 21.0, 'weight_factor': 1.0},
        'ch2_increment': {'volume': 16.0, 'weight_factor': 0.9}
    }
}

def calculate_molecular_volume(fragment: Dict) -> Optional[float]:
    """Calculate molecular volume for a single fragment"""
    try:
        # Get properties from database
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props:
            print(f"  No properties found for {fragment['name']}")
            return None
            
        # Get required properties
        heavy_atoms = props['heavy_atom_count']
        mw = props['molecular_weight']
        rotatable_bonds = props['rotatable_bond_count']
        h_donors = props['hydrogen_bond_donor_count']
        h_acceptors = props['hydrogen_bond_acceptor_count']
        
        if not all([heavy_atoms, mw]):
            print(f"  Missing required properties for {fragment['name']}")
            return None
            
        # Base volume calculation
        # Each atom type has different contribution
        # Values calibrated based on experimental data for [EMIM][BF4]
        # [EMIM][BF4] density ≈ 1240 kg/m³ at 298.15K
        # MW = 197.97 g/mol
        # Therefore volume ≈ 160 cm³/mol
        # With 13 heavy atoms total, that's about 12.3 cm³/mol per heavy atom
        base_per_heavy_atom = 12.3  # cm³/mol per heavy atom
        
        fragment_type = fragment['fragment_type'].lower()
        
        if fragment_type == 'cation':
            # Imidazolium and similar cations
            base_volume = heavy_atoms * base_per_heavy_atom
            charge_factor = 0.95  # Cations are more compact
            flexibility_factor = 1.0 + (0.005 * rotatable_bonds)  # Reduced impact
            
        elif fragment_type == 'anion':
            # Common IL anions
            base_volume = heavy_atoms * base_per_heavy_atom
            charge_factor = 1.05  # Anions slightly more diffuse
            flexibility_factor = 1.0 + (0.003 * rotatable_bonds)  # Reduced impact
            
            # Special case for fluorinated anions
            if 'fluor' in fragment['name'].lower():
                base_volume *= 1.02  # Slightly larger
                
        elif fragment_type == 'alkyl_chain':
            # Regular organic groups
            base_volume = heavy_atoms * base_per_heavy_atom
            charge_factor = 1.0
            flexibility_factor = 1.0 + (0.008 * rotatable_bonds)  # Reduced impact
            
        else:
            print(f"  Unknown fragment type: {fragment['fragment_type']}")
            return None
            
        # Molecular weight correction
        # Heavier atoms generally take more space
        mw_per_atom = mw / heavy_atoms
        mw_factor = 0.95 + (mw_per_atom / 200.0)  # Adjusted from 100.0
        
        # H-bonding affects packing
        h_bond_factor = 1.0 - (0.005 * (h_donors + h_acceptors))  # Reduced impact
        
        # Calculate final volume
        volume = (base_volume * 
                 charge_factor * 
                 flexibility_factor * 
                 mw_factor * 
                 h_bond_factor)
                 
        print(f"  Volume calculation for {fragment['name']}:")
        print(f"    Base volume: {base_volume:.1f}")
        print(f"    Charge factor: {charge_factor:.2f}")
        print(f"    Flexibility factor: {flexibility_factor:.2f}")
        print(f"    MW factor: {mw_factor:.2f}")
        print(f"    H-bond factor: {h_bond_factor:.2f}")
        print(f"    Final volume: {volume:.1f}")
        
        return volume
        
    except Exception as e:
        print(f"Error calculating molecular volume: {e}")
        return None

def estimate_fragment_density(fragment: Dict) -> Optional[float]:
    """
    Estimate density for a single fragment using molecular weight and volume
    Returns density in kg/m³
    """
    try:
        # First try to get density directly from database
        # density = get_estimated_density(fragment['name'], fragment['fragment_type'])
        # if density is not None:
        #     return density
            
        # If no density in database, calculate from molecular weight and volume
        molecular_volume = calculate_molecular_volume(fragment)
        if molecular_volume is None or molecular_volume <= 0:
            return None
            
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props or not props['molecular_weight']:
            return None
            
        # Convert molecular weight from g/mol to kg/mol
        mw_kg = props['molecular_weight'] / 1000
        
        # Convert volume to m³/mol
        volume_m3 = molecular_volume / 1000000
        
        # Calculate density in kg/m³
        density = mw_kg / volume_m3
        return density
        
    except Exception as e:
        print(f"Error estimating fragment density: {e}")
        return None

def calculate_ionic_liquid_density(ionic_liquid: Dict) -> Optional[float]:
    """
    Calculate density for a complete ionic liquid using an improved group contribution method
    Args:
        ionic_liquid: Dictionary containing cation, anion, and alkyl chain information
    Returns:
        Density in kg/m³
    """
    try:
        print(f"\nCalculating density for {ionic_liquid['cation'].get('name', 'Unknown')} {ionic_liquid['anion'].get('name', 'Unknown')} with {ionic_liquid['alkyl_chain'].get('name', 'Unknown')}:")
        
        # Calculate volume for each component
        cation_volume = calculate_molecular_volume(ionic_liquid['cation'])
        anion_volume = calculate_molecular_volume(ionic_liquid['anion'])
        alkyl_volume = calculate_molecular_volume(ionic_liquid['alkyl_chain'])
        
        print(f"\nComponent volumes:")
        print(f"  Cation: {f'{cation_volume:.1f}' if cation_volume is not None else 'None'} cm³/mol")
        print(f"  Anion: {f'{anion_volume:.1f}' if anion_volume is not None else 'None'} cm³/mol")
        print(f"  Alkyl: {f'{alkyl_volume:.1f}' if alkyl_volume is not None else 'None'} cm³/mol")
        
        if any(v is None for v in [cation_volume, anion_volume, alkyl_volume]):
            print("Error: One or more volumes are None")
            return None
            
        # Get molecular weights with better error handling
        try:
            cation_mw = float(ionic_liquid['cation'].get('molecular_weight', 0))
            anion_mw = float(ionic_liquid['anion'].get('molecular_weight', 0))
            alkyl_mw = float(ionic_liquid['alkyl_chain'].get('molecular_weight', 0))
        except (ValueError, TypeError) as e:
            print(f"Error converting molecular weights: {e}")
            return None
        
        # Print molecular weights for debugging
        print(f"\nMolecular weights:")
        print(f"  Cation: {cation_mw:.1f}")
        print(f"  Anion: {anion_mw:.1f}")
        print(f"  Alkyl: {alkyl_mw:.1f}")
        
        if cation_mw <= 0 or anion_mw <= 0 or alkyl_mw <= 0:
            print("Error: Invalid molecular weights (zero or negative)")
            return None
        
        # Calculate total mass
        total_mass = cation_mw + anion_mw + alkyl_mw
        print(f"Total mass: {total_mass:.1f}")
        
        # Base volume is sum of component volumes
        base_volume = cation_volume + anion_volume + alkyl_volume
        
        # Calculate interaction volume reduction
        # Ions pack more efficiently due to electrostatic interactions
        cation_charge = abs(float(ionic_liquid['cation'].get('charge', 1)))
        anion_charge = abs(float(ionic_liquid['anion'].get('charge', -1)))
        charge_interaction = cation_charge * anion_charge
        
        # Stronger charge interactions = tighter packing
        # Calibrated based on [EMIM][BF4]
        packing_factor = 0.98 - (0.01 * charge_interaction)
        
        # H-bonding also affects packing
        cation_h_donors = int(ionic_liquid['cation'].get('hydrogen_bond_donor_count', 0))
        cation_h_acceptors = int(ionic_liquid['cation'].get('hydrogen_bond_acceptor_count', 0))
        anion_h_donors = int(ionic_liquid['anion'].get('hydrogen_bond_donor_count', 0))
        anion_h_acceptors = int(ionic_liquid['anion'].get('hydrogen_bond_acceptor_count', 0))
        
        h_bond_pairs = min(cation_h_donors, anion_h_acceptors) + min(anion_h_donors, cation_h_acceptors)
        h_bond_factor = 1.0 - (0.005 * h_bond_pairs)
        
        # Final volume includes packing effects
        final_volume = base_volume * packing_factor * h_bond_factor
        
        print(f"\nVolume calculations:")
        print(f"  Base volume: {base_volume:.1f}")
        print(f"  Packing factor: {packing_factor:.3f}")
        print(f"  H-bond factor: {h_bond_factor:.3f}")
        print(f"  Final volume: {final_volume:.1f}")
        
        # Convert to density (kg/m³)
        # mass in g/mol, volume in cm³/mol
        # 1 g/cm³ = 1000 kg/m³
        density = (total_mass / final_volume) * 1000
        
        print(f"\nFinal density: {density:.1f} kg/m³")
        return density
        
    except Exception as e:
        print(f"Error calculating ionic liquid density: {e}")
        return None

def screen_fragments_by_density(fragments_data: Dict[str, List[Dict]], target_range: Tuple[float, float]) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated density
    Args:
        fragments_data: Dictionary of fragments organized by type (cation, anion, alkyl_chain)
        target_range: Tuple of (min_density, max_density) in kg/m³
    Returns:
        Dictionary of fragments organized by type that meet the target range
    """
    min_density, max_density = target_range
    
    # Initialize output dictionary with same structure as input
    screened_fragments = {frag_type: [] for frag_type in fragments_data.keys()}
    
    for frag_type, fragments in fragments_data.items():
        for fragment in fragments:
            estimated_density = estimate_fragment_density(fragment)
            
            if estimated_density is not None:
                # Use configurable margin for screening
                margin = 0.3  # Can be adjusted or passed as parameter
                if min_density * (1 - margin) <= estimated_density <= max_density * (1 + margin):
                    screened_fragments[frag_type].append({
                        **fragment,  # Preserve all original properties
                        'estimated_density': estimated_density
                    })
    
    return screened_fragments

def test_density_calculations(fragments_data: Optional[Dict[str, List[Dict]]] = None, 
                            combinations: Optional[List[Dict]] = None,
                            num_test_combinations: int = 3):
    """
    Test function to show detailed density calculations
    Args:
        fragments_data: Optional dictionary of fragments to test screening
        combinations: Optional list of ionic liquid combinations to test
        num_test_combinations: Number of combinations to test
    """
    print("\n=== Testing Density Calculations ===")
    
    if combinations:
        # Test provided combinations
        for i, il in enumerate(combinations[:num_test_combinations]):
            print(f"\nTesting combination {i+1}:")
            print(f"Cation: {il['cation']['name']}")
            print(f"Anion: {il['anion']['name']}")
            if 'alkyl_chain' in il:
                print(f"Alkyl chain: {il['alkyl_chain']['name']}")
                
            density = calculate_ionic_liquid_density(il)
            if density is not None:
                print(f"Calculated density: {density:.2f} kg/m³")
    
    if fragments_data:
        # Test fragment screening
        print("\nTesting fragment screening:")
        target_range = (800, 2000)  # kg/m³, can be parameterized
        screened = screen_fragments_by_density(fragments_data, target_range)
        
        for frag_type, fragments in screened.items():
            if fragments:  # Only print if there are fragments that passed screening
                print(f"\n{frag_type.capitalize()} fragments within density range {target_range}:")
                for frag in fragments:
                    print(f"  {frag['name']}: {frag['estimated_density']:.2f} kg/m³")

if __name__ == "__main__":
    print("\n=== Density Calculation Module Test ===")
    test_density_calculations()
