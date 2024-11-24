import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from models.shortList_frag import fragments
from utils.utils import get_fragment_properties

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
        heavy_atoms = props['heavy_atoms']
        mw = props['molecular_weight']
        rotatable_bonds = props['rotatable_bonds']
        h_donors = props['h_bond_donors']
        h_acceptors = props['h_bond_acceptors']
        
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
        
        if fragment['fragment_type'] == 'Cation':
            # Imidazolium and similar cations
            base_volume = heavy_atoms * base_per_heavy_atom
            charge_factor = 0.95  # Cations are more compact
            flexibility_factor = 1.0 + (0.005 * rotatable_bonds)  # Reduced impact
            
        elif fragment['fragment_type'] == 'Anion':
            # Common IL anions
            base_volume = heavy_atoms * base_per_heavy_atom
            charge_factor = 1.05  # Anions slightly more diffuse
            flexibility_factor = 1.0 + (0.003 * rotatable_bonds)  # Reduced impact
            
            # Special case for fluorinated anions
            if 'fluor' in fragment['name'].lower():
                base_volume *= 1.02  # Slightly larger
                
        else:  # Alkyl Chain
            # Regular organic groups
            base_volume = heavy_atoms * base_per_heavy_atom
            charge_factor = 1.0
            flexibility_factor = 1.0 + (0.008 * rotatable_bonds)  # Reduced impact
            
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
        cation_h_donors = int(ionic_liquid['cation'].get('h_bond_donors', 0))
        cation_h_acceptors = int(ionic_liquid['cation'].get('h_bond_acceptors', 0))
        anion_h_donors = int(ionic_liquid['anion'].get('h_bond_donors', 0))
        anion_h_acceptors = int(ionic_liquid['anion'].get('h_bond_acceptors', 0))
        
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

def screen_fragments_by_density(fragments: List[Dict], target_range: tuple) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated density
    Returns dictionary of fragments organized by type that meet the target range
    """
    min_density, max_density = target_range
    
    # Initialize dictionary for organized fragments
    fragments_data = {
        'cation': [],
        'anion': [],
        'alkyl_chain': []
    }
    
    print("\nScreening fragments by density...")
    
    for fragment in fragments:
        frag_type = fragment['fragment_type']
        if frag_type in fragments_data:
            estimated_density = estimate_fragment_density(fragment)
            
            if estimated_density is not None:
                # Check if estimated density contributes to target range
                # Use wider range for screening since this is just one component
                if min_density * 0.7 <= estimated_density <= max_density * 1.3:
                    fragments_data[frag_type].append({
                        'name': fragment['name'],
                        'smiles': fragment['smiles'],
                        'density': estimated_density,
                        'fragment_type': frag_type
                    })
    
    return fragments_data

def test_density_calculations():
    """Test function to show detailed density calculations"""
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
        'alkyl': {
            'name': 'Ethyl',
            'fragment_type': 'alkyl_chain',
            'molecular_weight': 29.1,
            'h_bond_donors': 0,
            'h_bond_acceptors': 0
        }
    }
    
    print("\nTest Ionic Liquid:")
    print(f"Cation: {test_il['cation']['name']}")
    print(f"Anion: {test_il['anion']['name']}")
    print(f"Alkyl: {test_il['alkyl']['name']}\n")
    
    # Calculate individual contributions
    print("Cation contribution:")
    print(f"  - Molecular Weight: {test_il['cation']['molecular_weight']} g/mol")
    cation_volume = calculate_molecular_volume(test_il['cation'])
    if cation_volume:
        print(f"  - Molecular Volume: {cation_volume:.1f} cm³/mol\n")
    else:
        print("  - Error: Could not calculate molecular volume\n")
    
    print("Anion contribution:")
    print(f"  - Molecular Weight: {test_il['anion']['molecular_weight']} g/mol")
    anion_volume = calculate_molecular_volume(test_il['anion'])
    if anion_volume:
        print(f"  - Molecular Volume: {anion_volume:.1f} cm³/mol\n")
    else:
        print("  - Error: Could not calculate molecular volume\n")
    
    print("Alkyl contribution:")
    print(f"  - Molecular Weight: {test_il['alkyl']['molecular_weight']} g/mol")
    alkyl_volume = calculate_molecular_volume(test_il['alkyl'])
    if alkyl_volume:
        print(f"  - Molecular Volume: {alkyl_volume:.1f} cm³/mol\n")
    else:
        print("  - Error: Could not calculate molecular volume\n")
    
    print("Calculating density for Ethylimidazolium Tetrafluoroborate with Ethyl:")
    density = calculate_ionic_liquid_density(test_il)
    
    print("\nComponent volumes:")
    print(f"  - Cation: {f'{cation_volume:.1f}' if cation_volume is not None else 'None'} cm³/mol")
    print(f"  - Anion: {f'{anion_volume:.1f}' if anion_volume is not None else 'None'} cm³/mol")
    print(f"  - Alkyl: {f'{alkyl_volume:.1f}' if alkyl_volume is not None else 'None'} cm³/mol")
    
    if density is not None:
        print(f"\nFinal Ionic Liquid Density: {density:.1f} kg/m³")
    else:
        print("\nError: Could not calculate final density")

if __name__ == "__main__":
    print("\n=== Density Calculation Module Test ===")
    test_density_calculations()
