"""
Density calculation module for ionic liquids using Group Contribution Methods.
Primary method: Modified Ye & Shreeve method with temperature correction (Gardas & Coutinho).
"""

from typing import Dict, Optional, List, Tuple
import math
from utils.utils import get_fragment_properties

# Constants for temperature correction (from Gardas & Coutinho)
ALPHA = 6.652e-4  # K^-1 (thermal expansion coefficient)
BETA = -5.919e-4  # MPa^-1 (isothermal compressibility)
T_REF = 298.15    # K (reference temperature)
P_REF = 0.1       # MPa (reference pressure)

# Atomic volumes (Å³) based on more recent literature
ATOMIC_VOLUMES = {
    'C': 16.35,   # Carbon (sp3)
    'C=': 14.7,   # Carbon (sp2)
    'C#': 13.9,   # Carbon (sp)
    'H': 6.71,    # Hydrogen
    'N': 14.39,   # Nitrogen (sp3)
    'N=': 13.2,   # Nitrogen (sp2)
    'O': 12.43,   # Oxygen
    'F': 13.31,   # Fluorine
    'P': 24.87,   # Phosphorus
    'S': 24.43,   # Sulfur
    'B': 18.32,   # Boron
    'Cl': 22.45,  # Chlorine
    'Br': 26.52,  # Bromine
}

# Ionic radii correction for halides (Å)
IONIC_RADII = {
    'Cl': 1.81,   # Chloride ion
    'Br': 1.96,   # Bromide ion
}

# Ring corrections for cyclic structures
RING_CORRECTIONS = {
    'aromatic_5': 0.816,  # 5-membered aromatic ring
    'aromatic_6': 0.896,  # 6-membered aromatic ring
    'aliphatic_5': 0.847, # 5-membered aliphatic ring
    'aliphatic_6': 0.915  # 6-membered aliphatic ring
}

# Correction factors for different chemical environments
ENVIRONMENT_FACTORS = {
    'aromatic': 0.97,     # Aromatic rings are more compact
    'conjugated': 0.98,   # Conjugated systems
    'charged': 0.95,      # Ionic species are more compact
    'h_bonding': 0.96,    # H-bonding reduces volume
    'halide': 1.15,       # Halide anions have larger effective volumes
}

def estimate_molecular_volume_from_formula(fragment: Dict) -> float:
    """
    Estimate molecular volume from molecular formula and structural features.
    Args:
        fragment: Dictionary containing fragment properties
    Returns:
        Estimated molecular volume in Å³
    """
    try:
        # Get properties from database
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props:
            print(f"Warning: No properties found for fragment {fragment['name']}")
            return 0.0
            
        # Get structural information
        smiles = props.get('smiles', '') or fragment.get('smiles', '')  # Try both locations
        name = fragment.get('name', '').lower()
        frag_type = fragment.get('fragment_type', '')
        
        if not smiles:
            print(f"Warning: No SMILES found for fragment {name}")
            return 0.0
            
        h_donors = int(props.get('hydrogen_bond_donor_count', 0))
        h_acceptors = int(props.get('hydrogen_bond_acceptor_count', 0))
        charge = float(props.get('charge', 0))
        heavy_atoms = int(props.get('heavy_atom_count', 0))
        
        # Special handling for alkyl chains
        if frag_type == 'alkyl_chain':
            # Count carbons in chain
            carbon_count = smiles.count('C')
            # Calculate volume based on chain length
            volume = carbon_count * (ATOMIC_VOLUMES['C'] + 2 * ATOMIC_VOLUMES['H'])
            # Add terminal hydrogens
            volume += ATOMIC_VOLUMES['H']
            return volume
            
        # Improved aromaticity and ring detection
        is_aromatic = ('c1' in smiles.lower() or 'n1' in smiles.lower() or 
                      'pyridinium' in name or 'imidazolium' in name)
        
        # Ring size detection
        ring_size = 0
        if 'c1ccc' in smiles.lower() or 'n1ccc' in smiles.lower():
            ring_size = 5
        elif 'c1cccc' in smiles.lower() or 'n1cccc' in smiles.lower() or 'pyridinium' in name:
            ring_size = 6
        elif 'imidazolium' in name:
            ring_size = 5
        
        # Base volume calculation
        volume = 0.0
        
        # Special handling for ionic fragments
        if frag_type == 'anion':
            # Handle halide anions
            if '[Cl-]' in smiles:
                radius = IONIC_RADII['Cl']
                volume = (4/3) * math.pi * (radius ** 3)
                volume *= ENVIRONMENT_FACTORS['halide']
                return volume
            elif '[Br-]' in smiles:
                radius = IONIC_RADII['Br']
                volume = (4/3) * math.pi * (radius ** 3)
                volume *= ENVIRONMENT_FACTORS['halide']
                return volume
            elif 'B' in smiles:  # Handle tetrafluoroborate
                volume = ATOMIC_VOLUMES['B'] + 4 * ATOMIC_VOLUMES['F']
                volume *= ENVIRONMENT_FACTORS['charged']
                return volume
        
        # Calculate volume based on SMILES pattern matching
        atoms = {
            'C': smiles.count('C') - smiles.count('Cl'),  # Exclude chlorine from carbon count
            'N': smiles.count('N'),
            'O': smiles.count('O'),
            'F': smiles.count('F'),
            'P': smiles.count('P'),
            'S': smiles.count('S'),
            'B': smiles.count('B'),
            'Cl': smiles.count('Cl'),
            'Br': smiles.count('Br'),
        }
        
        # Special handling for pyridinium and imidazolium rings
        if 'pyridinium' in name:
            # Base pyridinium ring volume (5 carbons + 1 nitrogen)
            volume = 5 * ATOMIC_VOLUMES['C='] + ATOMIC_VOLUMES['N=']
            # Add substituent volumes
            if 'methyl' in name:
                volume += ATOMIC_VOLUMES['C'] + 3 * ATOMIC_VOLUMES['H']
            elif 'ethyl' in name:
                volume += 2 * ATOMIC_VOLUMES['C'] + 5 * ATOMIC_VOLUMES['H']
            elif 'propyl' in name:
                volume += 3 * ATOMIC_VOLUMES['C'] + 7 * ATOMIC_VOLUMES['H']
            elif 'butyl' in name:
                volume += 4 * ATOMIC_VOLUMES['C'] + 9 * ATOMIC_VOLUMES['H']
        else:
            # Standard atom-based calculation for other fragments
            for atom, count in atoms.items():
                volume += ATOMIC_VOLUMES[atom] * count
            
            # Count hydrogens (approximate)
            explicit_h = smiles.count('H')
            implicit_h = heavy_atoms * 2  # Rough estimate
            volume += (explicit_h + implicit_h) * ATOMIC_VOLUMES['H']
        
        # Apply correction factors
        if charge != 0:
            volume *= ENVIRONMENT_FACTORS['charged']
        if h_donors + h_acceptors > 0:
            volume *= ENVIRONMENT_FACTORS['h_bonding']
        if is_aromatic:
            volume *= ENVIRONMENT_FACTORS['aromatic']
        
        return volume
        
    except Exception as e:
        print(f"Error calculating molecular volume: {str(e)}")
        return 0.0

def calculate_molecular_volume(fragment: Dict) -> float:
    """
    Calculate molecular volume for a fragment using its volume contributions.
    Args:
        fragment: Dictionary containing fragment properties
    Returns:
        Molecular volume in Å³
    """
    try:
        # First try to get pre-calculated volume from database
        if 'molecular_volume' in fragment:
            return float(fragment['molecular_volume'])
        
        # If not available, estimate from molecular formula and structure
        return estimate_molecular_volume_from_formula(fragment)
        
    except Exception as e:
        print(f"Error calculating molecular volume: {str(e)}")
        return 0.0

def calculate_density(cation: Dict, anion: Dict, alkyl_chain: Dict, functional_group: Dict = None,
                     temperature: float = T_REF, pressure: float = P_REF) -> Optional[float]:
    """
    Calculate ionic liquid density using modified Ye & Shreeve method.
    All calculations are done at standard conditions (298.15 K, 0.1 MPa).
    Temperature and pressure parameters are kept for API compatibility but not recommended for use.
    
    Args:
        cation: Dictionary containing cation properties
        anion: Dictionary containing anion properties
        alkyl_chain: Dictionary containing alkyl chain properties
        functional_group: Dictionary containing functional group properties (optional)
        temperature: Temperature in K (fixed at 298.15 K)
        pressure: Pressure in MPa (fixed at 0.1 MPa)
    
    Returns:
        Density in kg/m³ or None if calculation fails
    """
    try:
        # Use standard conditions
        temperature = T_REF  # 298.15 K
        pressure = P_REF     # 0.1 MPa
        
        # Calculate molecular weights with validation
        cation_props = get_fragment_properties(cation['name'], cation['fragment_type'])
        anion_props = get_fragment_properties(anion['name'], anion['fragment_type'])
        alkyl_props = get_fragment_properties(alkyl_chain['name'], alkyl_chain['fragment_type'])
        
        # Initialize functional group properties
        functional_group_props = None
        mw_functional_group = 0.0
        v_functional_group = 0.0
        
        # Process functional group if present
        if functional_group:
            functional_group_props = get_fragment_properties(functional_group['name'], functional_group['fragment_type'])
            if functional_group_props:
                mw_functional_group = float(functional_group_props['molecular_weight'])
                # If molecular weight is not available, estimate based on the group type
                if mw_functional_group <= 0:
                    group_name = functional_group['name'].lower()
                    if 'hydroxyl' in group_name:
                        mw_functional_group = 17.0  # -OH group
                    elif 'amino' in group_name:
                        mw_functional_group = 16.0  # -NH2 group
                    else:
                        mw_functional_group = 15.0  # Default value
        
        if not all([cation_props, anion_props, alkyl_props]):
            print("Warning: Could not get properties for all fragments")
            return None
            
        mw_cation = float(cation_props['molecular_weight'])
        mw_anion = float(anion_props['molecular_weight'])
        mw_alkyl = float(alkyl_props['molecular_weight'])
        
        if any(mw <= 0 for mw in [mw_cation, mw_anion, mw_alkyl]):
            print("Warning: Invalid molecular weight(s)")
            return None
            
        total_mw = mw_cation + mw_anion + mw_alkyl + mw_functional_group
        
        # Calculate molecular volumes with validation
        v_cation = calculate_molecular_volume(cation)
        v_anion = calculate_molecular_volume(anion)
        v_alkyl = calculate_molecular_volume(alkyl_chain)
        
        # Calculate functional group volume if present
        if functional_group:
            v_functional_group = calculate_molecular_volume(functional_group)
            # If volume calculation fails, estimate based on the group type
            if v_functional_group <= 0:
                group_name = functional_group['name'].lower()
                if 'hydroxyl' in group_name:
                    v_functional_group = 10.0  # Approximate volume for -OH group
                elif 'amino' in group_name:
                    v_functional_group = 15.0  # Approximate volume for -NH2 group
                else:
                    v_functional_group = 12.0  # Default value
        
        if any(v <= 0 for v in [v_cation, v_anion, v_alkyl]):
            print("Warning: Invalid molecular volume(s)")
            return None
            
        total_volume = v_cation + v_anion + v_alkyl + v_functional_group
        
        # Convert volume from Å³ to nm³
        total_volume_nm3 = total_volume / 1000
        
        # Calculate density with proper unit conversions
        # Convert from g/mol to kg/mol, nm³ to m³
        density = (total_mw / 1000) / (total_volume_nm3 * 1e-27 * 6.02214076e23)
        
        return density
        
    except Exception as e:
        print(f"Error calculating density: {str(e)}")
        return None

def validate_density(density: Optional[float], min_density: float = 700, max_density: float = 2200) -> bool:
    """
    Validate if the calculated density is within reasonable bounds for ionic liquids.
    Args:
        density: Calculated density in kg/m³
        min_density: Minimum acceptable density (default: 700 kg/m³)
        max_density: Maximum acceptable density (default: 2200 kg/m³)
    Returns:
        True if density is valid, False otherwise
    """
    if density is None:
        return False
    return min_density <= density <= max_density

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
    screened_fragments = {
        'cation': [],
        'anion': [],
        'alkyl_chain': []
    }
    
    # For each fragment type
    for frag_type, fragments in fragments_data.items():
        for fragment in fragments:
            # Create a test combination using this fragment
            test_combo = {
                'cation': fragments_data['cation'][0],
                'anion': fragments_data['anion'][0],
                'alkyl_chain': fragments_data['alkyl_chain'][0]
            }
            test_combo[frag_type] = fragment
            
            # Calculate density
            density = calculate_density(**test_combo)
            
            # Apply a margin (e.g., 40%) to the target range for coarse screening
            screening_min = min_density * 0.60 # Allow 40% below min target
            screening_max = max_density * 1.40 # Allow 40% above max target

            # Check if the density of the test combination is within the broadened range
            if density is not None and screening_min <= density <= screening_max:
                screened_fragments[frag_type].append(fragment)
    
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
    print("\nTesting density calculations...")
    
    if combinations:
        for combo in combinations[:num_test_combinations]:
            print(f"\nTesting combination: {combo['name']}")
            density = calculate_density(combo['cation'], combo['anion'], combo['alkyl_chain'])
            if density is not None:
                print(f"Calculated density: {density:.2f} kg/m³")
                print(f"Valid: {validate_density(density)}")
    
    if fragments_data:
        print("\nTesting fragment screening...")
        target_range = (800, 2000)  # Example target range
        screened = screen_fragments_by_density(fragments_data, target_range)
        for frag_type, fragments in screened.items():
            print(f"\n{frag_type}: {len(fragments)} fragments within density range")
            for frag in fragments[:3]:  # Show first 3 fragments
                print(f"  - {frag['name']}")

if __name__ == "__main__":
    print("\n=== Density Calculation Module Test ===")
    test_density_calculations()
