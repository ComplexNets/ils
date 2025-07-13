import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from .combine_fragments import get_filtered_fragments, combine_fragments
from utils.utils import get_fragment_properties

# UNIFAC group parameters for heat capacity calculation
# Values derived from literature data for common ionic liquids
unifac_params = {
    'IM+': {'A': 130.0, 'B': 0.15, 'C': -0.0002},  # Imidazolium core contribution
    'PY+': {'A': 125.0, 'B': 0.14, 'C': -0.0002},  # Pyridinium (scaled relative to IM+)
    'BF4-': {'A': 110.0, 'B': 0.12, 'C': -0.0001}, # From [EMIM][BF4] and [BMIM][BF4] data
    'PF6-': {'A': 120.0, 'B': 0.13, 'C': -0.0001}, # Literature data for PF6-based ILs
    'NTf2-': {'A': 140.0, 'B': 0.12, 'C': -0.0001}, # For bis(trifluoromethylsulfonyl)imide
    'SCN-': {'A': 90.0, 'B': 0.10, 'C': -0.0001},  # Literature data for thiocyanate ILs
    'DCA-': {'A': 95.0, 'B': 0.11, 'C': -0.0001}, # Literature data for dicyanamide ILs
    'CH2': {'A': 30.0, 'B': 0.05, 'C': -0.0001},   # Derived from alkyl chain length variations
    'CH3': {'A': 35.0, 'B': 0.05, 'C': -0.0001},   # Terminal methyl contribution
    # Functional group parameters
    'OH': {'A': 25.0, 'B': 0.04, 'C': -0.00005},   # Hydroxyl group
    'NH2': {'A': 30.0, 'B': 0.045, 'C': -0.00006}, # Amino group
    'COOH': {'A': 40.0, 'B': 0.05, 'C': -0.00007}, # Carboxylic acid group
    'CN': {'A': 35.0, 'B': 0.04, 'C': -0.00005},   # Nitrile group
    'SH': {'A': 28.0, 'B': 0.04, 'C': -0.00005},   # Thiol group
    'F': {'A': 18.0, 'B': 0.02, 'C': -0.00003},    # Fluoro group
    'Cl': {'A': 22.0, 'B': 0.025, 'C': -0.00004},  # Chloro group
    'Br': {'A': 26.0, 'B': 0.03, 'C': -0.00004},   # Bromo group
    'I': {'A': 30.0, 'B': 0.035, 'C': -0.00005},   # Iodo group
    'Ether': {'A': 20.0, 'B': 0.03, 'C': -0.00004} # Ether group
}

def estimate_fragment_heat_capacity(fragment: Dict) -> Optional[float]:
    """
    Estimate heat capacity for a single fragment
    Based on literature values for common ionic liquids:
    - [EMIM][BF4]: ~320 J/K·mol at 298.15K
    - [BMIM][PF6]: ~350 J/K·mol at 298.15K
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
        if fragment['fragment_type'] == 'cation':
            base_adjustment = 1.15
        elif fragment['fragment_type'] == 'anion':
            # Apply additional scaling for large anions like NTf2
            if fragment['name'] in ['bis(trifluoromethanesulfonyl)imide', 'NTf2']:
                base_adjustment = 0.95  # Reduced adjustment for NTf2
            else:
                base_adjustment = 1.05
        elif fragment['fragment_type'] == 'functional_group':
            base_adjustment = 1.0  # Default adjustment for functional groups
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
        print(f"    Total: {total_cp:.1f} J/K·mol") # Added unit here for clarity

        return total_cp

    except Exception as e:
        print(f"Error estimating fragment heat capacity: {e}")
        return None

def get_fragment_type(fragment_name: str, fragment_type: str) -> str:
    """Map fragment names to their UNIFAC types"""
    name_lower = fragment_name.lower()
    
    if fragment_type == 'cation':
        if 'imidazolium' in name_lower:
            return 'IM+'
        elif 'pyridinium' in name_lower:
            return 'PY+'
        return 'IM+'  # Default to imidazolium if unknown
        
    elif fragment_type == 'anion':
        if 'tetrafluoroborate' in name_lower or 'bf4' in name_lower:
            return 'BF4-'
        elif 'hexafluorophosphate' in name_lower or 'pf6' in name_lower:
            return 'PF6-'
        elif 'thiocyanate' in name_lower or 'scn' in name_lower:
            return 'SCN-'
        elif 'dicyanamide' in name_lower or 'dca' in name_lower:
            return 'DCA-'
        elif 'bis(trifluoromethylsulfonyl)imide' in name_lower or 'ntf2' in name_lower or 'tf2n' in name_lower:
            return 'NTf2-'
        elif 'triflate' in name_lower or 'trifluoromethanesulfonate' in name_lower:
            return 'NTf2-'  # Using NTf2 parameters as approximation
        return 'BF4-'  # Default to BF4 if unknown
    
    elif fragment_type == 'alkyl_chain':
        if name_lower.startswith('methyl'):
            return 'CH3'
        return 'CH2'  # Default to CH2 for other alkyl chains
        
    elif fragment_type == 'functional_group':
        if 'hydroxyl' in name_lower or 'hydroxy' in name_lower or 'alcohol' in name_lower:
            return 'OH'
        elif 'amino' in name_lower or 'amine' in name_lower:
            return 'NH2'
        elif 'carboxyl' in name_lower or 'acid' in name_lower:
            return 'COOH'
        elif 'nitrile' in name_lower or 'cyano' in name_lower:
            return 'CN'
        elif 'thiol' in name_lower or 'mercapto' in name_lower:
            return 'SH'
        elif 'fluoro' in name_lower:
            return 'F'
        elif 'chloro' in name_lower:
            return 'Cl'
        elif 'bromo' in name_lower:
            return 'Br'
        elif 'iodo' in name_lower:
            return 'I'
        elif 'ether' in name_lower or 'methoxy' in name_lower or 'ethoxy' in name_lower:
            return 'Ether'
        # Default case - return the most common functional group
        return 'OH'
        
    return 'CH2'  # Default case

def unifac_contribution(group_type: str, T: float = 298.15) -> float:
    """
    Calculate UNIFAC contribution with temperature dependence.
    These values are now calibrated for individual fragment contributions
    rather than representing a significant portion of the overall ionic liquid.
    
    Args:
        group_type: The UNIFAC group type identifier
        T: Temperature in K (default: 298.15K)
    Returns:
        Heat capacity contribution in J/K·mol
    """
    if group_type not in unifac_params:
        return 0.0
    
    # Calculate temperature deviation from reference
    dT = T - 298.15
    
    # Calculate the base contribution from parameters
    base_contrib = unifac_params[group_type]['A'] + unifac_params[group_type]['B'] * dT + unifac_params[group_type]['C'] * dT**2
    
    # Apply fragment-specific scaling factor based on group type
    # This helps ensure the contributions are appropriate for individual fragments
    # rather than representing whole ionic liquids
    if group_type in ['IM+', 'PY+']:  # Cation cores
        scaling = 1.0  # Base value already adjusted
    elif group_type in ['BF4-', 'PF6-', 'NTf2-', 'SCN-', 'DCA-']:  # Anions
        scaling = 1.0  # Base value already adjusted
    elif group_type in ['CH2', 'CH3']:  # Alkyl groups
        scaling = 1.0  # Base value already adjusted
    else:  # Functional groups
        scaling = 1.0  # Base value already adjusted
    
    return base_contrib * scaling

def add_alkyl_contributions(chain_length: int, T: float = 298.15) -> float:
    """Calculate heat capacity contribution from alkyl chain"""
    ch2_count = max(0, chain_length - 1)  # Number of CH2 groups
    ch3_count = 1  # Always one CH3 terminal group
    
    ch2_contrib = ch2_count * unifac_contribution('CH2', T)
    ch3_contrib = ch3_count * unifac_contribution('CH3', T)
    
    return ch2_contrib + ch3_contrib

def calculate_fragment_forward_cp(fragment: Dict, T: float = 298.15) -> Optional[float]:
    """
    Calculate heat capacity for a fragment using a balanced approach that combines
    UNIFAC group contributions with molecular property contributions, while avoiding
    double-counting effects.
    
    Args:
        fragment: Dictionary containing fragment properties
        T: Temperature in K (default: 298.15K)
    Returns:
        Heat capacity contribution in J/K·mol
    """
    try:
        # Get properties from database
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props:
            print(f"  No properties found for {fragment['name']}")
            return None
            
        # Get fragment type and properties
        frag_type = fragment.get('fragment_type', '')
        name = fragment.get('name', '').lower()
        
        # Get molecular properties for additional contributions
        mw = props.get('molecular_weight', 0)
        heavy_atoms = props.get('heavy_atom_count', 0)
        rotatable_bonds = props.get('rotatable_bond_count', 0)
        h_donors = props.get('hydrogen_bond_donor_count', 0)
        h_acceptors = props.get('hydrogen_bond_acceptor_count', 0)
        
        # Calculate total heat capacity using a balanced approach
        total_cp = 0.0
        
        # 1. Base contribution from UNIFAC parameters (primary contribution)
        unifac_type = get_fragment_type(fragment['name'], frag_type)
        if unifac_type:
            unifac_cp = unifac_contribution(unifac_type, T)
            total_cp += unifac_cp
            print(f"  UNIFAC contribution for {fragment['name']}: {unifac_cp:.1f} J/K·mol")
        
        # 2. Add a reduced molecular weight contribution (secondary contribution)
        # Using a smaller coefficient than in the reverse method to avoid double-counting
        if mw > 0:
            mw_contribution = 0.3 * mw  # Reduced coefficient
            total_cp += mw_contribution
            print(f"  MW contribution for {fragment['name']}: {mw_contribution:.1f} J/K·mol")
        
        # 3. Add a small contribution for conformational flexibility
        # This captures effects not fully represented in the UNIFAC parameters
        if rotatable_bonds > 0:
            flexibility_contribution = 1.5 * rotatable_bonds  # Reduced coefficient
            total_cp += flexibility_contribution
            print(f"  Flexibility contribution for {fragment['name']}: {flexibility_contribution:.1f} J/K·mol")
        
        # 4. Add a small contribution for hydrogen bonding capability
        # This is important for functional groups and affects heat capacity
        if h_donors > 0 or h_acceptors > 0:
            h_bond_contribution = 1.0 * h_donors + 0.5 * h_acceptors  # Reduced coefficients
            total_cp += h_bond_contribution
            print(f"  H-bond contribution for {fragment['name']}: {h_bond_contribution:.1f} J/K·mol")
        
        # If we have a total contribution, return it
        if total_cp > 0:
            print(f"  Total heat capacity for {fragment['name']}: {total_cp:.1f} J/K·mol")
            return total_cp
        
        # Fallback if all else fails: use Joback-like method based on heavy atoms
        if heavy_atoms > 0:
            joback_cp = 33.0 * heavy_atoms  # Average contribution per heavy atom
            print(f"  Fallback (Joback) heat capacity for {fragment['name']}: {joback_cp:.1f} J/K·mol")
            return joback_cp
        
        # Last resort fallback based on molecular weight
        if mw > 0:
            mw_fallback = 1.8 * mw
            print(f"  Fallback (MW-based) heat capacity for {fragment['name']}: {mw_fallback:.1f} J/K·mol")
            return mw_fallback
            
        return None
        
    except Exception as e:
        print(f"Error calculating fragment heat capacity: {e}")
        return None

def calculate_ionic_liquid_heat_capacity(ionic_liquid: Dict, temperature: float = 298.15) -> Optional[Dict[str, float]]:
    """
    Calculate molar and gravimetric heat capacity for a complete ionic liquid.
    Args:
        ionic_liquid: Dictionary containing cation, anion, alkyl_chain,
                      and optional functional_group information. Each should
                      have at least 'name' and 'fragment_type'.
        temperature: Temperature in K (default: 298.15K)
    Returns:
        A dictionary containing 'molar' (J/K·mol) and 'gravimetric' (J/K·g)
        heat capacities, or None if calculation fails.
    """
    try:
        cation_name = ionic_liquid.get('cation', {}).get('name', 'Unknown Cation')
        anion_name = ionic_liquid.get('anion', {}).get('name', 'Unknown Anion')
        alkyl_name = ionic_liquid.get('alkyl_chain', {}).get('name', 'Unknown Alkyl')
        fg_name = ionic_liquid.get('functional_group', {}).get('name') if ionic_liquid.get('functional_group') else None

        print(f"\nCalculating heat capacity for {cation_name} {anion_name} with {alkyl_name}" + (f" and {fg_name}" if fg_name else "") + ":")

        # --- Get Component Properties (Cp and MW) ---
        components = ['cation', 'anion', 'alkyl_chain']
        if fg_name:
            components.append('functional_group')

        component_cps = {}
        component_mws = {}
        missing_data = False

        for comp_type in components:
            fragment = ionic_liquid.get(comp_type)
            if not fragment:
                print(f"Error: Missing component data for {comp_type}")
                missing_data = True
                continue

            fragment['fragment_type'] = comp_type # Ensure type is set

            # Get Molar Cp using the forward design calculation method
            cp = calculate_fragment_forward_cp(fragment, temperature)
            component_cps[comp_type] = cp
            if cp is None:
                print(f"Error: Could not calculate heat capacity for {comp_type} '{fragment.get('name', 'N/A')}'")
                missing_data = True

            # Get Molecular Weight
            props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
            mw = props.get('molecular_weight') if props else None
            component_mws[comp_type] = mw
            if mw is None:
                print(f"Error: Could not retrieve molecular weight for {comp_type} '{fragment.get('name', 'N/A')}'")
                missing_data = True

        if missing_data:
            print("Error: Missing critical data for one or more components. Cannot calculate total heat capacity.")
            return None

        # --- Calculate Total Molar Cp ---
        cation_cp = component_cps.get('cation', 0)
        anion_cp = component_cps.get('anion', 0)
        alkyl_cp = component_cps.get('alkyl_chain', 0)
        functional_group_cp = component_cps.get('functional_group', 0)

        # Fallback for functional group Cp if calculation failed but MW was found
        if 'functional_group' in component_cps and component_cps['functional_group'] is None and component_mws.get('functional_group') is not None:
            fg_frag = ionic_liquid['functional_group']
            group_name = fg_frag.get('name', '').lower()
            
            # Map functional group names to UNIFAC parameters
            if 'hydroxyl' in group_name:
                functional_group_cp = unifac_contribution('OH', temperature)
            elif 'amino' in group_name:
                functional_group_cp = unifac_contribution('NH2', temperature)
            elif 'carboxyl' in group_name:
                functional_group_cp = unifac_contribution('COOH', temperature)
            elif 'nitrile' in group_name or 'cyano' in group_name:
                functional_group_cp = unifac_contribution('CN', temperature)
            elif 'thiol' in group_name or 'mercapto' in group_name:
                functional_group_cp = unifac_contribution('SH', temperature)
            elif 'fluoro' in group_name:
                functional_group_cp = unifac_contribution('F', temperature)
            elif 'chloro' in group_name:
                functional_group_cp = unifac_contribution('Cl', temperature)
            elif 'bromo' in group_name:
                functional_group_cp = unifac_contribution('Br', temperature)
            elif 'iodo' in group_name:
                functional_group_cp = unifac_contribution('I', temperature)
            elif 'ether' in group_name or 'methoxy' in group_name or 'ethoxy' in group_name:
                functional_group_cp = unifac_contribution('Ether', temperature)
            else:
                functional_group_cp = 20.0  # Default value
            component_cps['functional_group'] = functional_group_cp # Update the value
            print(f"  Using estimated functional group contribution: {functional_group_cp:.1f} J/K·mol")

        print(f"\nComponent Molar Heat Capacities (J/K·mol):")
        for comp_type, cp in component_cps.items():
             print(f"  {comp_type.capitalize()}: {f'{cp:.1f}' if cp is not None else 'N/A'}")

        # Summing components for base molar Cp
        base_molar_cp = sum(cp for cp in component_cps.values() if cp is not None)
        print(f"  Base sum of component contributions: {base_molar_cp:.1f} J/K·mol")

        # --- Add Interaction Terms ---
        # Retrieve H-bond donors/acceptors
        cation_props = get_fragment_properties(ionic_liquid['cation']['name'], 'cation') or {}
        anion_props = get_fragment_properties(ionic_liquid['anion']['name'], 'anion') or {}
        fg_props = {}
        if fg_name:
             fg_props = get_fragment_properties(ionic_liquid['functional_group']['name'], 'functional_group') or {}

        cation_h_donors = cation_props.get('hydrogen_bond_donor_count', 0)
        cation_h_acceptors = cation_props.get('hydrogen_bond_acceptor_count', 0)
        anion_h_donors = anion_props.get('hydrogen_bond_donor_count', 0)
        anion_h_acceptors = anion_props.get('hydrogen_bond_acceptor_count', 0)
        functional_group_h_donors = fg_props.get('hydrogen_bond_donor_count', 0)
        functional_group_h_acceptors = fg_props.get('hydrogen_bond_acceptor_count', 0)

        # Fallback for functional group H-bonding if not in props
        if fg_name and functional_group_h_donors == 0 and functional_group_h_acceptors == 0:
             group_name = ionic_liquid['functional_group'].get('name', '').lower()
             if 'hydroxyl' in group_name:
                 functional_group_h_donors = 1
                 functional_group_h_acceptors = 1
             elif 'amino' in group_name:
                 functional_group_h_donors = 2
                 functional_group_h_acceptors = 1
             # Add other functional groups if needed

        # 1. Enhanced H-bond contribution with cooperative effects
        h_bond_pairs = (
            min(cation_h_donors, anion_h_acceptors) +
            min(anion_h_donors, cation_h_acceptors) +
            min(functional_group_h_donors, anion_h_acceptors) + # Interaction with anion
            min(anion_h_donors, functional_group_h_acceptors) + # Interaction with anion
            min(functional_group_h_donors, cation_h_acceptors) + # Interaction with cation
            min(cation_h_donors, functional_group_h_acceptors)  # Interaction with cation
        )
        
        # Enhanced H-bond contribution with cooperative effects
        # Higher coefficient to account for network effects
        h_bond_contribution = 6.0 * h_bond_pairs
        
        # Add cooperative enhancement for multiple H-bond pairs
        if h_bond_pairs > 1:
            # Cooperative effect increases with number of H-bonds
            cooperative_factor = 1.0 + 0.2 * (h_bond_pairs - 1)
            h_bond_contribution *= cooperative_factor
            print(f"  H-bond cooperative factor: {cooperative_factor:.2f}")

        print(f"  H-bond pairs considered: {h_bond_pairs}")
        print(f"  H-bond contribution: {h_bond_contribution:.1f} J/K·mol")

        # 2. Add interfragment interaction contribution
        # This accounts for non-H-bond interactions between fragments
        # Based on the size and complexity of the fragments
        total_mw = sum(mw for mw in component_mws.values() if mw is not None)
        total_heavy_atoms = 0
        for comp_type in components:
            props = get_fragment_properties(ionic_liquid[comp_type]['name'], comp_type) or {}
            total_heavy_atoms += props.get('heavy_atom_count', 0)
        
        # Interfragment contribution based on system complexity
        # More complex systems have more potential for interactions
        interfragment_contribution = 0.15 * total_mw + 2.0 * total_heavy_atoms
        
        # Scale based on number of components
        interfragment_contribution *= (len(components) / 3.0)
        
        print(f"  Interfragment interaction contribution: {interfragment_contribution:.1f} J/K·mol")
        
        # 3. Add contribution for additional degrees of freedom
        # This accounts for vibrational, rotational modes not captured by fragment calculations
        # Based on system size and complexity
        dof_contribution = 0.1 * total_mw + 1.5 * total_heavy_atoms
        print(f"  Additional degrees of freedom contribution: {dof_contribution:.1f} J/K·mol")
        
        # 4. Add ionic interaction contribution
        # Specific to ionic liquids - accounts for coulombic interactions
        # Stronger for smaller ions (higher charge density)
        cation_mw = component_mws.get('cation', 0)
        anion_mw = component_mws.get('anion', 0)
        
        if cation_mw > 0 and anion_mw > 0:
            # Inverse relationship with size - smaller ions have stronger interactions
            ionic_contribution = 5000.0 / (cation_mw + anion_mw)
            ionic_contribution = min(ionic_contribution, 50.0)  # Cap at reasonable value
            print(f"  Ionic interaction contribution: {ionic_contribution:.1f} J/K·mol")
        else:
            ionic_contribution = 0.0
        
        # Calculate total molar heat capacity with all contributions
        total_molar_cp = base_molar_cp + h_bond_contribution + interfragment_contribution + dof_contribution + ionic_contribution
        
        print(f"\nHeat Capacity Contributions:")
        print(f"  Base component sum: {base_molar_cp:.1f} J/K·mol")
        print(f"  H-bond contribution: {h_bond_contribution:.1f} J/K·mol")
        print(f"  Interfragment contribution: {interfragment_contribution:.1f} J/K·mol")
        print(f"  Degrees of freedom contribution: {dof_contribution:.1f} J/K·mol")
        print(f"  Ionic interaction contribution: {ionic_contribution:.1f} J/K·mol")
        print(f"  Total: {total_molar_cp:.1f} J/K·mol")

        # --- Calculate Total Molar Mass ---
        print(f"\nComponent Molecular Weights (g/mol):")
        for comp_type, mw in component_mws.items():
             print(f"  {comp_type.capitalize()}: {f'{mw:.2f}' if mw is not None else 'N/A'}")
        print(f"  Total Molecular Weight: {total_mw:.2f} g/mol")

        # --- Calculate Gravimetric Cp ---
        gravimetric_cp = None
        if total_mw > 0:
            gravimetric_cp = total_molar_cp / total_mw
            print(f"\nFinal Heat Capacity:")
            print(f"  Molar: {total_molar_cp:.2f} J/K·mol")
            print(f"  Gravimetric: {gravimetric_cp:.3f} J/K·g")
        else:
            print("\nError: Cannot calculate gravimetric heat capacity due to zero or missing molecular weight.")
            print(f"  Molar Heat Capacity: {total_molar_cp:.2f} J/K·mol")

        return {
            'molar': total_molar_cp,
            'gravimetric': gravimetric_cp
        }

    except Exception as e:
        import traceback
        print(f"Error calculating ionic liquid heat capacity: {e}")
        traceback.print_exc()
        return None

def screen_fragments_by_heat_capacity(fragments_data: Dict[str, List[Dict]], target_range: tuple) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated heat capacity
    Args:
        fragments_data: Dictionary of fragments organized by type (cation, anion, alkyl_chain)
        target_range: Tuple of (min_capacity, max_capacity) in J/K·mol
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
                # Apply a margin (e.g., 40%) to the target range for coarse screening
                # Keep fragments if their estimated capacity is within the broadened range
                # This filters out only extreme outliers whose individual contribution is far from the target.
                screening_min = min_capacity * 0.60 # Allow 40% below min target
                screening_max = max_capacity * 1.40 # Allow 40% above max target

                if screening_min <= estimated_capacity <= screening_max:
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
            if 'functional_group' in il and il['functional_group']:
                 print(f"Functional group: {il['functional_group']['name']}")

            results = calculate_ionic_liquid_heat_capacity(il)
            if results:
                molar_cp = results.get('molar')
                gravimetric_cp = results.get('gravimetric')
                print(f"Calculated Molar Cp: {f'{molar_cp:.2f}' if molar_cp is not None else 'N/A'} J/K·mol")
                print(f"Calculated Gravimetric Cp: {f'{gravimetric_cp:.3f}' if gravimetric_cp is not None else 'N/A'} J/K·g")
            else:
                print("Calculation failed.")

    if fragments_data:
        # Test fragment screening
        print("\nTesting fragment screening:")
        target_range = (100, 400)  # J/K·mol, can be parameterized
        screened = screen_fragments_by_heat_capacity(fragments_data, target_range)

        for frag_type, fragments in screened.items():
            if fragments:  # Only print if there are fragments that passed screening
                print(f"\n{frag_type.capitalize()} fragments within heat capacity range {target_range}:")
                for frag in fragments:
                    print(f"  {frag['name']}: {frag['estimated_heat_capacity']:.2f} J/K·mol")

if __name__ == "__main__":
    print("\n=== Heat Capacity Calculation Module Test ===")
    test_heat_capacity_calculations()
