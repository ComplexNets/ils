"""
Hydrophobicity calculation module for ionic liquids using RDKit's Crippen method.
Primary method: Fragment-based logP calculation.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, Optional, List, Tuple
from rdkit import Chem
from rdkit.Chem import Crippen
from models.shortList_frag import fragments
from utils.utils import get_fragment_properties
from utils.rdkit_utils import get_rdkit_properties

# Constants for hydrophobicity calculations
DEFAULT_TEMP = 298.15  # K (reference temperature)
MIN_LOGP = -10.0  # Expanded minimum reasonable logP value
MAX_LOGP = 15.0  # Expanded maximum reasonable logP value

def calculate_fragment_logp(smiles: str) -> float:
    """
    Calculate the logP contribution of a fragment using RDKit's Crippen method.
    
    Args:
        smiles (str): SMILES representation of the fragment
        
    Returns:
        float: Calculated logP value for the fragment
        
    Raises:
        ValueError: If SMILES string cannot be parsed or logP calculation fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not parse SMILES string: {smiles}")
        
        return Crippen.MolLogP(mol)
    except Exception as e:
        raise ValueError(f"Error calculating logP for {smiles}: {str(e)}")

def calculate_ionic_liquid_hydrophobicity(combination: Dict) -> Optional[float]:
    """
    Calculate the hydrophobicity (logP) of an ionic liquid.
    
    Args:
        combination: Dictionary containing cation, anion, and alkyl_chain properties
        
    Returns:
        float: Calculated logP value, or None if calculation fails
    """
    try:
        cation = combination.get('cation', {})
        anion = combination.get('anion', {})
        alkyl_chain = combination.get('alkyl_chain', {})
        temperature = DEFAULT_TEMP
        
        if not all([cation, anion, alkyl_chain]):
            print(f"Missing component in combination: {combination}")
            return None
            
        # Get RDKit properties for each component
        cation_smiles = cation.get('smiles', '')
        anion_smiles = anion.get('smiles', '')
        alkyl_smiles = alkyl_chain.get('smiles', '')
        
        print(f"Processing SMILES - Cation: {cation_smiles}, Anion: {anion_smiles}, Alkyl: {alkyl_smiles}")
        
        cation_props = get_rdkit_properties(cation_smiles)
        anion_props = get_rdkit_properties(anion_smiles)
        alkyl_props = get_rdkit_properties(alkyl_smiles)
        
        if not all([cation_props, anion_props, alkyl_props]):
            print(f"Failed to get RDKit properties for one or more components:")
            print(f"Cation props: {cation_props}")
            print(f"Anion props: {anion_props}")
            print(f"Alkyl props: {alkyl_props}")
            return None
            
        # Extract logP values
        cation_logp = cation_props.get('MolLogP', 0.0)
        anion_logp = anion_props.get('MolLogP', 0.0)
        alkyl_logp = alkyl_props.get('MolLogP', 0.0)
        
        # Add functional group contribution if present
        functional_group_logp = 0.0
        if 'functional_group' in combination and combination['functional_group']:
            functional_group = combination['functional_group']
            functional_group_smiles = functional_group.get('smiles', '')
            
            if functional_group_smiles:
                print(f"Processing functional group SMILES: {functional_group_smiles}")
                functional_group_props = get_rdkit_properties(functional_group_smiles)
                
                if functional_group_props:
                    functional_group_logp = functional_group_props.get('MolLogP', 0.0)
                    print(f"Functional group LogP: {functional_group_logp}")
                else:
                    print(f"Failed to get RDKit properties for functional group: {functional_group_smiles}")
                    # Estimate based on functional group type
                    group_name = functional_group.get('name', '').lower()
                    if 'hydroxyl' in group_name or 'hydroxy' in group_name:
                        functional_group_logp = -1.0  # Hydroxyl groups decrease logP
                    elif 'amino' in group_name:
                        functional_group_logp = -0.5  # Amino groups slightly decrease logP
                    elif 'carboxyl' in group_name:
                        functional_group_logp = -1.5  # Carboxyl groups significantly decrease logP
                    elif 'fluoro' in group_name:
                        functional_group_logp = 0.5   # Fluoro groups increase logP
                    elif 'chloro' in group_name:
                        functional_group_logp = 1.0   # Chloro groups increase logP
                    elif 'bromo' in group_name:
                        functional_group_logp = 1.5   # Bromo groups increase logP
                    print(f"Estimated functional group LogP: {functional_group_logp}")
        
        print(f"LogP values - Cation: {cation_logp}, Anion: {anion_logp}, Alkyl: {alkyl_logp}, Functional Group: {functional_group_logp}")
        
        # Improved temperature correction with more accurate coefficients
        # Based on experimental data showing logP decreases with temperature for most ionic liquids
        # Temperature coefficient is now dependent on the base logP value
        base_logp = cation_logp + anion_logp + alkyl_logp + functional_group_logp
        
        # Different temperature correction factors for hydrophilic vs hydrophobic ILs
        if base_logp < 0:
            # Hydrophilic ILs - temperature has less effect
            temp_factor = 1.0 + 0.0015 * (temperature - DEFAULT_TEMP)
        else:
            # Hydrophobic ILs - temperature has more effect
            temp_factor = 1.0 + 0.0025 * (temperature - DEFAULT_TEMP)
        
        # Calculate total logP with temperature correction
        total_logp = base_logp * temp_factor
        
        print(f"Total LogP (before validation): {total_logp}")
        
        # Validate result
        if not validate_hydrophobicity(total_logp):
            print(f"Hydrophobicity validation failed for value: {total_logp}")
            return None
            
        return total_logp
        
    except Exception as e:
        print(f"Error in hydrophobicity calculation: {str(e)}")
        return None

def validate_hydrophobicity(logp: Optional[float], 
                          min_logp: float = MIN_LOGP,
                          max_logp: float = MAX_LOGP) -> bool:
    """
    Validate if the calculated logP is within reasonable bounds for ionic liquids.
    
    Args:
        logp: Calculated logP value
        min_logp: Minimum acceptable logP value (default: -10.0)
        max_logp: Maximum acceptable logP value (default: 15.0)
        
    Returns:
        bool: True if logP is valid, False otherwise
    """
    if logp is None:
        return False
        
    try:
        logp_value = float(logp)
        return min_logp <= logp_value <= max_logp
    except (TypeError, ValueError):
        return False

def screen_fragments_by_hydrophobicity(fragments_data: Dict[str, List[Dict]], 
                                     target_range: Tuple[float, float]) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated hydrophobicity.
    
    Args:
        fragments_data: Dictionary of fragments organized by type
        target_range: Tuple of (min_logp, max_logp)
        
    Returns:
        Dictionary of fragments organized by type that meet the target range
    """
    min_logp, max_logp = target_range
    screened_fragments = {
        'cation': [],
        'anion': [],
        'alkyl_chain': []
    }
    
    for frag_type, frags in fragments_data.items():
        for frag in frags:
            if 'smiles' not in frag:
                continue
                
            try:
                logp = calculate_fragment_logp(frag['smiles'])

                # Apply a margin (e.g., +/- 1.0 log unit) to the target range for coarse screening
                screening_min = min_logp - 1.0
                screening_max = max_logp + 1.0

                if screening_min <= logp <= screening_max:
                    screened_fragments[frag_type].append(frag)
            except ValueError: # Be more specific about the exception
                continue
                
    return screened_fragments

def test_hydrophobicity_calculations(fragments_data: Optional[Dict[str, List[Dict]]] = None,
                                   combinations: Optional[List[Dict]] = None,
                                   num_test_combinations: int = 3):
    """
    Test function to show detailed hydrophobicity calculations.
    
    Args:
        fragments_data: Optional dictionary of fragments to test screening
        combinations: Optional list of ionic liquid combinations to test
        num_test_combinations: Number of combinations to test
    """
    if fragments_data:
        print("\nTesting fragment screening...")
        target_range = (-2.0, 5.0)
        screened = screen_fragments_by_hydrophobicity(fragments_data, target_range)
        for frag_type, frags in screened.items():
            print(f"\n{frag_type}: {len(frags)} fragments in range {target_range}")
            
    if combinations:
        print("\nTesting hydrophobicity calculations...")
        for i, combo in enumerate(combinations[:num_test_combinations]):
            logp = calculate_ionic_liquid_hydrophobicity(combo)
            print(f"\nCombination {i+1}:")
            print(f"Name: {combo.get('name', 'Unknown')}")
            print(f"LogP: {logp:.2f}")
            print(f"Valid: {validate_hydrophobicity(logp)}")

if __name__ == "__main__":
    print("\n=== Hydrophobicity Calculation Module Test ===")
    test_hydrophobicity_calculations()
