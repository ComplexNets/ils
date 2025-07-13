"""
Viscosity calculation module for ionic liquids using QSPR (Quantitative Structure-Property Relationship).
Primary method: Molecular descriptor-based regression model for viscosity prediction.
Reference: Based on QSPR modeling techniques for ionic liquid viscosity estimation.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional, Tuple
import math
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from utils.utils import get_fragment_properties
from utils.rdkit_utils import get_rdkit_properties

# Constants for viscosity calculations
DEFAULT_TEMP = 298.15  # K (reference temperature, 25°C)
MIN_VISCOSITY = 0.0005  # Pa·s (minimum reasonable viscosity for ionic liquids)
MAX_VISCOSITY = 1000.0  # Pa·s (maximum reasonable viscosity for ionic liquids)
ACTIVATION_ENERGY = 30000  # J/mol (activation energy for ionic liquids)
GAS_CONSTANT = 8.314  # J/(mol·K) (universal gas constant)
CALIBRATION_FACTOR = 133.0  # Empirical calibration factor to match experimental data

# QSPR model coefficients derived from training data
QSPR_COEFFICIENTS = {
    'a': -3.0,      # Intercept term (adjusted to get base viscosity right)
    'b1': 0.005,    # MW coefficient (moderate impact of molecular weight)
    'b2': 0.1,      # logP coefficient (moderate)
    'b3': 0.01,     # TPSA coefficient (moderate)
    'b4': 0.2,      # HBD coefficient (moderate)
    'b5': 0.1,      # HBA coefficient (moderate)
    'b6': -0.05,    # NRB (Number of Rotatable Bonds) coefficient
}

def calculate_molecular_descriptors(smiles: str) -> Optional[Dict[str, float]]:
    """
    Calculate molecular descriptors for viscosity prediction using RDKit.
    
    Args:
        smiles: SMILES representation of the molecule
        
    Returns:
        Dictionary of molecular descriptors or None if calculation fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        # Calculate descriptors
        descriptors = {
            'mw': Descriptors.ExactMolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'rotatable': Descriptors.NumRotatableBonds(mol)
        }
        
        return descriptors
    except Exception as e:
        print(f"Error calculating molecular descriptors: {str(e)}")
        return None

def estimate_charge_interaction(cation_desc: Dict[str, float], anion_desc: Dict[str, float]) -> float:
    """
    Estimate the strength of charge interaction between cation and anion.
    
    Args:
        cation_desc: Molecular descriptors for cation
        anion_desc: Molecular descriptors for anion
        
    Returns:
        Estimated charge interaction factor
    """
    # Simple model: larger, more polar ions have stronger interactions
    cation_size = cation_desc['mw'] / 100.0  # Normalize by typical IL cation size
    anion_size = anion_desc['mw'] / 100.0
    
    cation_polarity = cation_desc['tpsa'] / 100.0  # Normalize by typical IL TPSA
    anion_polarity = anion_desc['tpsa'] / 100.0
    
    # Combine size and polarity effects
    interaction = (cation_size * anion_size) * (cation_polarity + anion_polarity)
    
    return min(1.0, interaction)  # Cap at 1.0 to prevent unrealistic values

def calculate_viscosity(cation: Dict, anion: Dict, alkyl_chain: Dict, functional_group: Dict = None, 
                     temperature: float = DEFAULT_TEMP) -> Optional[float]:
    """
    Calculate viscosity of an ionic liquid using group contribution method
    
    Args:
        cation: Dictionary containing cation information
        anion: Dictionary containing anion information
        alkyl_chain: Dictionary containing alkyl chain information
        functional_group: Dictionary containing functional group information (optional)
        temperature: Temperature in K (default: 298.15 K)
        
    Returns:
        Viscosity in Pa·s or None if calculation fails
    """
    try:
        # Get molecular descriptors for each component
        cation_desc = calculate_molecular_descriptors(cation.get('smiles', ''))
        anion_desc = calculate_molecular_descriptors(anion.get('smiles', ''))
        alkyl_desc = calculate_molecular_descriptors(alkyl_chain.get('smiles', '')) if alkyl_chain else None
        
        if not all([cation_desc, anion_desc]):
            return None
            
        # Combine descriptors
        total_mw = cation_desc['mw'] + anion_desc['mw']
        total_logp = cation_desc['logp'] + anion_desc['logp']
        total_tpsa = cation_desc['tpsa'] + anion_desc['tpsa']
        total_hbd = cation_desc['hbd'] + anion_desc['hbd']
        total_hba = cation_desc['hba'] + anion_desc['hba']
        total_rotatable = cation_desc['rotatable'] + anion_desc['rotatable']
        
        if alkyl_desc:
            total_mw += alkyl_desc['mw']
            total_logp += alkyl_desc['logp']
            total_tpsa += alkyl_desc['tpsa']
            total_hbd += alkyl_desc['hbd']
            total_hba += alkyl_desc['hba']
            total_rotatable += alkyl_desc['rotatable']
            
        if functional_group:
            functional_desc = calculate_molecular_descriptors(functional_group.get('smiles', ''))
            if functional_desc:
                total_mw += functional_desc['mw']
                total_logp += functional_desc['logp']
                total_tpsa += functional_desc['tpsa']
                total_hbd += functional_desc['hbd']
                total_hba += functional_desc['hba']
                total_rotatable += functional_desc['rotatable']
            
        # Calculate log(viscosity) using QSPR model formula
        # The QSPR model was originally trained for viscosity in cP at reference temperature (298.15K)
        log_viscosity = (
            QSPR_COEFFICIENTS['a'] +
            QSPR_COEFFICIENTS['b1'] * total_mw +
            QSPR_COEFFICIENTS['b2'] * total_logp +
            QSPR_COEFFICIENTS['b3'] * total_tpsa +
            QSPR_COEFFICIENTS['b4'] * total_hbd +
            QSPR_COEFFICIENTS['b5'] * total_hba +
            QSPR_COEFFICIENTS['b6'] * total_rotatable
        )
        
        # Convert log(viscosity) to viscosity in cP at reference temperature
        viscosity_cp_ref = math.exp(log_viscosity)
        
        # Apply temperature correction using Arrhenius equation
        # η(T) = η(T_ref) * exp[E_a/R * (1/T - 1/T_ref)]
        if temperature != DEFAULT_TEMP:
            temperature_factor = math.exp(
                (ACTIVATION_ENERGY / GAS_CONSTANT) * 
                ((1 / temperature) - (1 / DEFAULT_TEMP))
            )
            viscosity_cp = viscosity_cp_ref * temperature_factor
        else:
            viscosity_cp = viscosity_cp_ref
        
        # Apply calibration factor to match experimental data
        viscosity_cp = viscosity_cp * CALIBRATION_FACTOR
        
        # Convert from cP to Pa·s (1 cP = 0.001 Pa·s)
        viscosity = viscosity_cp * 0.001
        
        return viscosity
        
    except Exception as e:
        print(f"Error calculating viscosity: {str(e)}")
        return None

def validate_viscosity(viscosity: Optional[float], 
                      min_viscosity: float = MIN_VISCOSITY,
                      max_viscosity: float = MAX_VISCOSITY) -> bool:
    """
    Validate if the calculated viscosity is within reasonable bounds for ionic liquids.
    
    Args:
        viscosity: Calculated viscosity in Pa·s
        min_viscosity: Minimum acceptable viscosity (default: 0.0005 Pa·s)
        max_viscosity: Maximum acceptable viscosity (default: 1000.0 Pa·s)
        
    Returns:
        bool: True if viscosity is valid, False otherwise
    """
    if viscosity is None:
        return False
        
    try:
        viscosity_value = float(viscosity)
        return min_viscosity <= viscosity_value <= max_viscosity
    except (TypeError, ValueError):
        return False

def screen_fragments_by_viscosity(fragments_data: Dict[str, List[Dict]], 
                                target_range: Tuple[float, float]) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated viscosity contribution.
    
    Args:
        fragments_data: Dictionary of fragments organized by type
        target_range: Tuple of (min_viscosity, max_viscosity) in Pa·s
        
    Returns:
        Dictionary of fragments organized by type that meet the target range
    """
    min_viscosity, max_viscosity = target_range
    screened_fragments = {
        'cation': [],
        'anion': [],
        'alkyl_chain': []
    }
    
    # Create a minimal test IL for each fragment to estimate its contribution
    for frag_type, frags in fragments_data.items():
        for frag in frags:
            try:
                # Create a test IL using the fragment and minimal counter-ions
                test_il = {
                    'cation': {'smiles': 'C[N+](C)(C)C'} if frag_type != 'cation' else frag,  # Use TMA+ as test cation
                    'anion': {'smiles': '[BF4-]'} if frag_type != 'anion' else frag,  # Use BF4- as test anion
                    'alkyl_chain': {'smiles': 'CC'} if frag_type != 'alkyl_chain' else frag  # Use ethyl as test chain
                }
                
                viscosity = calculate_viscosity(
                    test_il['cation'],
                    test_il['anion'],
                    test_il['alkyl_chain']
                )

                # Apply a margin (e.g., 40%) to the target range for coarse screening
                screening_min = min_viscosity * 0.60 # Allow 40% below min target
                screening_max = max_viscosity * 1.40 # Allow 40% above max target

                # Check if the viscosity of the test combination is within the broadened range
                if viscosity is not None and screening_min <= viscosity <= screening_max:
                    screened_fragments[frag_type].append(frag)

            except Exception as e: # Catch specific exceptions if possible, or log the error
                print(f"Screening error for fragment {frag.get('name', 'Unknown')}: {e}")
                continue
                
    return screened_fragments

def test_viscosity_calculations(fragments_data: Optional[Dict[str, List[Dict]]] = None,
                              combinations: Optional[List[Dict]] = None,
                              num_test_combinations: int = 3):
    """
    Test function to show detailed viscosity calculations.
    
    Args:
        fragments_data: Optional dictionary of fragments to test screening
        combinations: Optional list of ionic liquid combinations to test
        num_test_combinations: Number of combinations to test
    """
    if fragments_data:
        print("\nTesting fragment screening...")
        target_range = (0.0005, 1000.0)  # Example viscosity range in Pa·s
        screened = screen_fragments_by_viscosity(fragments_data, target_range)
        for frag_type, frags in screened.items():
            print(f"\n{frag_type}: {len(frags)} fragments in range {target_range}")
            
    if combinations:
        print("\nTesting viscosity calculations...")
        for i, combo in enumerate(combinations[:num_test_combinations]):
            viscosity = calculate_viscosity(
                combo['cation'],
                combo['anion'],
                combo['alkyl_chain']
            )
            print(f"\nCombination {i+1}:")
            print(f"Name: {combo.get('name', 'Unknown')}")
            print(f"Viscosity: {viscosity:.4f} Pa·s")
            print(f"Valid: {validate_viscosity(viscosity)}")

if __name__ == "__main__":
    print("\n=== Viscosity Calculation Module Test ===")
    test_viscosity_calculations()
