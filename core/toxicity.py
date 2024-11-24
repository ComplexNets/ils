import sys
import os
import math
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from models.shortList_frag import fragments
from utils.utils import get_fragment_properties

# Group contribution parameters for toxicity calculation
toxicity_params = {
    'cation': {
        'Ethylimidazolium': {'base_tox': 0.35, 'weight_factor': 1.1},  # Adjusted for ~8-10 mM IC50
        '1-Butylpyridinium': {'base_tox': 0.32, 'weight_factor': 1.0}  # Adjusted for ~8-12 mM IC50
    },
    'anion': {
        'Tetrafluoroborate': {'base_tox': 0.25, 'weight_factor': 0.9},  # Less toxic
        'Bis(trifluoromethanesulfonyl)imide': {'base_tox': 0.30, 'weight_factor': 1.1}  # More toxic
    },
    'alkyl_chain': {
        'Ethyl': {'base_tox': 0.15, 'weight_factor': 0.9},  # Lower toxicity contribution
        'Butyl': {'base_tox': 0.20, 'weight_factor': 1.1}   # Higher toxicity for longer chain
    }
}

def calculate_logp_contribution(smiles: str) -> float:
    """
    Calculate normalized logP contribution for a molecule
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        float: Normalized logP contribution (0-1 range)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        # Calculate logP using RDKit's built-in method
        logp = Descriptors.MolLogP(mol)
        
        # Normalize logP to 0-1 range (typical logP values for ILs range from -2 to 4)
        norm_logp = (logp + 2) / 6
        return max(min(norm_logp, 1.0), 0.0)
        
    except Exception as e:
        print(f"Error calculating logP: {str(e)}")
        return 0.0

def normalize_to_ic50(toxicity_score: float) -> float:
    """
    Convert normalized toxicity score (0-1) to IC50 value in mM.
    Higher IC50 values indicate lower toxicity.
    
    Typical IC50 ranges for ionic liquids are roughly 0.1 mM to 100 mM.
    We use an exponential conversion to better represent the concentration scale.
    
    Args:
        toxicity_score (float): Normalized toxicity score (0-1)
        
    Returns:
        float: IC50 value in mM
    """
    # Convert toxicity score to IC50
    # Using exponential scale: 0.1 mM (very toxic) to 100 mM (less toxic)
    min_ic50 = 0.1  # mM
    max_ic50 = 100.0  # mM
    
    # Invert the toxicity score since higher IC50 means less toxic
    inverted_score = 1 - toxicity_score
    
    # Calculate IC50 using exponential scale
    log_min = math.log10(min_ic50)
    log_max = math.log10(max_ic50)
    log_ic50 = log_min + (log_max - log_min) * inverted_score
    
    return 10 ** log_ic50

def get_fragment_smiles(fragment_name: str, fragment_type: str) -> str:
    """Get SMILES string for a given fragment name and type."""
    for fragment in fragments:
        if fragment['name'] == fragment_name and fragment['fragment_type'] == fragment_type:
            return fragment['smiles']
    return ''

def calculate_ionic_liquid_toxicity(combination: Dict) -> Optional[Dict]:
    """
    Calculate the toxicity of an ionic liquid using group contribution method and logP.
    Returns IC50 value in mM, where higher values indicate lower toxicity.
    
    Args:
        combination (Dict): Dictionary containing cation, anion, and alkyl_chain components
        
    Returns:
        Optional[Dict]: Dictionary containing IC50 value and component contributions
    """
    try:
        # Get component names
        cation_name = combination['cation']['name']
        anion_name = combination['anion']['name']
        alkyl_name = combination['alkyl_chain']['name']
        
        # Calculate base toxicity for each component
        components = {}
        total_tox = 0.0
        total_weight = 0.0
        
        # Cation contribution
        if cation_name in toxicity_params['cation']:
            cation_params = toxicity_params['cation'][cation_name]
            cation_smiles = get_fragment_smiles(cation_name, 'cation')
            logp_contrib = calculate_logp_contribution(cation_smiles) if cation_smiles else 0.0
            components['cation'] = {
                'base_tox': cation_params['base_tox'],
                'logp': logp_contrib
            }
            total_tox += cation_params['base_tox'] * cation_params['weight_factor']
            total_weight += cation_params['weight_factor']
        
        # Anion contribution
        if anion_name in toxicity_params['anion']:
            anion_params = toxicity_params['anion'][anion_name]
            anion_smiles = get_fragment_smiles(anion_name, 'anion')
            logp_contrib = calculate_logp_contribution(anion_smiles) if anion_smiles else 0.0
            components['anion'] = {
                'base_tox': anion_params['base_tox'],
                'logp': logp_contrib
            }
            total_tox += anion_params['base_tox'] * anion_params['weight_factor']
            total_weight += anion_params['weight_factor']
        
        # Alkyl chain contribution
        if alkyl_name in toxicity_params['alkyl_chain']:
            alkyl_params = toxicity_params['alkyl_chain'][alkyl_name]
            alkyl_smiles = get_fragment_smiles(alkyl_name, 'alkyl_chain')
            logp_contrib = calculate_logp_contribution(alkyl_smiles) if alkyl_smiles else 0.0
            components['alkyl_chain'] = {
                'base_tox': alkyl_params['base_tox'],
                'logp': logp_contrib
            }
            total_tox += alkyl_params['base_tox'] * alkyl_params['weight_factor']
            total_weight += alkyl_params['weight_factor']
        
        if total_weight == 0:
            print("Error: No valid components found for toxicity calculation")
            return None
            
        # Calculate average toxicity score (0-1 scale, higher is more toxic)
        avg_tox = total_tox / total_weight
        
        # Convert to IC50 value (higher is less toxic)
        ic50 = normalize_to_ic50(avg_tox)
        
        return {
            'ic50_mm': ic50,
            'components': components
        }
        
    except Exception as e:
        print(f"Error calculating toxicity: {str(e)}")
        return None

def get_toxicity_range(combinations: List[Dict]) -> Dict[str, float]:
    """
    Calculate the range of toxicity values for a list of combinations.
    
    Args:
        combinations (List[Dict]): List of ionic liquid combinations
        
    Returns:
        Dict[str, float]: Dictionary containing min and max toxicity values
    """
    toxicity_values = []
    
    for combo in combinations:
        tox = calculate_ionic_liquid_toxicity(combo)
        if tox is not None:
            toxicity_values.append(tox['ic50_mm'])
    
    if not toxicity_values:
        return {'min': 0.0, 'max': 1.0}
        
    return {
        'min': min(toxicity_values),
        'max': max(toxicity_values)
    }

def screen_fragments_by_toxicity(fragments_data: Dict[str, List[Dict]], target_range: tuple) -> Dict[str, List[Dict]]:
    """
    Screen fragments based on estimated IC50 toxicity values.
    
    Args:
        fragments_data (Dict): Dictionary of fragments organized by type (cation, anion, alkyl_chain)
        target_range (tuple): Tuple of (min_ic50, max_ic50) in mM
    Returns:
        Dict[str, List[Dict]]: Dictionary of fragments that meet criteria, organized by type
    """
    min_ic50, max_ic50 = target_range
    
    # Initialize output dictionary with same structure as input
    screened_fragments = {frag_type: [] for frag_type in fragments_data.keys()}
    
    for frag_type, fragments in fragments_data.items():
        for fragment in fragments:
            # Create a simple test combination with just this fragment
            test_combo = {
                'cation': 'Ethylimidazolium',  # Default reference cation
                'anion': 'Tetrafluoroborate',   # Default reference anion
                'alkyl_chain': 'Ethyl'          # Default reference alkyl chain
            }
            
            # Replace the corresponding component with our test fragment
            test_combo[frag_type] = fragment['name']
            
            # Calculate toxicity for this combination
            result = calculate_ionic_liquid_toxicity(test_combo)
            if result is not None:
                ic50 = result['ic50_mm']
                
                # Use configurable margin for screening
                margin = 0.2  # Can be adjusted or passed as parameter
                if (min_ic50 * (1 - margin) <= ic50 <= max_ic50 * (1 + margin)):
                    screened_fragments[frag_type].append({
                        **fragment,  # Preserve all original properties
                        'estimated_ic50': ic50,
                        'toxicity_components': result['components'][frag_type]
                    })
    
    return screened_fragments

def test_toxicity_calculations(fragments_data: Optional[Dict[str, List[Dict]]] = None,
                             combinations: Optional[List[Dict]] = None,
                             num_test_combinations: int = 3):
    """
    Test function to show detailed toxicity calculations
    
    Args:
        fragments_data: Optional dictionary of fragments to test screening
        combinations: Optional list of ionic liquid combinations to test
        num_test_combinations: Number of combinations to test
    """
    print("\n=== Testing Toxicity Calculations ===")
    
    if combinations:
        # Test provided combinations
        for i, il in enumerate(combinations[:num_test_combinations]):
            print(f"\nTesting combination {i+1}:")
            print(f"Cation: {il['cation']['name']}")
            print(f"Anion: {il['anion']['name']}")
            if 'alkyl_chain' in il:
                print(f"Alkyl chain: {il['alkyl_chain']['name']}")
            
            simple_combo = {
                'cation': il['cation']['name'],
                'anion': il['anion']['name'],
                'alkyl_chain': il['alkyl_chain']['name']
            }
            
            result = calculate_ionic_liquid_toxicity(simple_combo)
            if result:
                print(f"Calculated IC50: {result['ic50_mm']:.2f} mM")
                print("Component contributions:")
                for component, values in result['components'].items():
                    print(f"  {component}: base_tox={values['base_tox']:.3f}, logP={values['logp']:.3f}")
    
    if fragments_data:
        # Test fragment screening
        print("\nTesting fragment screening:")
        target_range = (5.0, 100.0)  # mM, can be parameterized
        screened = screen_fragments_by_toxicity(fragments_data, target_range)
        
        for frag_type, fragments in screened.items():
            if fragments:  # Only print if there are fragments that passed screening
                print(f"\n{frag_type.capitalize()} fragments within IC50 range {target_range}:")
                for frag in fragments:
                    print(f"  {frag['name']}: {frag['estimated_ic50']:.2f} mM")
                    print(f"    Base toxicity: {frag['toxicity_components']['base_tox']:.3f}")
                    print(f"    LogP contribution: {frag['toxicity_components']['logp']:.3f}")
    
    # Print summary of all valid combinations
    print("\nSummary of All Valid Combinations:")
    toxicity_values = []
    for combo in combinations:
        simple_combo = {
            'cation': combo['cation']['name'],
            'anion': combo['anion']['name'],
            'alkyl_chain': combo['alkyl_chain']['name']
        }
        result = calculate_ionic_liquid_toxicity(simple_combo)
        if result:
            toxicity_values.append((result['ic50_mm'], combo))
    
    if toxicity_values:
        print(f"Number of valid combinations: {len(toxicity_values)}")
        ic50_values = [t[0] for t in toxicity_values]
        print(f"IC50 range: {min(ic50_values):.2f} - {max(ic50_values):.2f} mM")
        print(f"Average IC50: {sum(ic50_values)/len(ic50_values):.2f} mM")
        
        print("\nLeast toxic combinations (highest IC50, top 5):")
        sorted_combinations = sorted(toxicity_values, key=lambda x: x[0], reverse=True)
        for ic50, combo in sorted_combinations[:5]:
            print(f"  {combo['cation']['name']} + {combo['anion']['name']} + {combo['alkyl_chain']['name']}: {ic50:.2f} mM")

if __name__ == "__main__":
    print("\n=== Toxicity Calculation Module Test ===")
    # Create some test combinations
    test_combinations = [
        {
            'cation': {'name': 'EMIM'},
            'anion': {'name': 'BF4'},
            'alkyl_chain': {'name': 'ethyl'}
        },
        {
            'cation': {'name': 'BMIM'},
            'anion': {'name': 'PF6'},
            'alkyl_chain': {'name': 'butyl'}
        },
        {
            'cation': {'name': 'HMIM'},
            'anion': {'name': 'Cl'},
            'alkyl_chain': {'name': 'hexyl'}
        }
    ]
    test_toxicity_calculations(combinations=test_combinations)
