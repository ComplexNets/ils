import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from models.shortList_frag import fragments
from utils.utils import get_fragment_properties

# Group contribution parameters for toxicity calculation
toxicity_params = {
    'cation': {
        'imidazolium': {'base_tox': 0.45, 'weight_factor': 1.2},  # Base toxicity for imidazolium core
        'pyridinium': {'base_tox': 0.40, 'weight_factor': 1.1}    # Slightly less toxic than imidazolium
    },
    'anion': {
        'tetrafluoroborate': {'base_tox': 0.30, 'weight_factor': 1.0},
        'hexafluorophosphate': {'base_tox': 0.35, 'weight_factor': 1.1},
        'bis(trifluoromethylsulfonyl)imide': {'base_tox': 0.40, 'weight_factor': 1.2}
    },
    'alkyl_chain': {
        'base': {'base_tox': 0.20, 'weight_factor': 1.0},
        'ch2_increment': {'base_tox': 0.05, 'weight_factor': 1.1}  # Toxicity increases with chain length
    }
}

def calculate_ionic_liquid_toxicity(combination: Dict) -> Optional[float]:
    """
    Calculate the toxicity of an ionic liquid using group contribution method.
    Returns toxicity value between 0 and 1, where higher values indicate higher toxicity.
    
    Args:
        combination (Dict): Dictionary containing cation, anion, and alkyl_chain components
        
    Returns:
        Optional[float]: Calculated toxicity value or None if calculation fails
    """
    try:
        # Get component properties
        cation = next((f for f in fragments if f['name'] == combination['cation']), None)
        anion = next((f for f in fragments if f['name'] == combination['anion']), None)
        alkyl_chain = next((f for f in fragments if f['name'] == combination['alkyl_chain']), None)
        
        if not all([cation, anion, alkyl_chain]):
            print(f"Error: Missing component properties for toxicity calculation")
            return None
            
        # Calculate cation contribution
        cation_params = toxicity_params['cation'].get(cation['type'], {'base_tox': 0.3, 'weight_factor': 1.0})
        cation_tox = cation_params['base_tox'] * cation_params['weight_factor']
        
        # Calculate anion contribution
        anion_params = toxicity_params['anion'].get(anion['type'], {'base_tox': 0.25, 'weight_factor': 1.0})
        anion_tox = anion_params['base_tox'] * anion_params['weight_factor']
        
        # Calculate alkyl chain contribution
        chain_length = alkyl_chain.get('length', 1)
        chain_params = toxicity_params['alkyl_chain']
        chain_tox = (chain_params['base']['base_tox'] + 
                    chain_params['ch2_increment']['base_tox'] * (chain_length - 1))
        chain_tox *= chain_params['base']['weight_factor']
        
        # Combine contributions with interaction terms
        total_tox = (cation_tox * 0.4 +  # Cation typically has largest impact
                    anion_tox * 0.3 +    # Anion has moderate impact
                    chain_tox * 0.3)     # Chain length has significant impact
                    
        # Add synergistic effects
        interaction_factor = 1.0
        if chain_length > 4:  # Longer chains tend to increase overall toxicity
            interaction_factor += 0.1 * (chain_length - 4)
            
        total_tox *= interaction_factor
        
        # Normalize to 0-1 range
        total_tox = min(max(total_tox, 0.0), 1.0)
        
        return total_tox
        
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
            toxicity_values.append(tox)
    
    if not toxicity_values:
        return {'min': 0.0, 'max': 1.0}
        
    return {
        'min': min(toxicity_values),
        'max': max(toxicity_values)
    }
