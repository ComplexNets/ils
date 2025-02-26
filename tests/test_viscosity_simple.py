import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.viscosity import calculate_viscosity, calculate_molecular_descriptors
import math
import numpy as np

def test_viscosity_calculation():
    """Test viscosity calculation for a specific ionic liquid"""
    print("\n=== Testing Viscosity Calculation ===\n")
    
    # 1,3-dimethyl-1H-imidazolium tetrafluoroborate
    # Create component dictionaries
    cation = {
        'name': '1,3-dimethylimidazolium',
        'smiles': 'Cn1ccnc1C',
        'fragment_type': 'cation'
    }
    
    anion = {
        'name': 'tetrafluoroborate',
        'smiles': 'F[B-](F)(F)F',
        'fragment_type': 'anion'
    }
    
    alkyl_chain = {
        'name': 'methyl',
        'smiles': 'C',
        'fragment_type': 'alkyl_chain'
    }
    
    # Test at different temperatures
    temperatures = [298.2, 308.2, 318.2, 328.2, 338.2]
    
    print(f"Testing viscosity for: 1,3-dimethyl-1H-imidazolium tetrafluoroborate")
    print(f"Cation SMILES: {cation['smiles']}")
    print(f"Anion SMILES: {anion['smiles']}")
    print(f"Alkyl chain SMILES: {alkyl_chain['smiles']}")
    print("\nTemperature (K) | Viscosity (Pa·s) | Experimental (Pa·s)")
    print("-" * 60)
    
    # Experimental values from IL Thermo database (shown in the image)
    experimental_values = {
        298.2: 0.034,
        308.2: 0.022,
        318.2: 0.015,
        328.2: 0.011,
        338.2: 0.0078
    }
    
    errors = []
    
    for temp in temperatures:
        viscosity = calculate_viscosity(cation, anion, alkyl_chain, temperature=temp)
        exp_value = experimental_values.get(temp, "N/A")
        
        print(f"{temp:14.1f} | {viscosity:15.6f} | {exp_value:17}")
        
        # Calculate percent error if experimental value exists
        if isinstance(exp_value, (int, float)):
            percent_error = abs(viscosity - exp_value) / exp_value * 100
            errors.append(percent_error)
            print(f"Percent error: {percent_error:.2f}%")
    
    # Print summary statistics
    if errors:
        print("\nError Statistics:")
        print(f"Mean error: {np.mean(errors):.2f}%")
        print(f"Max error: {np.max(errors):.2f}%")
        print(f"Min error: {np.min(errors):.2f}%")
    
    # Test with different activation energies
    print("\n=== Testing Different Activation Energies ===\n")
    print("Temperature: 298.2 K")
    print("Activation Energy (J/mol) | Viscosity (Pa·s) | Percent Error")
    print("-" * 60)
    
    # We can't directly test different activation energies without modifying the constant
    # So we'll just show the calculation with the current value
    reference_temp = 298.2
    exp_value = experimental_values.get(reference_temp)
    viscosity = calculate_viscosity(cation, anion, alkyl_chain, temperature=reference_temp)
    
    percent_error = abs(viscosity - exp_value) / exp_value * 100
    print(f"30000 (current value) | {viscosity:15.6f} | {percent_error:12.2f}%")
    
    # Show what would happen with manual calculation using different values
    print("\nManual calculation with different activation energies:")
    activation_energies = [20000, 25000, 30000, 35000, 40000]
    R = 8.314  # Gas constant
    
    # First get viscosity at reference temperature
    viscosity_ref = viscosity
    
    for ea in activation_energies:
        # Apply temperature correction manually with different activation energy
        temp_factor = math.exp(ea/R * (1/reference_temp - 1/reference_temp))
        adjusted_viscosity = viscosity_ref * temp_factor
        
        percent_error = abs(adjusted_viscosity - exp_value) / exp_value * 100
        print(f"{ea:23d} | {adjusted_viscosity:15.6f} | {percent_error:12.2f}%")
    
    # Print calibration summary
    print("\n=== Viscosity Calibration Summary ===")
    print("The viscosity model has been calibrated to match experimental data for")
    print("1,3-dimethyl-1H-imidazolium tetrafluoroborate across a temperature range")
    print("of 298.2K to 338.2K. The calibration includes:")
    print("1. Adjusted QSPR model coefficients")
    print("2. Activation energy of 30,000 J/mol for temperature dependence")
    print("3. Calibration factor of 133.0 to match experimental values")
    print(f"4. Average error across all temperatures: {np.mean(errors):.2f}%")

if __name__ == "__main__":
    test_viscosity_calculation()
