"""
Test script to check ILThermo database connection
"""
import requests
import json

def test_ilthermo_connection(il_name):
    """Test connection to ILThermo database"""
    print(f"Testing ILThermo connection for: {il_name}")
    
    base_url = "https://ilthermo.boulder.nist.gov/ILT2/ilsearch"
    
    # Standardize the name: ensure proper spacing and lowercase
    search_name = ' '.join(il_name.lower().split())
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'application/json'
    }
    
    try:
        params = {'cmp': search_name, 'orderby': 'T', 'output': 'json'}
        print(f"Request URL: {base_url}?cmp={search_name}&orderby=T&output=json")
        
        response = requests.get(base_url, params=params, headers=headers, timeout=10)
        print(f"Response status code: {response.status_code}")
        
        if response.status_code == 200:
            try:
                # Parse JSON response
                data = response.json()
                print(f"Response JSON: {json.dumps(data, indent=2)}")
                
                # Check if we got any results
                if 'res' in data and data['res']:
                    print(f"✅ Found '{search_name}' in ILThermo database!")
                    return True
                else:
                    print(f"❌ '{search_name}' not found in ILThermo database.")
                    return False
            except json.JSONDecodeError:
                print(f"❌ Error parsing JSON response: {response.text[:200]}...")
                return False
        else:
            print(f"❌ HTTP error: {response.status_code}")
            return False
            
    except Exception as e:
        print(f"❌ Error accessing IL Thermo: {e}")
        return False

if __name__ == "__main__":
    # Test with a known IL in ILThermo
    test_ilthermo_connection("1-ethyl-3-methylimidazolium tetrafluoroborate")
    
    # Test with a slightly different formatting
    test_ilthermo_connection("1-ethyl-3-methyl imidazolium tetrafluoroborate")
    
    # Test with a known IL not in our list
    test_ilthermo_connection("1-propyl-3-methylimidazolium tetrafluoroborate")
