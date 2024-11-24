import os
import pymysql
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Get database credentials from the .env file
DB_HOST = os.getenv("DB_HOST")
DB_USER = os.getenv("DB_USER")
DB_PASSWORD = os.getenv("DB_PASSWORD")
DB_NAME = os.getenv("DB_NAME")

def print_table_structure(cursor, table_name):
    """Print the structure of a table"""
    print(f"\n=== {table_name} Structure ===")
    cursor.execute(f"DESCRIBE {table_name}")
    columns = cursor.fetchall()
    for col in columns:
        print(f"  - {col['Field']}: {col['Type']} ({col['Null']} null) {col['Key']} key")

def print_fragment_data(cursor, table_name):
    """Print fragment data grouped by type"""
    # Get all fragments
    cursor.execute(f"SELECT * FROM {table_name}")
    fragments = cursor.fetchall()
    
    if not fragments:
        print(f"No data found in {table_name}")
        return
        
    # Group by type
    fragments_by_type = {}
    for frag in fragments:
        frag_type = frag.get('Type', 'Unknown')
        if frag_type not in fragments_by_type:
            fragments_by_type[frag_type] = []
        fragments_by_type[frag_type].append(frag)
    
    # Print fragments by type
    for frag_type, frags in fragments_by_type.items():
        print(f"\n{frag_type}s ({len(frags)}):")
        for frag in frags:
            print(f"\n  {frag['Name']}:")
            for key, val in frag.items():
                if key != 'Name':
                    print(f"    {key}: {val}")

try:
    # Connect to the database
    connection = pymysql.connect(
        host=DB_HOST,
        user=DB_USER,
        password=DB_PASSWORD,
        database=DB_NAME,
        cursorclass=pymysql.cursors.DictCursor
    )

    with connection.cursor() as cursor:
        # Fetch all tables in the database
        cursor.execute("SHOW TABLES")
        tables = cursor.fetchall()

        print("\n=== Database Tables ===")
        for table_entry in tables:
            table_name = list(table_entry.values())[0]  # Get the table name
            print(f"- {table_name}")
            
            # Print table structure
            print_table_structure(cursor, table_name)
            
            # Print fragment data
            print_fragment_data(cursor, table_name)
            print("\n" + "="*50)

except pymysql.MySQLError as e:
    print(f"Error connecting to the database: {e}")
finally:
    if connection:
        connection.close()
