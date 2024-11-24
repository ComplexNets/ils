#!/usr/bin/env python3
import os
import sys
from utils.utils import connect_to_database

def get_tables(cursor):
    """Get all table names from the database"""
    cursor.execute("SHOW TABLES")
    return [table[0] for table in cursor.fetchall()]

def get_table_structure(cursor, table_name):
    """Get table structure"""
    cursor.execute(f"DESCRIBE {table_name}")
    columns = cursor.fetchall()
    
    print(f"\n=== {table_name} Structure ===")
    print("Columns:")
    for col in columns:
        print(f"  - {col[0]}: {col[1]} ({col[2]}) {'PRIMARY KEY' if col[3] == 'PRI' else ''}")

def explore_table(cursor, table_name):
    """Explore table data"""
    # Get table structure
    get_table_structure(cursor, table_name)
    
    # Get all rows with all columns
    cursor.execute(f"SELECT * FROM {table_name}")
    rows = cursor.fetchall()
    
    # Get column names
    cursor.execute(f"DESCRIBE {table_name}")
    columns = [col[0] for col in cursor.fetchall()]
    
    print(f"\n=== {table_name} Data ===")
    print(f"Found {len(rows)} rows")
    
    # Group by Type
    data_by_type = {}
    for row in rows:
        row_dict = dict(zip(columns, row))
        frag_type = row_dict.get('Type', 'Unknown')
        if frag_type not in data_by_type:
            data_by_type[frag_type] = []
        data_by_type[frag_type].append(row_dict)
    
    # Print data by type
    for frag_type, fragments in data_by_type.items():
        print(f"\n{frag_type}s ({len(fragments)}):")
        for frag in fragments:
            print(f"  - {frag['Name']}:")
            for col, val in frag.items():
                if col != 'Name':
                    print(f"    {col}: {val}")

def main():
    """Main function to explore database tables"""
    conn, cursor = connect_to_database()
    if not conn:
        print("Failed to connect to database")
        return
        
    try:
        # Get all tables
        tables = get_tables(cursor)
        print("\nFound tables:", tables)
        
        # Explore each table
        for table in tables:
            explore_table(cursor, table)
            
    except Exception as e:
        print(f"Error exploring database: {e}")
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    main()
