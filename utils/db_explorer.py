#!/usr/bin/env python3
import sys
import os
import pymysql
from dotenv import load_dotenv

def explore_database():
    """Explore the database structure"""
    try:
        # Load environment variables
        load_dotenv()
        
        # Database connection parameters
        db_params = {
            'host': os.getenv('DB_HOST'),
            'user': os.getenv('DB_USER'),
            'password': os.getenv('DB_PASSWORD'),
            'db': os.getenv('DB_NAME'),
            'charset': 'utf8mb4',
            'cursorclass': pymysql.cursors.DictCursor
        }
        
        # Connect to database
        connection = pymysql.connect(**db_params)
        cursor = connection.cursor()
        
        try:
            # Get all tables
            cursor.execute("SHOW TABLES")
            tables = cursor.fetchall()
            print("\nTables in database:")
            for table in tables:
                table_name = list(table.values())[0]  # Get the table name from the dict
                print(f"\n=== Table: {table_name} ===")
                
                # Get columns for each table
                cursor.execute(f"DESCRIBE {table_name}")
                columns = cursor.fetchall()
                print("Columns:")
                for col in columns:
                    print(f"  {col['Field']}: {col['Type']}")
                    
                # Get sample data
                cursor.execute(f"SELECT * FROM {table_name} LIMIT 1")
                sample = cursor.fetchone()
                if sample:
                    print("\nSample row:")
                    for key, value in sample.items():
                        print(f"  {key}: {value}")
                        
        finally:
            cursor.close()
            connection.close()
            
    except Exception as e:
        print(f"Error exploring database: {str(e)}")

if __name__ == "__main__":
    explore_database()
