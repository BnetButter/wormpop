import sqlite3
import sys
import os
import csv

def execute_sql(database, sqlscript):
    # Check if the files exist
    if not os.path.isfile(database):
        print(f"Database file '{database}' does not exist.")
        return
    if not os.path.isfile(sqlscript):
        print(f"SQL script file '{sqlscript}' does not exist.")
        return

    # Read the SQL script
    with open(sqlscript, 'r') as file:
        sql_query = file.read()

    # Connect to the SQLite database
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    try:
        # Execute the SQL script
        cursor.execute(sql_query)
        results = cursor.fetchall()

        # Get column names
        column_names = [description[0] for description in cursor.description]

        # Use csv.writer to print to stdout
        csv_writer = csv.writer(sys.stdout)
        
        # Write the headers

        # Write the rows vertically
        for row in results:
            for col_name, value in zip(column_names, row):
                csv_writer.writerow([col_name, value])
            csv_writer.writerow([])  # Blank line between rows for better readability

    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
    
    finally:
        # Close the connection
        conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <database> <sqlscript>")
    else:
        database = sys.argv[1]
        sqlscript = sys.argv[2]
        execute_sql(database, sqlscript)
