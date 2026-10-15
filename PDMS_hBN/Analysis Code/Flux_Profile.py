import pandas as pd
import os


def process_flux_file(input_filename, output_filename, min_timestep=52300000, max_timestep=56300000):
    """
    Reads a LAMMPS output file, filters the data by timestep, and saves it to an Excel file.

    Args:
        input_filename (str): The name of the input file (e.g., 'flux_1700000.out').
        output_filename (str): The name of the output Excel file (e.g., 'Flux.xlsx').
        min_timestep (int): The minimum timestep to include in the output.
        max_timestep (int): The maximum timestep to include in the output.
    """
    if not os.path.exists(input_filename):
        print(f"Error: The input file '{input_filename}' was not found.")
        return

    try:
        # Find the number of comment lines to skip
        skiprows_count = 0
        with open(input_filename, 'r') as f:
            for line in f:
                if line.strip().startswith('#'):
                    skiprows_count += 1
                else:
                    break

        # Read the data into a pandas DataFrame.
        # We specify 'header=None' because the data file does not have a header row.
        # The 'names' parameter is used to manually assign column names.
        df = pd.read_csv(
            input_filename,
            sep='\s+',
            skiprows=skiprows_count,
            header=None,  # No header row in the data
            names=['TimeStep', 'f_hot_rescale', 'f_cold_rescale1', 'f_cold_rescale2']
        )

        # Filter the DataFrame to keep only rows within the desired timestep range
        filtered_df = df[(df['TimeStep'] >= min_timestep) & (df['TimeStep'] <= max_timestep)]

        if filtered_df.empty:
            print(f"No data found in the specified range [{min_timestep} - {max_timestep}].")
            return

        # Write the filtered DataFrame to an Excel file
        filtered_df.to_excel(output_filename, index=False, engine='openpyxl')

        print(f"Successfully filtered and converted data from '{input_filename}' to '{output_filename}'")
        print(f"Only data in the timestep range [{min_timestep} - {max_timestep}] was included.")
        print("The first few rows of the filtered data are:")
        print(filtered_df.head())

    except Exception as e:
        print(f"An error occurred: {e}")


# The name of your input file
input_file = 'Rectangular_50%_flux.out'

# The name of the output Excel file
output_file = 'new_try_Rectangular_50%_flux1ns.xlsx'

# Call the function with the desired timestep range
# You can change these values to filter for a different range
process_flux_file(input_file, output_file)