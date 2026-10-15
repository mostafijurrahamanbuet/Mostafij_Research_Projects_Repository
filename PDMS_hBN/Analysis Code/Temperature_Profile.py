def extract_and_save_data_with_average_summary(input_file_path, output_file_path, start_step, end_step, chunk_count):
    """
    Extracts and saves detailed temperature data and a final summary of averages.

    Args:
        input_file_path (str): The path to the data file.
        output_file_path (str): The path to the output file.
        start_step (int): The starting timestep.
        end_step (int): The ending timestep.
        chunk_count (int): The total number of chunks.
    """
    chunk_temps = {chunk: [] for chunk in range(1, chunk_count + 1)}

    with open(input_file_path, 'r') as file:
        current_step = None
        for line in file:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit():
                current_step = int(parts[0])
                if current_step < start_step or current_step > end_step:
                    current_step = None
                continue

            if current_step is not None and len(parts) == 5:
                try:
                    chunk_num = int(parts[0])
                    temperature = float(parts[3])
                    chunk_temps[chunk_num].append(temperature)
                except (ValueError, IndexError):
                    continue

    all_chunk_averages = {}

    with open(output_file_path, 'w') as out_file:
        for chunk in range(1, chunk_count + 1):
            out_file.write(f"chunk : {chunk}\n\n")

            temperatures = chunk_temps.get(chunk, [])
            if temperatures:
                avg_temp = sum(temperatures) / len(temperatures)
                out_file.write(f"Average Temperature: {avg_temp:.3f}\n\n")
                all_chunk_averages[chunk] = avg_temp

                for temp in temperatures:
                    out_file.write(f"\t\t\t\t{temp:.3f}\n")
            out_file.write("-" * 20 + "\n")

        # Write the final summary table at the end of the file
        out_file.write("\n\n" + "=" * 30 + "\n")
        out_file.write("Summary of Average Temperatures\n")
        out_file.write("=" * 30 + "\n")
        out_file.write("Chunk Number\t\tAverage Temp\n")
        out_file.write("-" * 30 + "\n")

        for chunk_num, avg_temp in all_chunk_averages.items():
            out_file.write(f"{chunk_num}\t\t\t{avg_temp:.3f}\n")


# Example usage with your data file and desired output file
input_file_name = '50%_Rectangular_tmp_profile.data'
output_file_name = 'another_try_50%_Rectangular_tmp_profile_x1ns.txt'
start = 52300000
end = 56300000
chunks = 89

extract_and_save_data_with_average_summary(input_file_name, output_file_name, start, end, chunks)
