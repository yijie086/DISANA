import os
from collections import defaultdict

def calculate_total_charge(filepaths):
    """
    Reads one or more space-separated text files, skips header lines, and calculates
    the total sum of 'gated charge' and 'ungated charge' from the columns.

    Args:
        filepaths (list): A list of file paths (strings) to process.

    Returns:
        tuple: A tuple (total_gated_charge, total_ungated_charge) containing the sums.
    """
    total_charges = defaultdict(float)

    # Define the indices for the columns we want to sum (0-indexed)
    # Based on the file structure:
    # 0: run number
    # 1: gated charge
    # 2: ungated charge
    GATED_CHARGE_INDEX = 1
    UNGATED_CHARGE_INDEX = 2

    print(f"--- Processing {len(filepaths)} Files ---")

    for filepath in filepaths:
        # Check if the file exists before trying to open it
        if not os.path.exists(filepath):
            print(f"Error: File not found at {filepath}. Skipping.")
            continue

        file_gated_sum = 0.0
        file_ungated_sum = 0.0
        line_count = 0
        
        try:
            with open(filepath, 'r') as f:
                # Iterate through lines
                for line in f:
                    # Skip the first two header/separator lines
                    if line_count < 2:
                        line_count += 1
                        continue

                    # Split the line by whitespace
                    parts = line.split()

                    # We expect exactly 3 columns of data: run_number, gated_charge, ungated_charge
                    if len(parts) >= 3:
                        try:
                            # Convert the charge values to float and sum them up
                            gated_charge = float(parts[GATED_CHARGE_INDEX])
                            ungated_charge = float(parts[UNGATED_CHARGE_INDEX])
                            
                            file_gated_sum += gated_charge
                            file_ungated_sum += ungated_charge
                            
                        except ValueError:
                            # Catch lines where conversion to float fails (e.g., malformed data)
                            print(f"Warning: Skipped line in {filepath} due to non-numeric data: '{line.strip()}'")
                            
                    line_count += 1
            
            # Add the file's sums to the overall totals
            total_charges['gated'] += file_gated_sum
            total_charges['ungated'] += file_ungated_sum
            
            print(f"Processed '{filepath}':")
            print(f"  Gated Charge Sum:   {file_gated_sum:,.3f}")
            print(f"  Ungated Charge Sum: {file_ungated_sum:,.3f}")
            
        except Exception as e:
            print(f"An unexpected error occurred while reading {filepath}: {e}")
            
    # Return the final totals
    return total_charges['gated'], total_charges['ungated']


if __name__ == "__main__":
    # List the file names to be processed. These must be in the same directory as the script.
    files_to_analyze = [
        "charge_sp18_outb.txt",
        "charge_sp19.txt"
    ]

    gated_total, ungated_total = calculate_total_charge(files_to_analyze)

    print("\n" + "="*40)
    print("           OVERALL TOTALS")
    print("="*40)
    print(f"Total Gated Charge (All Runs):   {gated_total:,.3f}")
    print(f"Total Ungated Charge (All Runs): {ungated_total:,.3f}")
    print("="*40)