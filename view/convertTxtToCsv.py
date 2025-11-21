import pandas as pd

# Input and output file paths
input_txt = "data.txt"
output_csv = "data.csv"

# Read the space-separated .txt file into a pandas DataFrame
df = pd.read_csv(input_txt, delim_whitespace=True)

# Save the DataFrame to a .csv file
df.to_csv(output_csv, index=False)

print(f"CSV file saved to: {output_csv}")
