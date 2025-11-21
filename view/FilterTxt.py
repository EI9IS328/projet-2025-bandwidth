import pandas as pd

# Input and output file paths
input_txt = "res.txt"
output_txt = "data.txt"

# Read the data into a pandas DataFrame, ensuring we handle space separation correctly
df = pd.read_csv(input_txt, delim_whitespace=True)

# Print the first few rows and column names to inspect the data
print("Column names:", df.columns)
print("First few rows of the data:\n", df.head())

# Strip any leading/trailing whitespace from the column names
df.columns = df.columns.str.strip()

# Filter only the columns: Step, X, Y, Z, pnGlobal
try:
    filtered_df = df[['Step', 'X', 'Y', 'Z', 'pnGlobal']]

    # Function to format the values
    def format_value(val):
        if val == 0:
            return "0"
        else:
            return f"{val:.6e}"  # Use scientific notation with 6 significant figures

    # Apply the formatting function to the 'pnGlobal' column
    filtered_df['pnGlobal'] = filtered_df['pnGlobal'].apply(format_value)

    # Save the filtered data back to a new .txt file with space separation
    filtered_df.to_csv(output_txt, index=False, header=True, sep=' ', float_format='%.6e')

    print(f"Filtered data saved to: {output_txt}")
except KeyError as e:
    print(f"Error: {e}")
    print("Check if column names are correct. Available columns are:", df.columns)
