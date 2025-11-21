import pandas as pd

# Input CSV file and output VTK file paths
input_csv = "data.csv"
output_vtk = "data.vtk"

# Read the CSV into a pandas DataFrame
df = pd.read_csv(input_csv)

# Extract the columns of interest: Step, X, Y, Z, pnGlobal
x = df['X']
y = df['Y']
z = df['Z']
pnGlobal = df['pnGlobal']

# Write the VTK file in PolyData format
with open(output_vtk, "w") as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("Converted CSV to VTK\n")
    f.write("ASCII\n")
    f.write("DATASET POLYDATA\n")
    f.write(f"POINTS {len(df)} float\n")

    # Write the points (X, Y, Z)
    for xi, yi, zi in zip(x, y, z):
        f.write(f"{xi} {yi} {zi}\n")

    # Write the point data (pnGlobal) as scalars
    f.write(f"\nPOINT_DATA {len(df)}\n")
    f.write("SCALARS pnGlobal float 1\n")
    f.write("LOOKUP_TABLE default\n")

    for pi in pnGlobal:
        f.write(f"{pi}\n")

print(f"VTK file saved to: {output_vtk}")
