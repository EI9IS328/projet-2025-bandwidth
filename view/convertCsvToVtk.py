import pandas as pd

csv_file = "data.csv"
vtk_file = "data.vtk"

# Load CSV
df = pd.read_csv(csv_file)

# Expect columns: x, y, z, pnglobal
x = df["x"]
y = df["y"]
z = df["z"]
p = df["pnglobal"]

with open(vtk_file, "w") as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("Point data exported from CSV\n")
    f.write("ASCII\n")
    f.write("DATASET POLYDATA\n")
    f.write(f"POINTS {len(df)} float\n")

    # write coordinates
    for xi, yi, zi in zip(x, y, z):
        f.write(f"{xi} {yi} {zi}\n")

    # write per-point scalar
    f.write(f"\nPOINT_DATA {len(df)}\n")
    f.write("SCALARS pnglobal float 1\n")
    f.write("LOOKUP_TABLE default\n")

    for pi in p:
        f.write(f"{pi}\n")

print("Saved:", vtk_file)
