import csv
import math

filename = "../build/results.csv"   # nom de ton fichier

values = []

with open(filename, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        try:
            values.append(float(row["pnGlobal"]))
        except ValueError:
            pass

if not values:
    print("Aucune valeur trouvée.")
else:
    min_val = min(values)
    max_val = max(values)
    mean_val = sum(values) / len(values)

    print("Statistiques pnGlobal :")
    print(f"  Min  : {min_val:.6e}")
    print(f"  Max  : {max_val:.6e}")
    print(f"  Mean : {mean_val:.6e}")
