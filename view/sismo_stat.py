import glob
import csv
import statistics
import os

def analyze_sismo_files(directory="../build"):
    # Cherche tous les fichiers du type sismo*.csv
    pattern = os.path.join(directory, "recev_results_*.csv")
    files = sorted(glob.glob(pattern))

    if not files:
        print("Aucun fichier sismo*.csv trouvé.")
        return

    print(f"{len(files)} fichiers trouvés.\n")


    for filename in files:
        values = []

        with open(filename, "r", newline="") as f:
            reader = csv.DictReader(f)

            for row in reader:
                try:
                    values.append(float(row["varnp1"]))
                except (KeyError, ValueError):
                    continue

        if not values:
            print(f"{filename} : aucune donnée valide")
            continue

        vmin = min(values)
        vmax = max(values)
        vmean = statistics.mean(values)
       # vvar = statistics.variance(values)


        print(f"Fichier : {filename}")
        print(f"  min     = {vmin:e}")
        print(f"  max     = {vmax:e}")
        print(f"  moyenne = {vmean:e}")
        #print(f"  variance = {vvar:e}")

        print("-" * 40)

if __name__ == "__main__":
    analyze_sismo_files()
