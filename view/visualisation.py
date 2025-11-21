import os
import pandas as pd
import matplotlib.pyplot as plt


def plot_from_csv():
    csv_path = "data.csv"
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Fichier introuvable: {csv_path}")

    df = pd.read_csv(csv_path)

    expected_cols = {"Step", "pnglobal"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"Colonnes manquantes dans le CSV : {missing}")

    # trié le step
    df = df.sort_values("Step").reset_index(drop=True)

    # Garder que les colonnes Step et pnglobal et supprimer les lignes nulles
    df = df.dropna(subset=["Step", "pnglobal"])

    # Création du graphe
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(df["Step"], df["pnglobal"], marker='x', linestyle='--', label='pnglobal')
    ax.set_xlabel("Step")
    ax.set_ylabel("pnglobal") 
    ax.grid(True)
    ax.legend(loc='upper right')

    plt.title("Analyse temporelle du pnglobal")
    plt.tight_layout()

    output_path = "plot.png"
    fig.savefig(output_path, dpi=150)

    plt.show()
    plt.close(fig)


plot_from_csv()
