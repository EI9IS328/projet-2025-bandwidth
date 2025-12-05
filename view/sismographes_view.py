import numpy as np
import matplotlib.pyplot as plt

# Charger le fichier
data = np.loadtxt("datas.txt", skiprows=1)

steps = data[:, 0]
pnglobal = data[:, 4]
X = data[:, 1]
Y = data[:, 2]
Z = data[:, 3]

plt.figure(figsize=(8,5))
plt.plot(steps, pnglobal, marker='o', linestyle='-', color='b')

# Ajouter les valeurs comme étiquettes à côté des points
for x, y in zip(steps, pnglobal):
    plt.text(x, y, "", fontsize=8, ha='left', va='bottom')

plt.xlabel("Step")
plt.ylabel("pnglobal")
plt.title("Evolution de pnglobal en fonction des Steps sur" + f"(X={X[1]:.0f}, Y={Y[1]:.0f}, Z={Z[1]:.0f})")
plt.grid(True)
plt.tight_layout()
plt.show()
