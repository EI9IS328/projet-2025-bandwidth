import csv
import matplotlib.pyplot as plt

def size_of(row):
    return int(row["ex"]) * int(row["ey"]) * int(row["ez"])

x, overhead = [], []
with open("semproxy_times.csv", newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        total = float(row["elapsed_total_time_s"])
        compute = float(row["elapsed_compute_time_s"])
        x.append(size_of(row))
        overhead.append(total - compute)

x, overhead = zip(*sorted(zip(x, overhead)))

plt.plot(x, overhead, marker="o")
plt.xlabel("Problem size (ex × ey × ez)")
plt.ylabel("Overhead time (s) = total − compute")
plt.title("SEMPROXY — Overhead vs problem size")
plt.grid(True)
plt.tight_layout()
plt.savefig("overhead_vs_size.png", dpi=200)
plt.show()