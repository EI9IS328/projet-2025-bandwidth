import csv
import matplotlib.pyplot as plt

def problem_size(ex, ey, ez):
    return ex * ey * ez

sizes = []
compute_times = []

with open("semproxy_times.csv", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        ex = int(row["ex"])
        ey = int(row["ey"])
        ez = int(row["ez"])
        sizes.append(problem_size(ex, ey, ez))
        compute_times.append(float(row["elapsed_compute_time_s"]))

#sort by size
sizes, compute_times = zip(*sorted(zip(sizes, compute_times)))

plt.plot(sizes, compute_times, marker="o")
plt.xlabel("Problem size (ex × ey × ez)")
plt.ylabel("Compute time (s)")
plt.title("SEMPROXY – Compute time vs problem size")
plt.grid(True)
plt.savefig("compute_time_vs_size.png", dpi=200)
plt.show()