import csv
import matplotlib.pyplot as plt
import seaborn

with open("top1percent.csv", "r") as f:
    reader = csv.DictReader(f)
    centers = []
    pmins = []
    pmaxes = []
    deltap = []
    slopes = []
    for row in reader:
        centers.append(float(row["center"]))
        pmins.append(float(row["Pmin"]))
        pmaxes.append(float(row["Pmax"]))
        deltap.append(float(row["Pmax-Pmin"]))
        slopes.append(float(row["slope"]))

plt.hist(centers, 50)
plt.xlabel("ML Estimate of Cline Center (km)")
plt.ylabel("Count")
plt.show()

# plt.figure()
plt.hist(pmins, 50)
plt.xlabel("ML estimate of P min")
plt.ylabel("Count")
plt.show()

# plt.figure()
plt.hist(pmaxes, 50)
plt.xlabel("ML estimate of P max")
plt.ylabel("Count")
plt.show()

# plt.figure()
plt.hist(deltap, 50)
plt.xlabel("ML Estimate of $\Delta$P")
plt.ylabel("Count")
plt.show()

# plt.figure()
plt.hist(slopes, 50, range=(0, .35))
plt.xlabel("ML Estimate of Cline Slope")
plt.ylabel("Count")
plt.show()
