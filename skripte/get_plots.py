import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Example: python get_plots.py /path/to/summary_csv")
    sys.exit(1)

csv_folder = sys.argv[1]

plot_dir = os.path.join(csv_folder, "plots")
os.makedirs(plot_dir, exist_ok=True)


pattern = os.path.join(csv_folder, "*_summary.csv")
csv_files = glob.glob(pattern)

if not csv_files:
    print("No *_summary.csv file found.")
    sys.exit(0)

for csv_file in csv_files:
    df = pd.read_csv(csv_file)

    
    df = df.dropna(subset=["% Diff:", "ipTM:"])

    df = df.sort_values("% Diff:")

    x = df["% Diff:"]
    y = df["ipTM:"]

    if x.empty or y.empty:
        print(f"Skipped no fli found: {csv_file}")
        continue

    
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, linewidth=1.5, marker='o', markersize=5, color="#8A2BE2")

    
    plt.xlim(0, 100)
    plt.ylim(0, 1)
    plt.xticks(range(0, 101, 10))
    plt.yticks([i/10 for i in range(11)])

    
    plt.xlabel("% sequence divergence", fontsize=12)  
    plt.ylabel("ipTM score", fontsize=12)            
 #The axis labelling can be adjusted here if necessary.   
    
    title = os.path.basename(csv_file).replace("_summary.csv", "")
    plt.title(title)

    
    plt.grid(True, linestyle="--", alpha=0.3)

    
    out_png = os.path.join(plot_dir, f"{title}_iptm_vs_divergence.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"Plots saved: {out_png}")

print("DONE: All Plots and weighted means saved.")

