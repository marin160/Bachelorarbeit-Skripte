import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print("Example: python get_plots_and_means_with_txt.py /path/to/summary_csvs path/to/txt with gen name")
    sys.exit(1)

csv_folder = sys.argv[1]
txt_path =sys.argv[2]

plot_dir = os.path.join(csv_folder, "plots")
os.makedirs(plot_dir, exist_ok=True)


pattern = os.path.join(csv_folder, "*_summary.csv")
csv_files = glob.glob(pattern)

means = pd.DataFrame({
    "Name:": pd.Series(dtype="str"),
    "Gene name:": pd.Series(dtype="str"),
    "Weighted mean:": pd.Series(dtype="float"),
    "Weighted standarddeviation:": pd.Series(dtype="float"),
    "Mean:": pd.Series(dtype="float"),
    "Pseudo mean:": pd.Series(dtype="float"),
    "ipTM > 0.7:":pd.Series(dtype="int"), 
    "Number of Homologs > 80% of initial value:": pd.Series(dtype="int"), 
    "Mean of the highest 4:": pd.Series(dtype="float")
})

if not csv_files:
    print("No *_summary.csv file found.")
    sys.exit(0)

df_txt = pd.read_csv(txt_path, sep='\t')
num_to_gene = dict(zip(df_txt.iloc[:,0].astype(str).str.lower(),
                             df_txt.iloc[:,1].astype(str)))

for csv_file in csv_files:
    df = pd.read_csv(csv_file)

    
    df = df.dropna(subset=["% Diff:", "ipTM:"])

    df = df.sort_values("% Diff:")

    x = df["% Diff:"]
    y = df["ipTM:"]
    z = df["% Identity:"]

    weighted_mean = (y*z).sum() / z.sum()

    weighted_standarddeviation = ((z * (y - weighted_mean)**2).sum() / z.sum())**0.5

    name = df.iloc[0, 0]

    unweighted_mean= y.sum()/len(y)
    #unweighted average of all data

    mean_plus_one= (y.sum()+0)/(len(y)+1)
    #unweighted average of all data + adding one ipTM=0 value

    better07= (y >0.7).sum()
    #number of homologs with iptm scores better than 0.7

    initial_iptm= y.loc[z[z == 100].index[0]]
    better80per= (y > initial_iptm*0.8).sum()
    #number of homologs with scores better than 80% of initial value
    
    top_4_iptm= y.nlargest(4)
    top4_mean = (top_4_iptm.sum()+0.22*(4-len(top_4_iptm)))/4
    #unweighted average of the 4 best homologs, if less than 4, replace the missing scores by 0.22
    
    title = os.path.basename(csv_file).replace("_summary.csv", "")

    num = title.split("_")[0].lower()
    gene_name = num_to_gene.get(num, None)
    
    new_row = pd.DataFrame([{"Name:": name, 
                             "Gene name:":gene_name,
                             "Weighted mean:": weighted_mean, 
                             "Weighted standarddeviation:": weighted_standarddeviation, 
                             "Mean:": unweighted_mean, 
                             "Pseudo mean:": mean_plus_one, 
                             "ipTM > 0.7:": better07, 
                             "Number of Homologs > 80% of initial value:":better80per, 
                             "Mean of the highest 4:": top4_mean}])
    means = pd.concat([means, new_row], ignore_index=True)

    

    if x.empty or y.empty:
        print(f"Skipped {csv_file}, no ipTM / % Identity Values found")
        continue

    
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, linewidth=1.5, marker='o', markersize=5, color="#8A2BE2", label=name)
    plt.axhline(y=weighted_mean, color="green", linestyle="dotted", label="Weighted mean")
    
    #Sigma bands (dont know if useful or not)
    #plt.axhspan(weighted_mean - 2*weighted_standarddeviation, weighted_mean + 2*weighted_standarddeviation, color='royalblue', alpha=0.1, label='2 Sigma band')
    #plt.axhspan(weighted_mean - weighted_standarddeviation, weighted_mean + weighted_standarddeviation, color='royalblue', alpha=0.15, label='1 Sigma band')
    
    plt.legend(loc="lower right")

    plt.xlim(0, 100)
    plt.ylim(0, 1)
    plt.xticks(range(0, 101, 10))
    plt.yticks([i/10 for i in range(11)])

    
    plt.xlabel("% sequence divergence", fontsize=12)  
    plt.ylabel("ipTM score", fontsize=12)            
    #The axis labelling can be adjusted here if necessary.   
    
    

    if gene_name is not None:
        plot_title = gene_name
    else:
        plot_title = title

    plt.title(plot_title)

    
    plt.grid(True, linestyle="--", alpha=0.3)

    
    out_png = os.path.join(plot_dir, f"{plot_title}_iptm_vs_divergence.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"Plots saved: {out_png}")


means.to_csv("weighted_means.csv", index=False)
print("DONE: All Plots and weighted means saved.")

