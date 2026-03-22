import sys
import json
import os
import pandas as pd
import glob
from pathlib import Path


if len(sys.argv) !=5:
    print("ERROR:NUMBER OF ARGUMENTS SHOULD BE  EXACTLY 5")
    print("python get_iptm.py path/to/AF_models path/to/best_hit_tab.csv path/to/iptm_start_scores.csv output.csv  ")
    sys.exit(1)



maindir= sys.argv[1]
best_hit_table_csv= sys.argv[2]
itpm_score_start=sys.argv[3]
output_csv= sys.argv[4]



iptm_scores=[]
homolog_names=[]



for subdir in os.listdir(maindir):
    subdir_path= os.path.join(maindir, subdir)


    if not os.path.isdir(subdir_path):
        continue


    json_string = os.path.join(subdir_path, "*_summary_confidences.json")
    json_paths = glob.glob(json_string)


    if not json_paths:
        continue
   
    json_file = json_paths[0]


    with open(json_file) as f:
        data = json.load(f)
        iptm_score = data['iptm']
   
    homolog_name = os.path.basename(json_file)
    homolog_name = homolog_name.replace("_summary_confidences.json", "")

    homolog_name_short = homolog_name.replace("fold_", "").replace("_human_model_0_a", "").replace("_model_a", "")
    homolog_names.append(homolog_name_short)
    iptm_scores.append(iptm_score)
    

df = pd.DataFrame({"Name:": homolog_names, "ipTM:": iptm_scores})
best_hit_df = pd.read_csv(best_hit_table_csv)

def norm_name(s):
    if pd.isna(s):
        return ""
    s = str(s).strip().lower()
    for ext in [".fasta", ".fa", ".faa", ".fna", ".csv", ".txt"]:
        if s.endswith(ext):
            s = s[:-len(ext)]
    s = s.replace("fold_", "").replace("_summary_confidences", "")
    s = s.replace("_human_model_0_a", "").replace("_model_a", "")
    s = s.replace(" ", "").replace("|", "_")
    return s


df["name_lower"] = df["Name:"].apply(norm_name)
best_hit_df["Homolog Fasta:_lower"] = best_hit_df["Homolog Fasta:"].apply(norm_name)


merged = df.merge(
    best_hit_df[[
        "Complex name:", 
        "Homolog Fasta:_lower", 
        "% Identity:", 
        "Uniprot Accession Number:"
    ]],
    left_on="name_lower",
    right_on="Homolog Fasta:_lower",
    how="left"
)
print("Unmatched:", merged["Complex name:"].isna().sum(), "of", len(merged))

merged = merged.drop(columns=["name_lower", "Homolog Fasta:_lower"])

merged["% Identity:"] = merged["% Identity:"].fillna(0)
merged["% Diff:"] = 100 - merged["% Identity:"]

column_order = ["Name:", "Complex name:", "Uniprot Accession Number:", "% Identity:", "% Diff:", "ipTM:"]
merged = merged.reindex(columns=column_order).sort_values("Name:")

start_iptm_df=pd.read_csv(itpm_score_start)
start_iptm_df["Uniprot Accession Number:"]= None
start_iptm_df ["% Identity:"]=100
start_iptm_df["% Diff:"]=0.0
start_iptm_df = start_iptm_df[column_order]

merged_with_start_iptm= pd.concat([merged, start_iptm_df], ignore_index=True)
merged_with_start_iptm= merged_with_start_iptm.sort_values("Name:")

output_dir = Path(output_csv).parent
output_dir.mkdir(parents=True, exist_ok=True)

merged_with_start_iptm["Complex name:"] = merged_with_start_iptm["Complex name:"].str.replace("_model_A", "")

for query_name, sub_df in merged_with_start_iptm.groupby("Complex name:"):
    out_file = output_dir / f"{query_name}_summary.csv"
    sub_df.to_csv(out_file, index=False)
print("DONE: All Query-Tables are saved")
