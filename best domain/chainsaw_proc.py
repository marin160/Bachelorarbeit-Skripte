#Das soll den Chainsaw output lesen und daraus die erste und vierte Spalte extrahieren und das als eine neue tsv Datei speichern, hier werden die Bereiche der Domänen bestimmt
import os
import sys
import pandas as pd

target_path= "temp/domains.tsv"


def get_tsv_data(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t")
    df_column = df.iloc[:,[0,4]]
    print(df_column)
    return df_column

def arguments(): 
    if len(sys.argv) !=2: 
        print("ERROR:NUMBER OF ARGUMENTS SHOULD BE  EXACTLY 2")
   
    tsvpath = sys.argv[1]

    df_out= get_tsv_data(tsvpath)
    df_out.columns = ["name", "domain_cutoffs"]
    

    df_orig = pd.read_csv(tsvpath, sep="\t") 
    nan_mask = df_out["domain_cutoffs"].isna()
    for name in df_out[nan_mask]["name"]:
        nres = df_orig[df_orig["chain_id"] == name]["nres"].iloc[0]
        print(f"{name}: nres={nres} → 1-{nres}")
        df_out.loc[df_out["name"] == name, "domain_cutoffs"] = f"1-{nres}"

    
    print(df_out)

    os.makedirs("temp", exist_ok= True)
    df_out.to_csv(target_path, sep="\t", index=False)
    print(f"domains.tsv SAVED IN {target_path}")

if __name__ == "__main__":
    arguments()