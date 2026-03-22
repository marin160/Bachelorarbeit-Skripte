import os
import glob
import subprocess
import pandas as pd
from pathlib import Path

acc = []
names = []
homologs = []
evalues = []
fasta = []

def extract_alphafold_seq(accession, filename):
    cmd = f"fetch nr50:{accession} > {filename}.fasta"
    subprocess.run(cmd, shell=True, check=True)



def blast_search(fasta_file, acc, names, homologs, evalues, fasta, db='nr50', evalue='1E-5', num_threads=12, min_cov=0.8, max_slen=1500):

    base_name = Path(fasta_file).stem
    base_name = base_name.replace("fold_", "")
    base_name = base_name.replace("_human_model_0_A", "")
    
    tab1 = f"{base_name}_tab1.txt"
    tab2 = f"{base_name}_tab2.txt"
    
    cmd = [
        'blastp',
        '-query', fasta_file,
        '-db', db,
        '-num_threads', str(num_threads),
        '-evalue', evalue,
        '-outfmt', '6 qseqid sseqid pident qlen slen length evalue',
        '-out', tab1
    ]
    subprocess.run(cmd, check=True)
    
    columns = ['qseqid', 'sseqid', 'pident', 'qlen', 'slen', 'length', 'evalue']
    df = pd.read_csv(tab1, sep='\t', header=None, names=columns)
    
    #Filter: length >= 80% qlen and slen < 1500
    df_filtered = df[(df['length'] >= min_cov * df['qlen']) & (df['slen'] < max_slen)].copy()
    
    df_filtered = df_filtered.sort_values('evalue')

    
    df_filtered.to_csv(tab2, sep='\t', index=False, header=True)
    
    os.remove(tab1)
    
    print(f"Processed {fasta_file}: {len(df_filtered)} hits saved to {tab2}")

    bins= [
        (95,99,'H99'),
        (90, 95,'H95'),
        (85,90,'H90'),
        (80, 85,'H85'),
        (75,80,'H80'),
        (70, 75,'H75'),
        (65,70,'H70'),
        (60, 65,'H65'),
        (55, 60,'H60'),
        (50,55,'H55'),
        (45, 50,'H50'),
        (40, 45,'H45'),
        (35,40, 'H40'),
        (30, 35,'H35'),
        (25, 30,'H30'),
        (20,25,'H25'), ]
    #bins ist als tupel verpackt, tupel = (lower, upper und suffix)
    
    os.makedirs("homolog_fasta", exist_ok=True)

    for lower, upper, suffix in bins:
        pident_hits= df_filtered[(df_filtered['pident']>=lower) & (df_filtered ['pident']< upper)]
    
        if len(pident_hits) >0: 
            best_hit = pident_hits.head(1) 
            accession = best_hit["sseqid"].iloc[0].split("|")[1]
            acc.append(accession)
            names.append(base_name)
            homologs.append(best_hit["pident"].iloc[0])
            evalues.append(best_hit["evalue"].iloc[0])
            fasta.append(f"{base_name}_{suffix}.fasta")
            
            fsta = f"homolog_fasta/{base_name}_{suffix}"
            extract_alphafold_seq(accession, fsta)
            modify_fasta(fsta + ".fasta")

def modify_fasta(input_path):
    with open(input_path, 'r') as f:
        content= f.read()
    
    lines= content.splitlines()
    if not lines or not lines[0].startswith('>'):
        print(f"Skipped {input_path}: Fasta file is broken.")
        return
    
    modified_lines= [lines[0]]+[line.replace('B','D').replace('Z','E').replace('U','C').replace('O','X')for line in lines[1:]]
    modified_content='\n'.join(modified_lines)

    with open(input_path, 'w') as f:
        f.write(modified_content)


def main():
    fasta_files = glob.glob("*A.fasta")
    
    if not fasta_files:
        print("No .fasta files found in the current directory.")
        return
    

    for fasta_file in fasta_files:
        blast_search(fasta_file, acc, names, homologs, evalues, fasta)

    best_hit_table = pd.DataFrame({
        "Complex name:": names,
        "Uniprot Accession Number:": acc,
        "% Identity:": homologs,
        "E-Value:": evalues,
        "Homolog Fasta:": fasta
    })
    best_hit_table.to_csv("best_hit_table.csv", index= False)
    print("The best_hit_table.csv saved.")

if __name__ == "__main__":
    main()

