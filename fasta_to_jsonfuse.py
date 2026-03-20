import os
import glob
import subprocess

os.makedirs("homolog_fasta/jsons", exist_ok=True)
ubiquitin = "/database/AF3_input/helpers/ubraw.json"

for fasta_file in glob.glob("homolog_fasta/*.fasta"):
     print("FILE:", fasta_file)
     name = os.path.splitext(os.path.basename(fasta_file))[0]
     print("NAME:", name)


     subprocess.run(["fasta2json",fasta_file], check=True)


     json_file = f"homolog_fasta/{name}.json"


     cmd=[
        "jsonfuse",
        "--json1", json_file,
        "--id1", "A",
        "--json2", ubiquitin,
        "--id2", "A",
        "--seed", "1",
        "--name", name,
        "--output", f"homolog_fasta/jsons/{name}.json"
        ]
     
     print("CMD:", cmd)
     subprocess.run(cmd, check= True)
     print(f"saved{name}.json in /jsons")