import sys
import json
import os
import pandas as pd
import glob

if len(sys.argv)!= 3:
    print("ERROR: NUMBER OF ARGUMENTS SHOULD BE EXACTLY 3")
    print("python iptm_start_models path/to/models iptm_start.csv")
    sys.exit(1)

inputdir= sys.argv[1]
output_csv = sys.argv[2]

table_rows=[]

for subdir in os.listdir(inputdir):
    subdir_path= os.path.join(inputdir, subdir)


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
   
    basename = os.path.basename(json_file)
    name = basename.replace("_summary_confidences.json", "")

    complex_name= subdir

    table_rows.append({
        "Name:": name, 
        "Complex name:": complex_name,
        "ipTM:": iptm_score })

if not table_rows:
    print("No summary_confidences.json found !")
    sys.exit(1)


df = pd.DataFrame(table_rows)
df.to_csv(output_csv, index=False)
print("DONE: All Query-Tables are saved")

