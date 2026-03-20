import sys
import os
import shutil
import zipfile
from pymol import cmd, stored

def arguments():
    if len(sys.argv) !=3:
        print("ERROR:NUMBER OF ARGUMENTS SHOULD BE  EXACTLY 3")
        
    input_dir= sys.argv[1]
    program_id= sys.argv[2]

    if program_id not in ("pymol","plip", "afb"):
        print("ERROR: program_id must be pymol, plip or afb") 
        
    if os.path.isdir(input_dir): 
        print("DIRECTORY EXISTS")
    else:
        print ("ERROR: DIRECTORY DOES NOT EXIST")
        
        
    return input_dir, program_id



#Hier werden Ordner erstellt die für die nachfolgenden arbeiten gebraucht werden 

def make_directories(program_id):
    os.makedirs("temp/2chains", exist_ok=True)
    
    if program_id == "plip": 
        os.makedirs("temp/plip_in", exist_ok=True)
    if program_id == "afb": 
         os.makedirs("temp", exist_ok=True)




 
#Diese Funktion ist dafür da um die Komplexe in ihre Ketten zu zerteilen. Das findet in einer PyMol umgebung statt
def split_chains(input_path, output_dir):
    if not os.path.isfile(input_path):
        print(f"NOT FOUND: {input_path}")
        return
    if not input_path.endswith((".pdb", ".cif")):
        return
    
    cmd.reinitialize() 
    cmd.load(input_path) 
    #Erstellt eine leere PyMol Session und Lädt den Komplex in PyMol

    complex_name = os.path.splitext(os.path.basename(input_path))[0]
    #Extrahiert den Dateinamen (Bspl: KleinerWauWau.pdb -> "KleinerWauWau")

    chains = cmd.get_chains(complex_name)
    # holt alle Ketten_IDs (A, B, ...) aus der complex_name PyMol Datei
     
    if not chains:
        print(f"ERROR: DOMAIN NOT FOUND IN {input_path}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    #Falls die Datei mit dem Namen bereits existiert wird mit der nächsten weitergemacht ohne die alte zu überschreiben
     
    for chain in chains:
        if chain != "":
            chain_id = chain
        else:
            chain_id = "NO_CHAIN"

        chain_complex_name = f"{complex_name}{chain_id}"
        
        cmd.create(chain_complex_name, f"{complex_name} and chain {chain}")
        #Nimmt alle Aminosäuren die im Komplex UND AUCH in der Kette vorkommen und speichert sie mit dem Namen chain_complex_name ab

        out_file_pdb = os.path.join(output_dir, f"{chain_complex_name}.pdb")
        cmd.save(out_file_pdb, chain_complex_name)
        print(f"Saved: {out_file_pdb}")
    cmd.delete("all")

#Diese Funktion zeigt was passieren soll wenn pymol ausgewählt wird
def id_pymol(input_dir): 
    output_dir = os.path.join("temp", "2chains")
    for files in os.listdir(input_dir):
        input_path= os.path.join(input_dir, files)
        if not os.path.isfile(input_path):
            continue
            #alles wird übersprungen was kein datei ist
            
        #ausführung des codes:
        split_chains(input_path, output_dir)

          
def id_plip(input_dir):
    output_dir = os.path.join("temp", "2chains")
    for files in os.listdir(input_dir):
        input_path= os.path.join(input_dir, files)
        if not os.path.isfile(input_path):
            continue
            #die getrennten ketten müssen auch in plip_in gepsichert werden
            
        #ausführung des codes:
        split_chains(input_path, output_dir)
    prep_plip(input_dir)

#pymol enviroment
def prep_plip(input_dir):
    
    for files in os.listdir(input_dir):
        input_path= os.path.join(input_dir, files)
        if not os.path.isfile(input_path):
            continue

        base, ext = os.path.splitext(files)
        ext= ext.lower()
        #das soll den Dateinamen einmal vor den punkt trennen und alles davor nimmt und einmal mit punkt alles danach
        #beispielsweise fold_ptgisxubb_human_model_0.pdb, base wäer hier: fold_ptgisxubb_human_model_0  ext: .pdb
        # ext.lower() ist nur eine vorsichtmaß namen sodass das .pdb oder welcher datei auch immer auf jeden fall klein geschrieben ist, wichtig für die nächsten schritte

        if ext== ".pdb":
            plip_in= "temp/plip_in"
            target_dir = os.path.join(plip_in, files)
            shutil.copy2(input_path, target_dir)
            print("PDB FILES COPIED TO", target_dir)

        #pymol aufrufen und cif datei in pdb umwandeln
        elif ext == ".cif": 
            plip_in= "temp/plip_in"
            cmd.reinitialize() 
            cmd.load(input_path)
            complex_name = f"{base}.pdb" 
            out_file_pdb = os.path.join(plip_in, complex_name)
            cmd.save(out_file_pdb, base)
            print(f"SUCCESFULLY CONVERTED TO .pdb")
            cmd.delete("all")

def afb_model_0_extractor(input_dir):
    temp_dir = "temp"
    os.makedirs(temp_dir, exist_ok=True)
    #input_dir ist das Hauptverzeichnis 
    #input_path ist der vollständige Pfad zu einem einzlenen eintrag innerhalb des input_dir Verzeichnisses.

    for sub_dir in os.listdir(input_dir):
        input_path = os.path.join(input_dir, sub_dir)
    
        if not sub_dir.endswith(".zip"):
            continue
        
        found = False

        with zipfile.ZipFile(input_path, 'r') as z:
            files = z.namelist()
            for file in files:
                if "model_0" in file and file.endswith(".cif"):
                    file_content = z.read(file)
                    full_path = os.path.join(temp_dir, os.path.basename(file))
                    with open(full_path, 'wb') as f:
                        f.write(file_content)
                    found = True
                    break #wichtig Sobald wir EINE haben, hören wir bei DIESER Zip auf.
        if not found:
            print(f"\n NO model_0 FOUND in {input_path}")

    print(f"\nDONE WITH EXTRACTION OF model_0 INTO {temp_dir}.")

def id_afb(input_dir):
    temp_dir = "temp"
    output_dir = os.path.join("temp", "2chains")
    afb_model_0_extractor(input_dir)
    for file in os.listdir(temp_dir):
        print(file)
        input_path= os.path.join(temp_dir, file)
        print(input_path)
        if not os.path.isfile(input_path): 
            continue 
        split_chains(input_path, output_dir)
 
   
if __name__ == "__main__":
    input_dir, program_id = arguments()
    print("Input_dir:", input_dir)
    print("Program_ID:", program_id)

    make_directories(program_id)

    if program_id == "pymol": 
        id_pymol(input_dir) 
    elif program_id == "plip": 
        id_plip(input_dir) 
    elif program_id == "afb":
        id_afb(input_dir)


        
    