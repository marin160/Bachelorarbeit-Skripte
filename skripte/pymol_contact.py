# das skript soll mithilfe von pymol die Kontaktstellen zwischen den beiden Proteinen finden  
import sys
import os
import pandas as pd
from pymol import cmd, stored

#Diese Funktion ist dazu da um aus den Ketten die Kontakt Reste zu ermitteln
def find_contact(filename):
    cmd.reinitialize()
    cmd.load(filename, "complex")
    cmd.split_chains("complex")
    #Erstellt eine leere PyMol Session und Lädt den Komplex in PyMol. Außerdem werden die Komplexe in ihre Ketten aufgeteilt
    
    cmd.select("contactsA", "(chain A) within 5 of (chain B)")
    stored.listA = []
    #Kontakt-Reste von A zu B werden ausgewählt und in einer Liste gespeichert
    
    cmd.iterate("contactsA and name CA", "stored.listA.append(int(resi))")
    print(f"{filename} - A->B: {sorted(set(stored.listA))}")
    #Man erhält eine Liste aller Residues-Nummern der Kontaktstellen. Die Zahlen werden durch den sorted Kommand bei print alle von klein nach groß geordnet.
    #cmd.iterate(selection, expression) geht alle Atome in der Auswahl durch
    
   
    cmd.select("contactsB", "(chain B) within 5 of (chain A)")
    stored.listB = []
    #Kontakt-Reste von B zu A werden ausgewählt und in einer Liste gespeichert
    
    cmd.iterate("contactsB and name CA", "stored.listB.append(int(resi))")
    print(f"{filename} - B->A: {sorted(set(stored.listB))}")
    return stored.listA, stored.listB

def arguments(): 
    if len(sys.argv) !=2: 
        print("ERROR:NUMBER OF ARGUMENTS SHOULD BE  EXACTLY 2")

    input_dir = sys.argv[1]
    
    protein_names=[]
    contact_list = []

    for file in os.listdir(input_dir):
        if not file.endswith((".pdb", ".cif")):
            continue
    

        full_path= os.path.join(input_dir, file)

        contactsA, contactsB= find_contact(full_path)

        filename_base = os.path.splitext(file)[0]

        protein_names.append(filename_base + "_A")
        contact_list.append(",".join(str(x) for x in contactsA))
       # (str(x) for x in contactA) wandelt alle Elemente die in contactA enhlaten sind in strings um
        #","join(str(x) for x in contactA) druch die join Funktion werden alle string zusammen geführt, die einzlenne strings sind "," getrennt

        protein_names.append(filename_base + "_B")
        contact_list.append(",".join(str(x)for x in contactsB))


    complex_contacts= pd.DataFrame()
    complex_contacts["name"]= protein_names
    complex_contacts["contacts"]= contact_list

    os.makedirs("temp", exist_ok= True)
    #einfach nur als Sicherheit damit falls der Ordner aus ungeklärlichen gründne nicht mehr exsistiert das Skript ausgeführt werden kann


    complex_contacts.to_csv("temp/contacts.tsv", sep= "\t", index= False)


if __name__ == "__main__":
    arguments()
