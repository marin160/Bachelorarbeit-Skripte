# Hier soll die Domäne gefunden werden die die meisten Kontaktsellen beinhlatet und von dieser soll dann die Seqeunz als FASTA ausgegeben werden
import os
import sys
import pandas as pd
from pymol import cmd, stored


#Hier sind Funktionen definiert
#Funktion die den cuttingstring in die Unterdomänen trennt
def split_domainranges(cutting_string):
    domains= cutting_string.split(',')
    #Teilt z.B diesen String '5-51,86-443_588-659,800-832' an jedem Komma
    domain_list = []
    for domain in domains: 
        subranges = domain.split('_')
        subdomain= []
        #Hier wird jede Subdomäne durch einen '_' geteeilt wie z.B 86-443_588-659
        for sub in subranges:
            start, end = map(int, sub.split('-'))
            #'-' extrahiert den Startpunkt und das Ende von jeder Subdomäne
            subdomain.append((start,end))
        domain_list.append(subdomain)
       #1 Teilergebnis:  [(5-51)], [(86-443),(588-659)],[(800-832)] drei Domänen
    return domain_list

#Die Funktion dient dazu um zu schauen wieviele Kontakt Reste sich in einer Domäne befinden und zählt diese   
def count_contacts_in_domains(contact_residues, domain_list):
    result_count =[]
    for domain in domain_list: 
        count= 0
        for start, end in domain:
            for contact_residue in contact_residues: 
                if start <= contact_residue <= end: 
                    count +=1 
                    #Es wird geschaut für jede Nummer des Kontakt-Rests ob dieser zwischem dem Start und Endpunkt einer Domäne liegt 
                    #Beispiel: contact_residues =[37,,44,62] domain=(86,443)
                    
        result_count.append(count)
    return result_count
    # 2 Teilergebnis: Es wurde eine Liste erstellt die, die Anzahl der Kontakte pro Domäne hat
    
#Diese Funktion schaut welche Domäne die  beste ist, beduetet welche Domäne die meisten Kontaktstellen aufweist
def best_domain(contact_residues, cutting_string):
    domain_list = split_domainranges(cutting_string)
    counts= count_contacts_in_domains(contact_residues, domain_list)
    # Es wird zuerst die Funktion split_domainranges aufgerufen damit die Domänstrukutur vorhanden ist
    # danach wird die count_contacts_in_domain aufgerufenum die Kontakte der Reste pro Domäne zählen zu können
    # counts= [2,0,0,0]
    if max(counts) == 0:
        print("ERROR: NO CONTACT DOMAIN FOUND")
        return None
    else: 
        best_index= counts.index(max(counts))
        print(f"Contact amount: {counts[best_index]}")
        print(f"Best Domain: {best_index+1}")
        return domain_list[best_index], counts[best_index], best_index+1
        # es wird nach der höchsten anzahl an Kontakten gesucht. Falls keine gefunden werden, gibt es eine Fehlermeldung
        # falls es kontaktstellen in den Domänen gibt, wird die Domäne mit den meisten Kontakten wiedergegeben und die Anzahl wieviel wird auch angegeben

#contact_A = [31, 32, 34, 37, 40, 78, 79]
#split_A = "5-51,86-443_588-659,800-832"
#print(best_domain(contact_A, split_A))
# ergebniss mit diesem Beispiel wäre ([(5, 51)], 5) bedeutet die beste Domain ist die von 5-51 und es gibt 5 Kontaktstellen in dieser (31,32,34,37, 40)


#Diese Funktion extrahiert die Seqeunz der besten Domäne und gibt die Seqeunz als FASTA Format aus
def domain_to_expression(domain):
    # domain: z.B. [(5,7)] oder [(86,443),(588,659)]
    partial_expressions = []
    for subdomain in domain:
        start_residue = subdomain[0]
        end_residue   = subdomain[1]

        if start_residue == end_residue:
            partial_expressions.append(f"resi {start_residue}") 
        else:
            partial_expressions.append(f"resi {start_residue}-{end_residue}") 

    resi_expression = " + ".join(partial_expressions)  
    return resi_expression
    #diese Funktion erstell einen String der als PyMol Kommand genutzt werden kann um den Bereich der Domäne auswählen zu können




#Diese Funktion soll bei der ausfürhung der Seqeunzen ab Anfang und am Ende 10 Aminosäuren mehr dran hängen.
def domain_extension(domain, extend_by=10):
    domain_start= min(start for start, end in domain)
    domain_end= max(end for start, end in domain)
    new_start = max(1, domain_start - extend_by)
    new_end   = domain_end + extend_by
    return [(new_start, new_end)]


#Diese funktion diesnt dazu um .fasta dateien zu erstellen und die seqeunz in diese abzuspeichern  
def get_fasta(input):
    path_to_con_tsv="temp/contacts.tsv"
    path_to_domain_tsv= "temp/domains.tsv"
    structure_folder = input
    target_folder = "fasta_files"
    
    contacts_df = pd.read_csv("temp/contacts.tsv", sep="\t")  # TSV laden [web:9]
    contdict = {} 
    for index, row in contacts_df.iterrows(): 
        name = row["name"]            
        contacts_str = str(row["contacts"])
        contacts_list = [int(x) for x in contacts_str.split(",")]  
        contdict[name] = contacts_list

    domains_df = pd.read_csv("temp/domains.tsv", sep="\t") 
    domains_dict = {}
    for index, row in domains_df.iterrows():
        raw_name = row["name"]                  
        name = raw_name[:-1] + "_" + raw_name[-1] 
        domains_dict[name] = row["domain_cutoffs"]
    
    best_domains_path = os.path.join(target_folder, "best_domains.txt")
    os.makedirs(target_folder, exist_ok=True)
    best_domains_file = open(best_domains_path, "w")

    best_domains_file.write("#name\tbest_domain_extended\tbest_domain_original\tcontact_count\tdomain_index\n")

    string_out= ""

    for name in domains_dict:  # z.B. fold_muc1xubb_human_model_0A
        base = name[:-2]             # fold_muc1xubb_human_model_0
        chain_id = name[-1]            # 'A' oder 'B' 
    
        cif_path = os.path.join(structure_folder, base + ".cif")
        pdb_path = os.path.join(structure_folder, base + ".pdb")

        if os.path.exists(cif_path):
            struct_path = cif_path
        elif os.path.exists(pdb_path):
            struct_path = pdb_path
        else:
            print(f"ERROR: No CIF/PDB-Data found for {name}. It will be skipped")
            continue
   
        
        if name not in contdict:
            print(f"ERROR: No contacts found for {name} in contacts.tsv. It will be skipped")
            continue



        
        cmd.reinitialize()
        cmd.load(struct_path, base)

        cmd.select("chain_sel", f"{base} and chain {chain_id}")

        cutting_string = domains_dict[name]
        result_domain, contact_count, best_index = best_domain(contdict[name], cutting_string)

        extended_domain = domain_extension(result_domain, extend_by=10)

        selection_string = f"chain_sel and ({domain_to_expression(extended_domain)})"
        fasta_str = cmd.get_fastastr(selection_string)

        data_name = f"{name}.fasta"
        #data_name = f"{name}"
        fasta_file_string =fasta_str
        comp_path = os.path.join(target_folder, data_name)

        os.makedirs(target_folder, exist_ok=True)

        with open(comp_path, "w") as f:
            f.write(fasta_file_string)

        print("FASTA-Files saved unter:", comp_path)

        string_out = string_out + f"Name of the complex: {base} chain: {chain_id}  The best domain is the {best_index}.  Domain range {result_domain} and number of contacts {contact_count}\n"
        string_out = string_out + f"Sequence:\n{fasta_file_string}\n"
        report_path = os.path.join(target_folder, "best_domains_report.txt")
        with open(report_path, "w") as rep:
            rep.write(string_out)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("ERROR:python domain_extract.py <structure_folder>")
        sys.exit(1)
     
    structure_folder = sys.argv[1]
    if not os.path.isdir(structure_folder):
        print(f"ERROR: {structure_folder} is not a directory")
        sys.exit(1)
     
    get_fasta(structure_folder)

