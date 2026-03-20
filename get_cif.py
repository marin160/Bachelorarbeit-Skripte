import os
import shutil
import sys

def copy_cif_from_candidates(candidates_root):
    # Ordner, in dem das Skript liegt (Ziel)
    script_dir = os.path.dirname(os.path.abspath(__file__))

    kopiert = 0
    for root, dirs, files in os.walk(candidates_root):
        for file in files:
            if file.lower().endswith(".cif"):
                src = os.path.join(root, file)
                dst = os.path.join(script_dir, file)

                # Duplikat-Handling: _1, _2 etc.
                if os.path.exists(dst):
                    base, ext = os.path.splitext(file)
                    i = 1
                    while True:
                        new_name = f"{base}_{i}{ext}"
                        dst = os.path.join(script_dir, new_name)
                        if not os.path.exists(dst):
                            break
                        i += 1

                shutil.copy2(src, dst)  # Nur kopieren!
                print(f"Kopiert: {os.path.basename(src)} -> {os.path.basename(dst)}")
                kopiert += 1

    print(f"\nFertig! {kopiert} .cif-Dateien in '{script_dir}' kopiert.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        candidates_root = input("Pfad zum 'canidates'-Ordner eingeben: ").strip()
    else:
        candidates_root = sys.argv[1]

    if not os.path.exists(candidates_root):
        print(f"Fehler: '{candidates_root}' existiert nicht!")
    else:
        copy_cif_from_candidates(candidates_root)
