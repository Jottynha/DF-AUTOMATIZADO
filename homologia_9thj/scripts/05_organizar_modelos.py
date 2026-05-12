import glob, shutil, os

os.makedirs("resultados/modelos", exist_ok=True)
for file in glob.glob("target.B*.pdb"):
    destino = os.path.join("resultados/modelos", file)
    shutil.move(file, destino)
    print("Movido:", file, "->", destino)

