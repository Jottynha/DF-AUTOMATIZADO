# Automatização de Desenho de Fármacos

Base inicial para automatizar o download de proteínas e ligantes a partir do PDB. O script tenta baixar a estrutura do RCSB e, se não encontrar, usa uma cópia local do workspace quando existir.

## Como Executar
```bash
python3 src/base.py --pdb-id 9THJ --chain A --output-dir resultados/9thj
```

### Saídas geradas
- `resultados/9thj/raw/9thj.cif` ou `resultados/9thj/raw/9thj.pdb`: estrutura bruta baixada/recuperada
- `resultados/9thj/receptor/9thj_A.pdb`: cadeia extraída da proteína
- `resultados/9thj/ligands/*.sdf`: ligantes detectados ou informados manualmente

### Ligantes manuais

Se quiser informar ligantes específicos:
```bash
python3 src/base.py --pdb-id 9THJ --chain A --output-dir resultados/9thj --ligand-id 9TN
```

