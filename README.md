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

## Docking

Depois de baixar o receptor e o ligante com `base.py`, você pode executar o docking com `src/docking.py`:

```bash
python3 src/docking.py --pdb-id 9THJ --chain A --output-dir resultados/9thj --ligand-id A1JV9
```

### Saídas do docking

- `resultados/9thj/docking/<ligante>/receptor.pdbqt`: receptor preparado para o Vina
- `resultados/9thj/docking/<ligante>/ligand.pdbqt`: ligante preparado para o Vina
- `resultados/9thj/docking/<ligante>/config.txt`: configuração usada no Vina
- `resultados/9thj/docking/<ligante>/out.pdbqt`: poses geradas pelo Vina
- `resultados/9thj/docking/<ligante>/log.txt`: saída textual da execução

