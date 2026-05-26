# <div align="center">Automatização de Desenho de Fármacos</div>

<div align="center">
    <img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
    <img src="https://img.shields.io/badge/AutoDock_Vina-2E86C1?style=for-the-badge&logo=databricks&logoColor=white" alt="Vina">
    <img src="https://img.shields.io/badge/AutoDock_Dock6-FF6B6B?style=for-the-badge&logo=molecular&logoColor=white" alt="Dock6">
    <img src="https://img.shields.io/badge/Open_Babel-8E44AD?style=for-the-badge&logo=atom&logoColor=white" alt="Open Babel">
</div>

<div align="center">
    <p><strong>Pipeline de pré-processamento, docking molecular (Vina & Dock6) e análise de resultados.</strong></p>
</div>

---

## Visão Geral

Este projeto automatiza as etapas principais de modelagem molecular:

| Etapa | Descrição |
|-------|-----------|
| 1. | Download e preparação da proteína (PDB/mmCIF) |
| 2. | Extração de cadeia e ligante |
| 3. | Conversão para PDBQT (Open Babel) |
| 4. | Docking com AutoDock Vina ou Dock6 |
| 5. | Cálculo de RMSD (Kabsch) |
| 6. | Identificação da melhor pose |

---

## Módulos

### Pré-processamento (`src/base.py`)

**Função:** Download e preparação de proteína/ligante do PDB

```bash
python3 src/base.py \
  --pdb-id 9THJ \
  --chain A \
  --ligand-id A1JV9 \
  --output-dir saida/
```

**Recursos:**
- Baixa proteínas do PDB em CIF ou PDB
- Extrai cadeia específica em PDB
- Extrai ligante de referência da estrutura
- Fallback para arquivos locais (`homologia_9thj/data/`, `Backup/`, `Documentos/`)
- Detecta automaticamente ligantes na estrutura
- Download de SDF ideal do PDB

**Saída:**
- `saida/raw/9thj.cif` — estrutura bruta
- `saida/receptor/9thj_A.pdb` — receptor extraído
- `saida/ligands/A1JV9_ref.pdb` — ligante de referência
- `saida/ligands/A1JV9.sdf` — ligante ideal

---

### Docking - AutoDock Vina

#### Módulo: `src/vina/redocking.py`

**Função:** Redocking completo (pré-processamento + docking + RMSD)

```python
from src.vina import redocking

best_pose, rmsd = redocking.preprocess_and_dock(
    pdb_id="9THJ",
    chain="A",
    ligand_id="A1JV9",
    output_dir="saida/",
    vina_exe="vina"
)
print(f"Melhor pose: {best_pose}")
print(f"RMSD: {rmsd:.3f} Å")
```

**Recursos:**
- Usa **apenas** pré-processamento de `src/base.py`
- Converte PDB → PDBQT (via Open Babel)
- Calcula caixa de docking automaticamente (baseada no ligante)
- Executa AutoDock Vina
- Calcula RMSD com alinhamento Kabsch
- Extrai e retorna a melhor pose

#### Pipeline: `pipeline/vina_pipeline.py`

**Uso:**
```bash
python3 pipeline/vina_pipeline.py 9THJ A A1JV9 saida/
```

**Output no terminal:**
```
Usando receptor pré-preparado do Backup: ...
Melhor pose: saida/vina_out/best_pose_9thj_A_A1JV9.pdbqt
RMSD da melhor pose: 3.775 Å
```

---

### Docking - AutoDock Dock6

#### Pipeline: `pipeline/dock6_unificado.py`

**Uso:**
```bash
python3 pipeline/dock6_unificado.py redocking 9THJ A A1JV9 saida/ \
  --receptor-ms Backup/dock6/rec.ms \
  --dock6-param-dir Backup/dock6
```

**Módulos suportados:**
- `src/dock6/docking.py` — executa docking
- `src/dock6/analysis.py` — análise de resultados

---

## Estrutura de Saída

Resultado típico do pipeline Vina:

```
saida/
├── raw/
│   └── 9thj.cif                          # Estrutura bruta do PDB
├── receptor/
│   ├── 9thj_A.pdb                        # Receptor extraído (PDB)
│   └── 9thj_A.pdbqt                      # Receptor preparado (PDBQT)
├── ligands/
│   ├── A1JV9_ref.pdb                     # Ligante referência (extraído)
│   ├── A1JV9_ref.pdbqt                   # Ligante preparado (PDBQT)
│   └── A1JV9.sdf                         # Ligante ideal (PDB)
└── vina_out/
    ├── 9thj_A_A1JV9_out.pdbqt            # Todas as poses do Vina
    └── best_pose_9thj_A_A1JV9.pdbqt      # Melhor pose
```

---

## Dependências

### Python

```bash
pip install numpy
```

### Ferramentas externas

- **`obabel`** (Open Babel) — conversão PDB ↔ PDBQT
- **`vina`** (AutoDock Vina 1.2.3+) — docking molecular
- **`dock6`** (AutoDock Dock6) — docking alternativo

### Instalação em Ubuntu/Debian

```bash
sudo apt-get install openbabel autodock-vina autodock
```

### Nota sobre processamento de receptor

- Para **9THJ**, usa receptor pré-preparado do `Backup/vina/receptor.pdbqt`
- Para outros PDBs, gera PDBQT via Open Babel (pode precisar de `meeko` ou `AutoDockTools` para máxima compatibilidade)

---

## Formato de Saída - RMSD

O RMSD é calculado usando **alinhamento Kabsch (SVD)**:

1. Centraliza estruturas referência e pose
2. Computa matriz de covariância
3. Encontra melhor rotação via SVD (usando numpy)
4. Calcula RMSD após alinhamento ótimo

**Fallback:** Se numpy não estiver disponível, calcula RMSD sem rotação ótima (valor aproximado por cima).

---

## Exemplos

### 1. Pré-processamento apenas

```bash
python3 src/base.py \
  --pdb-id 9THJ \
  --chain A \
  --ligand-id A1JV9 \
  --output-dir saida/
```

Gera estrutura + ligante em `saida/receptor/` e `saida/ligands/`.

### 2. Redocking completo (Vina)

```bash
python3 pipeline/vina_pipeline.py 9THJ A A1JV9 saida/
```

Output no terminal:
```
Usando receptor pré-preparado do Backup: Backup/vina/receptor.pdbqt
Melhor pose: saida/vina_out/best_pose_9thj_A_A1JV9.pdbqt
RMSD da melhor pose: 3.775 Å
```

### 3. Redocking completo (Dock6)

```bash
python3 pipeline/dock6_unificado.py redocking 9THJ A A1JV9 saida/ \
  --receptor-ms Backup/dock6/rec.ms \
  --dock6-param-dir Backup/dock6
```

---

## Diretórios do Projeto

| Diretório | Conteúdo |
|-----------|----------|
| `src/` | Código principal (base.py, vina/, dock6/) |
| `src/base.py` | Pré-processamento (download, extração) |
| `src/vina/` | Módulo de redocking com Vina |
| `src/dock6/` | Módulo de redocking com Dock6 |
| `pipeline/` | Scripts de pipeline executáveis |
| `homologia_9thj/` | Modelagem por homologia (MODELLER) |
| `Backup/` | Arquivos pré-processados de referência |
| `Documentos/` | Documentação e referências |

---

## Notas Técnicas

### Processamento de Receptor (Vina)

O módulo `src/vina/redocking.py` usa o receptor pré-preparado do Backup quando disponível para 9THJ. Isso contorna problemas de compatibilidade com tipos de átomos PDBQT gerados por Open Babel.

Para máxima compatibilidade com outros PDBs, recomenda-se usar `meeko` ou `AutoDockTools`:

```bash
# Com meeko (alternativa moderna)
mk_prepare_receptor.py -i receptor.pdb -o receptor.pdbqt

# Ou com AutoDockTools (clássico)
prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
```

### Cálculo de RMSD

O RMSD entre pose e ligante de referência é calculado:

1. Extraindo coordenadas de ambos os arquivos PDB/PDBQT
2. Alinhando estruturalmente (Kabsch/SVD)
3. Calculando distância RMSD após alinhamento ótimo

**Limitações:**
- Se a ordem de átomos diferir (ex: hidrogênios adicionados), faz ajuste automático
- Para análise rigorosa, mapear átomos por nome/elemento

---

## Autores e Licença

Projeto de Automatização em Desenho de Fármacos - 2026

---
