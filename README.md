# <div align="center">Automatização de Desenho de Fármacos</div>

<div align="center">
    <img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
    <img src="https://img.shields.io/badge/AutoDock_Vina-2E86C1?style=for-the-badge&logo=databricks&logoColor=white" alt="Vina">
    <img src="https://img.shields.io/badge/Open_Babel-8E44AD?style=for-the-badge&logo=atom&logoColor=white" alt="Open Babel">
</div>

<div align="center">
    <p><strong>Pipeline de docking molecular, redocking e análise de resultados.</strong></p>
</div>

---

## Visão Geral

Este projeto automatiza as etapas principais de modelagem molecular:

<table align="center">
    <tr>
        <td>1.</td>
        <td>Download e preparação da proteína</td>
    </tr>
    <tr>
        <td>2.</td>
        <td>Preparação do ligante</td>
    </tr>
    <tr>
        <td>3.</td>
        <td>Execução do docking</td>
    </tr>
    <tr>
        <td>4.</td>
        <td>Cálculo de RMSD</td>
    </tr>
    <tr>
        <td>5.</td>
        <td>Análise dos resultados</td>
    </tr>
</table>

---

## Módulos

<table>
    <tr>
        <th>Módulo</th>
        <th>Função</th>
        <th>Exemplo</th>
    </tr>
    <tr>
        <td><code>src/base.py</code></td>
        <td>Baixa estruturas e extrai cadeia/ligantes</td>
        <td><code>python3 src/base.py --pdb-id 9THJ --chain A --output-dir resultados/9thj</code></td>
    </tr>
    <tr>
        <td><code>src/docking.py</code></td>
        <td>Prepara PDBQT e executa o Vina</td>
        <td><code>python3 src/docking.py --pdb-id 9THJ --chain A --output-dir resultados/9thj --ligand-id A1JV9</code></td>
    </tr>
    <tr>
        <td><code>src/analysis.py</code></td>
        <td>Extrai score, calcula RMSD e gera tabela</td>
        <td><code>python3 src/analysis.py --docking-dir resultados/9thj/docking/A1JV9 --ligand-crystal resultados/9thj/ligands/A1JV9_crystal.pdb</code></td>
    </tr>
</table>

### `src/base.py`

- Baixa proteínas e ligantes do PDB.
- Faz fallback para arquivos locais do workspace.

### `src/docking.py`

- Prepara receptor e ligante em `PDBQT`.
- Calcula a caixa de docking automaticamente ou manualmente.

### `src/analysis.py`

- Lê o log do Vina.
- Calcula RMSD com alinhamento Kabsch.
- Gera saída em JSON ou tabela.

---

## Pipeline Automatizado

O script `scripts/entrega2_redocking.py` executa o fluxo completo:

```bash
python3 scripts/entrega2_redocking.py 9THJ A A1JV9 resultados/entrega2
```

### Etapas executadas

- Download da estrutura
- Preparação da proteína
- Download do ligante cristalográfico
- Docking
- Cálculo de RMSD
- Geração de resumo final

---

## Saídas Esperadas

- `resultados/9thj/raw/` — estrutura bruta
- `resultados/9thj/receptor/` — receptor preparado
- `resultados/9thj/ligands/` — ligantes baixados
- `resultados/9thj/docking/` — resultados do docking
- `resultados_redocking.json` — resumo em JSON

---

## Dependências

```bash
pip install numpy
```

Ferramentas externas:

- `obabel` — conversão de formatos
- `vina` — AutoDock Vina

---

