"""Módulos organizados para execução de DOCK6.

Camadas principais:
- prepare: preparação química/geométrica das entradas.
- docking: execução de baixo nível do DOCK6.
- analysis: extração de score/RMSD.
- redocking: fluxo completo para redocking.
- dock_library: triagem de bibliotecas de ligantes.
"""

__all__ = [
    "analysis",
    "dock_library",
    "docking",
    "prepare",
    "redocking",
]
