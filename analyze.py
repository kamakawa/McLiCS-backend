#!/usr/bin/env python3
"""
analyze.py — Script de análise dos resultados do MClist-backend
================================================================
Uso:
    python3 analyze.py                        # lê po.dat na pasta atual
    python3 analyze.py resultados/po.dat      # lê arquivo específico
    python3 analyze.py po.dat -o grafico.png  # salva em arquivo específico
    python3 analyze.py po.dat --show          # exibe o gráfico interativo
    python3 analyze.py po.dat --txt           # imprime tabela de estatísticas

Saída:
    - Gráfico com S(x), E(x) e barras de erro (variância)
    - Detecção automática da transição de fase (modos thermal e electric)
    - Tabela de estatísticas no terminal (com --txt)
"""

import sys
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator


# ── Leitura do arquivo ────────────────────────────────────────────────────────

def detect_mode(header_line):
    """Detecta o modo de evolução pelo cabeçalho do po.dat."""
    h = header_line.strip().lower()
    if h.startswith("t "):
        return "thermal", "Temperatura (T)", "T"
    elif h.startswith("e "):
        return "electric", "Campo Elétrico (E)", "E_field"
    else:
        return "step", "Passo (índice)", "step"


def read_po_dat(file_path):
    """
    Lê po.dat e retorna um DataFrame com as colunas:
      x, S, varS, E, varE
    onde x é temperatura, campo elétrico ou índice de passo, dependendo do modo.
    """
    if not os.path.isfile(file_path):
        print(f"  [ERRO] Arquivo não encontrado: {file_path}")
        sys.exit(1)

    mode, xlabel, xcol = "thermal", "Temperatura (T)", "T"
    rows = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Linha de cabeçalho
            if not line[0].replace('.','').replace('-','').isdigit():
                mode, xlabel, xcol = detect_mode(line)
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                rows.append([float(p) for p in parts[:5]])
            except ValueError:
                continue

    if not rows:
        print("  [ERRO] Nenhum dado válido encontrado em po.dat.")
        sys.exit(1)

    df = pd.DataFrame(rows, columns=["x", "S", "varS", "E", "varE"])
    df["errS"] = np.sqrt(df["varS"].clip(lower=0))
    df["errE"] = np.sqrt(df["varE"].clip(lower=0))

    if mode == "thermal":
        df = df.sort_values("x").reset_index(drop=True)

    return df, mode, xlabel, xcol


# ── Detecção de transição de fase ────────────────────────────────────────────

def find_transition(df):
    """
    Estima a temperatura (ou campo) de transição de fase como o ponto de
    maior variação absoluta de dS/dx. Retorna None se não for detectável.
    """
    if len(df) < 4:
        return None
    dS = np.abs(np.gradient(df["S"].values, df["x"].values))
    idx = np.argmax(dS)
    # Só reporta se a variação for significativa
    if dS[idx] < 0.01:
        return None
    return df["x"].iloc[idx]


# ── Estatísticas ──────────────────────────────────────────────────────────────

def print_stats(df, xlabel):
    """Imprime uma tabela de estatísticas descritivas no terminal."""
    sep = "─" * 52
    print(f"\n{sep}")
    print(f"  Estatísticas  —  {len(df)} pontos")
    print(sep)
    print(f"  {'':12s}  {'S':>10s}  {'E':>10s}")
    print(f"  {'Mínimo':12s}  {df['S'].min():10.4f}  {df['E'].min():10.4f}")
    print(f"  {'Máximo':12s}  {df['S'].max():10.4f}  {df['E'].max():10.4f}")
    print(f"  {'Média':12s}  {df['S'].mean():10.4f}  {df['E'].mean():10.4f}")
    print(f"  {'Desvio pad.':12s}  {df['S'].std():10.4f}  {df['E'].std():10.4f}")

    Tx = find_transition(df)
    if Tx is not None:
        print(f"\n  Transição de fase estimada em x ≈ {Tx:.4f}")
    print(sep + "\n")


# ── Gráfico ───────────────────────────────────────────────────────────────────

def plot(df, mode, xlabel, output_path, show):
    """Gera o gráfico de S(x) e E(x) com barras de erro."""

    # Paleta
    C_S    = "#1565C0"   # azul escuro  — parâmetro de ordem
    C_E    = "#B71C1C"   # vermelho     — energia
    C_err  = "#90CAF9"   # azul claro   — erro de S
    C_errE = "#EF9A9A"   # vermelho claro — erro de E
    C_Tx   = "#F57F17"   # laranja      — transição

    fig = plt.figure(figsize=(11, 7))
    fig.patch.set_facecolor("white")

    gs = gridspec.GridSpec(2, 1, hspace=0.08, height_ratios=[3, 1])
    ax_main = fig.add_subplot(gs[0])
    ax_err  = fig.add_subplot(gs[1], sharex=ax_main)

    x = df["x"].values

    # ── Painel principal: S e E ──────────────────────────────────────────────
    ax_main.set_facecolor("#FAFAFA")
    ax_main.grid(True, linestyle="--", linewidth=0.5, alpha=0.6, color="#CCCCCC")

    # S (eixo esquerdo)
    ax_main.plot(x, df["S"], color=C_S, linewidth=2, marker="o",
                 markersize=4, label="S  (parâmetro de ordem)", zorder=3)
    ax_main.fill_between(x,
                         df["S"] - df["errS"],
                         df["S"] + df["errS"],
                         alpha=0.18, color=C_err)
    ax_main.set_ylabel("Parâmetro de Ordem  S", color=C_S, fontsize=12)
    ax_main.tick_params(axis="y", labelcolor=C_S)
    ax_main.set_ylim(-0.05, 1.05)

    # E (eixo direito)
    ax2 = ax_main.twinx()
    ax2.plot(x, df["E"], color=C_E, linewidth=2, marker="s",
             markersize=4, linestyle="--", label="E  (energia média)", zorder=3)
    ax2.fill_between(x,
                     df["E"] - df["errE"],
                     df["E"] + df["errE"],
                     alpha=0.15, color=C_errE)
    ax2.set_ylabel("Energia Média  E", color=C_E, fontsize=12)
    ax2.tick_params(axis="y", labelcolor=C_E)

    # Linha de transição
    Tx = find_transition(df)
    if Tx is not None:
        ax_main.axvline(Tx, color=C_Tx, linewidth=1.5,
                        linestyle=":", label=f"Transição ≈ {Tx:.3f}", zorder=2)

    # Legenda unificada
    lines_a, labels_a = ax_main.get_legend_handles_labels()
    lines_b, labels_b = ax2.get_legend_handles_labels()
    ax_main.legend(lines_a + lines_b, labels_a + labels_b,
                   loc="best", fontsize=10, framealpha=0.9)

    ax_main.set_title("Resultados da Simulação de Monte Carlo — MClist",
                      fontsize=13, pad=10)
    plt.setp(ax_main.get_xticklabels(), visible=False)

    # ── Painel inferior: variância de S ──────────────────────────────────────
    ax_err.set_facecolor("#FAFAFA")
    ax_err.grid(True, linestyle="--", linewidth=0.5, alpha=0.6, color="#CCCCCC")
    ax_err.bar(x, df["varS"], width=(x[-1]-x[0])/(len(x)+1) if len(x) > 1 else 0.01,
               color=C_S, alpha=0.5, label="Var(S)")
    ax_err.set_ylabel("Var(S)", fontsize=10, color=C_S)
    ax_err.tick_params(axis="y", labelcolor=C_S)
    ax_err.set_xlabel(xlabel, fontsize=12)
    ax_err.xaxis.set_minor_locator(AutoMinorLocator())

    if Tx is not None:
        ax_err.axvline(Tx, color=C_Tx, linewidth=1.5, linestyle=":")

    fig.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"  Gráfico salvo em: {output_path}")

    if show:
        plt.show()

    plt.close()


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analisa e plota os resultados do MClist-backend (po.dat).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input", nargs="?", default="po.dat",
        help="Caminho para o arquivo po.dat (padrão: ./po.dat)"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Nome do arquivo de saída do gráfico (padrão: <input>.png)"
    )
    parser.add_argument(
        "--show", action="store_true",
        help="Exibe o gráfico interativo além de salvar"
    )
    parser.add_argument(
        "--txt", action="store_true",
        help="Imprime tabela de estatísticas no terminal"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    df, mode, xlabel, xcol = read_po_dat(args.input)

    # Nome de saída descritivo se não for especificado pelo usuário
    if args.output:
        output = args.output
    else:
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output = f"mclist_{mode}_{timestamp}.png"

    print(f"\n  MClist Analyzer")
    print(f"  {'─'*40}")
    print(f"  Entrada : {args.input}")
    print(f"  Saída   : {output}")
    print(f"  Modo    : {mode}")
    print(f"  Pontos  : {len(df)}")
    print(f"  {xcol} de {df['x'].min():.4f} até {df['x'].max():.4f}\n")

    if args.txt:
        print_stats(df, xlabel)

    plot(df, mode, xlabel, output, args.show)

    Tx = find_transition(df)
    if Tx is not None:
        print(f"  Transição de fase estimada em {xcol} ≈ {Tx:.4f}")

    print()


if __name__ == "__main__":
    main()