import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

DATA_FILE = "po.dat"
OUTPUT_PLOT_FILE = "grafico.png"

def read_po_dat(file_path):
    """
    Lê o arquivo po.dat gerado pelo MClist, ignorando linhas de metadados
    e cabeçalho, e retorna um DataFrame com colunas T, S e E.
    """
    print(f"Lendo dados de {file_path}...")

    data_lines = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                # Ignora linhas de metadados e o cabeçalho
                if not line or line.startswith('[') or 'T S varS E varE' in line:
                    continue

                parts = line.split()
                # Aceita apenas linhas com pelo menos 5 colunas (T, S, varS, E, varE)
                if len(parts) >= 5:
                    try:
                        T = float(parts[0])
                        S = float(parts[1])
                        E = float(parts[3])
                        data_lines.append([T, S, E])
                    except ValueError:
                        continue

    except FileNotFoundError:
        print(f"❌ Erro: Arquivo {file_path} não encontrado.")
        return None

    if not data_lines:
        print("❌ Erro: Nenhuma linha de dados válida encontrada no arquivo.")
        return None

    df = pd.DataFrame(data_lines, columns=['T', 'S', 'E'])
    print(f"✅ Leitura concluída. {len(df)} pontos carregados.")
    print(df.head())
    return df


def plot_thermodynamics(df):
    """
    Plota o Parâmetro de Ordem (S) e a Energia Média (E) em função da Temperatura (T).
    """
    # Corrigido: temperatura em ordem crescente
    df = df.sort_values(by='T', ascending=True).reset_index(drop=True)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set_title("Validação Termodinâmica da Simulação de Monte Carlo", fontsize=14)
    ax1.set_xlabel("Temperatura (T)", fontsize=12)

    # Parâmetro de Ordem (S)
    ax1.set_ylabel("Parâmetro de Ordem Médio (S)", color='blue', fontsize=12)
    ax1.plot(df['T'], df['S'], marker='o', linestyle='-', color='blue', label='S vs T')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Energia Média (E)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Energia Média (E)", color='red', fontsize=12)
    ax2.plot(df['T'], df['E'], marker='s', linestyle='--', color='red', label='E vs T')
    ax2.tick_params(axis='y', labelcolor='red')

    # Legendas combinadas
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

    fig.tight_layout()
    plt.savefig(OUTPUT_PLOT_FILE, dpi=300)
    print(f"\n✅ Gráfico salvo como '{OUTPUT_PLOT_FILE}'")

def main():
    df_data = read_po_dat(DATA_FILE)
    if df_data is not None:
        plot_thermodynamics(df_data)

if __name__ == "__main__":
    main()
