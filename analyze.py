import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

DATA_DIR = "."

def read_simulation_parameters(file_path="input_parameters.txt"):
    params = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or '=' not in line:
                    continue
                key, value = line.split('=', 1)
                key, value = key.strip(), value.strip()
                if key == "Ti":
                    params['Ti'] = float(value)
                elif key == "Tf":
                    params['Tf'] = float(value)
                elif key == "dT":
                    params['dT'] = float(value)
    except FileNotFoundError:
        print("⚠️ input_parameters.txt não encontrado. Usando padrão.")
        params['Ti'], params['Tf'], params['dT'] = 1.2, 0.0, -0.02
    return params

def compute_S_from_Q(df):
    """Calcula o S global a partir de nx,ny,nz usando o tensor Q."""
    if not all(c in df.columns for c in ('nx','ny','nz')):
        return None
    vecs = df[['nx','ny','nz']].to_numpy(dtype=float)
    N = vecs.shape[0]
    if N == 0:
        return None
    Q = np.zeros((3,3), dtype=float)
    for u in vecs:
        norm = np.linalg.norm(u)
        if norm > 0:
            u = u / norm
        Q += 1.5 * np.outer(u, u) - 0.5 * np.eye(3)
    Q /= N
    vals = np.linalg.eigvalsh(Q)
    return float(vals[-1])  # maior autovalor

def main():
    sim_params = read_simulation_parameters()
    if not all(k in sim_params for k in ('Ti','Tf','dT')):
        print("❌ Erro: parâmetros Ti, Tf, dT faltando.")
        return

    files = sorted(glob.glob(os.path.join(DATA_DIR, "director_field_*.csv")))
    if not files:
        print("❌ Nenhum arquivo director_field*.csv encontrado.")
        return

    temps = []
    S_values = []

    for k, file_path in enumerate(files):
        df = pd.read_csv(file_path)
        temperature = sim_params['Ti'] + k * sim_params['dT']

        # respeita intervalo [Ti, Tf]
        if sim_params['dT'] < 0 and temperature < sim_params['Tf']:
            break
        if sim_params['dT'] > 0 and temperature > sim_params['Tf']:
            break

        S_Q = compute_S_from_Q(df)
        if S_Q is None:
            print(f"⚠️ Não foi possível calcular S em {file_path}. Pulando.")
            continue

        temps.append(temperature)
        S_values.append(S_Q)

    if not temps:
        print("❌ Nenhum dado válido coletado.")
        return

    # ordena por temperatura crescente
    temps, S_values = zip(*sorted(zip(temps, S_values)))

    print("\nResultados (Temperatura, S):")
    for t, s in zip(temps, S_values):
        print(f"{t:.6f}\t{s:.6f}")

    # Plot
    plt.figure(figsize=(10,6))
    plt.plot(temps, S_values, marker='o', linestyle='-')
    plt.xlabel("Temperatura (T)")
    plt.ylabel("Parâmetro de Ordem Médio (S)")
    plt.title("Transição de Fase — S vs T")
    plt.grid(True)
    plt.axhline(0.5, linestyle='--', linewidth=0.8)

    # eixo x invertido se dT < 0
    if sim_params['dT'] < 0:
        plt.gca().invert_xaxis()

    plt.tight_layout()
    plt.savefig("grafico.png")
    print("✅ Gráfico salvo como grafico.png")

if __name__ == "__main__":
    main()
