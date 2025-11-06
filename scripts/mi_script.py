
import pandas as pd
import networkx as nx


def load_network(path):

    # Leemos el archivo por tabuladores o por espacios
    try:
        df = pd.read_csv(path, sep="\t", header=None)
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", header=None)

    # 2 columnas -> red sin pesos
    if df.shape[1] == 2:
        df.columns = ["u", "v"]
        G = nx.from_pandas_edgelist(df, "u", "v")
    # 3 columnas -> red con pesos
    else:
        df.columns = ["u", "v", "w"]
        G = nx.from_pandas_edgelist(df, "u", "v", edge_attr="w")

    return G


if __name__ == "__main__":
    network_path = "data/network_guild.txt"

    try:
        G = load_network(network_path)
        print("[OK] Red cargada desde:", network_path)
        print(f"[INFO] Nodos: {G.number_of_nodes()} | Aristas: {G.number_of_edges()}")
    except FileNotFoundError:
        print(f"[ERROR] No se encontró el archivo de red: {network_path}")
    except Exception as e:
        print(f"[ERROR] No se pudo cargar la red: {e}")
