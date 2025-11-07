
import pandas as pd
import networkx as nx


def load_network(path):

    # Intentamos leer con cabecera primero
    df = pd.read_csv(path, sep="\t")

    if {"protein1_hugo", "protein2_hugo"}.issubset(df.columns):
        # Caso de la red STRING en HUGO
        df = df[["protein1_hugo", "protein2_hugo"]]
        df.columns = ["u", "v"]
        G = nx.from_pandas_edgelist(df, "u", "v")
        return G

    # Leemos el archivo por tabuladores o por espacios
    df = pd.read_csv(path, sep=None, engine="python", header=None)

    # 2 columnas -> red sin pesos
    if df.shape[1] == 2:
        df.columns = ["u", "v"]
        G = nx.from_pandas_edgelist(df, "u", "v")
    # 3 columnas -> red con pesos
    else:
        df.columns = ["u", "v", "w"]
        G = nx.from_pandas_edgelist(df, "u", "v", edge_attr="w")

    return G

def load_seeds(path):
    # Lee un archivo de semillas y devuelve una lista de strings
    seeds = []
    with open(path) as f:
        for line in f:
            gene = line.strip()
            if gene:       
                seeds.append(gene)
    return seeds

if __name__ == "__main__":
    network_path = network_path = "data/string_network_filtered_hugo-400.tsv"
    seeds_path = "data/genes_seed.txt"


    # Cargar la red
    try:
        G = load_network(network_path)
        print("[OK] Red cargada desde:", network_path)
        print(f"[INFO] Nodos: {G.number_of_nodes()} | Aristas: {G.number_of_edges()}")
    except FileNotFoundError:
        print(f"[ERROR] No se encontró el archivo de red: {network_path}")
        raise SystemExit(1)
    except Exception as e:
        print(f"[ERROR] No se pudo cargar la red: {e}")
        raise SystemExit(1)

    # Cargar las semillas
    try:
        seeds_raw = load_seeds(seeds_path)
        print(f"[OK] Semillas cargadas desde: {seeds_path}")
        print(f"[INFO] Semillas leídas: {len(seeds_raw)} -> {seeds_raw}")
    except FileNotFoundError:
        print(f"[ERROR] No se encontró el archivo de semillas: {seeds_path}")
        raise SystemExit(1)
    except Exception as e:
        print(f"[ERROR] No se pudieron cargar las semillas: {e}")
        raise SystemExit(1)
    
    # Comprobar cuántas semillas hay en la red
    seeds_in_network = [s for s in seeds_raw if s in G]
    missing_seeds = [s for s in seeds_raw if s not in G]

    print(f"[INFO] Semillas presentes en la red: {len(seeds_in_network)} / {len(seeds_raw)}")
    if missing_seeds:
        print(f"[WARN] Estas semillas no están en la red: {missing_seeds}")