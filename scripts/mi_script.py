
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

def rwr_propagation(G, seeds, restart=0.5, max_iter=50, tol=1e-6):
    nodes = list(G.nodes())
    scores = {n: 0.0 for n in nodes}

    # Asignar puntuación inicial (distribución uniforme en las semillas)
    for s in seeds:
        if s in scores:
            scores[s] = 1.0 / len(seeds)

    # Guardamos la distribución inicial
    init_scores = scores.copy()

    # Iteraciones del RWR
    for _ in range(max_iter):
        new_scores = {n: 0.0 for n in nodes}

        # Difundir la señal a los vecinos
        for n in nodes:
            neighbors = list(G.neighbors(n))
            if not neighbors:
                continue
            spread = scores[n] / len(neighbors)
            for nb in neighbors:
                new_scores[nb] += (1 - restart) * spread

        # Añadir parte del reinicio (vuelta a las semillas)
        for n in nodes:
            new_scores[n] += restart * init_scores.get(n, 0.0)

        # Comprobar convergencia
        diff = sum(abs(new_scores[n] - scores[n]) for n in nodes)
        scores = new_scores
        if diff < tol:
            break

    return scores

def diamond_like(G, seeds, k=50):
    # Módulo inicial = semillas
    module = [s for s in seeds if s in G]
    module_set = set(module)

    # Guardamos el orden
    order = {}
    step = 1
    for g in module:
        order[g] = step
        step += 1

    # Ir agregando hasta llegar a k genes
    while len(module_set) < k:
        best_node = None
        best_count = -1

        # Candidatos = nodos de la red que no están en el módulo
        for node in G.nodes():
            if node in module_set:
                continue
            # Contamos cuántos vecinos tiene dentro del módulo actual
            count = sum(1 for nb in G.neighbors(node) if nb in module_set)
            if count > best_count:
                best_count = count
                best_node = node

        # Si no hay nadie conectado, paramos
        if best_node is None or best_count == 0:
            break

        # Añadimos el mejor
        module_set.add(best_node)
        order[best_node] = step
        step += 1

    return order

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

    # Ejecutar propagación RWR
    scores = rwr_propagation(G, seeds_in_network, restart=0.5, max_iter=50)
    print("\n[OK] Propagación completada (RWR).")
    top5 = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:5]
    print("[INFO] Top 5 nodos más puntuados:")
    for g, s in top5:
        print(f"   {g}: {s:.6f}")

    # DIAMOnD simplificado
    diamond_order = diamond_like(G, seeds_in_network, k=30)
    print("\n[OK] Expansión estilo DIAMOnD completada. Primeros genes del módulo:")
    # ordenamos por rank
    first10 = sorted(diamond_order.items(), key=lambda x: x[1])[:10]
    for g, rank in first10:
        print(f"  {rank:02d}. {g}")