import random


def convert_directed_to_undirected_weighted(input_file, output_file, weight_range=(1, 10)):
    edge_set = set()
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue  # Skip comments and blank lines
            u, v = map(int, line.strip().split())
            edge = tuple(sorted((u, v)))
            edge_set.add(edge)

    # Generate weights and sort edges
    weighted_edges = [(u, v, random.randint(*weight_range)) for u, v in edge_set]
    weighted_edges.sort()  # Sort by (u, v)

    with open(output_file, 'w') as outfile:
        for u, v, weight in weighted_edges:
            outfile.write(f"{u} {v} {weight}\n")

    print(f"Converted and sorted {len(weighted_edges)} edges.")

# Example usage:
convert_directed_to_undirected_weighted("directed_graph.txt", "undirected_weighted_graph.txt")