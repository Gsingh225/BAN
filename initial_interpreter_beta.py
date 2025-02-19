#!/usr/bin/env python3

import re
import networkx as nx
import matplotlib.pyplot as plt

def parse_molecule(text):
    """
    Parses the bracket-and-arrow text format into a nested Python dict.
    
    Example input (Diethylamine):
    N [ 
      bottom -> lp,
      left   -> -H,
      top    -> -C [
          left  -> -H,
          right -> -H,
          top   -> -C [
              left  -> -H,
              right -> -H,
              top   -> -H
          ]
      ],
      right  -> -C [
          left  -> -H,
          right -> -H,
          top   -> -C [
              left  -> -H,
              right -> -H,
              top   -> -H
          ]
      ]
    ]

    Returns a dict of the form:
    {
      'atom': 'N',
      'charge': None,
      'substituents': [
        {
          'orientation': 'bottom',
          'bond_order': None,
          'child': {
            'atom': 'lp',
            'charge': None,
            'substituents': []
          }
        },
        ...
      ]
    }
    """

    # Pre-process the text to separate tokens
    prepped = text
    for tok in ["[", "]", ",", "->"]:
        prepped = prepped.replace(tok, f" {tok} ")
    tokens = prepped.split()

    def consume_atom_and_block(token_list):
        """
        Recursively parse something like:
           N [ top -> -C [ ... ], bottom -> lp ]
        """
        if not token_list:
            raise ValueError("Ran out of tokens while expecting an atom token!")

        # First token might be "N", "N(+1)", "C", "lp", "H", etc.
        atom_token = token_list.pop(0)
        match_atom = re.match(r'^([A-Za-z]+)(\([+-]?\d+\))?$', atom_token)
        if not match_atom:
            raise ValueError(f"Unexpected token for atom: {atom_token}")

        base_atom = match_atom.group(1)       # e.g. "N", "C", "lp" ...
        maybe_charge = match_atom.group(2)    # e.g. "(+1)" or None
        charge = maybe_charge.strip("()") if maybe_charge else None

        node = {
            'atom': base_atom,
            'charge': charge,
            'substituents': []
        }

        # If next token is '[', parse the bracket block
        if token_list and token_list[0] == '[':
            token_list.pop(0)  # consume '['
            while token_list and token_list[0] != ']':
                if token_list[0] == ',':
                    token_list.pop(0)  # skip commas

                if not token_list:
                    raise ValueError("Unexpected end inside bracket block.")

                # orientation or label
                orientation = token_list.pop(0)

                # If the next token is '->', consume it
                if token_list and token_list[0] == '->':
                    token_list.pop(0)

                # Possibly a bond symbol at the start of the next token (e.g. '-C' or '=O' or '#N')
                bond_order = None
                if token_list:
                    peek = token_list[0]
                    bond_match = re.match(r'^([-=#])(\S.*)$', peek)
                    if bond_match:
                        bond_order = bond_match.group(1)  # '-', '=', '#'
                        remainder = bond_match.group(2)   # e.g. 'C', 'N(+1)', ...
                        # Replace the token with remainder so we can parse that as an atom token
                        token_list[0] = remainder

                # Recursively parse child
                child = consume_atom_and_block(token_list)

                node['substituents'].append({
                    'orientation': orientation,
                    'bond_order': bond_order,
                    'child': child
                })

            # We expect a closing ']'
            if not token_list or token_list[0] != ']':
                raise ValueError("Missing ']' after bracket block.")
            token_list.pop(0)  # consume ']'

        return node

    top_node = consume_atom_and_block(tokens)

    if tokens:
        raise ValueError(f"Extra tokens left after parse: {tokens}")

    return top_node


def build_nx_graph(node, g=None):
    """
    Recursively build a NetworkX graph from the nested parse dict.
    Each atom or 'lp' is a graph node; each substituent is an edge.

    Returns: (G, node_id)
      where G is the final NetworkX graph,
            node_id is the integer ID for this node in the graph.
    """
    if g is None:
        g = nx.Graph()

    # Create a new node ID
    node_id = len(g.nodes)

    # Build a label that includes any charge
    label = node['atom']
    if node['charge']:
        label += f"({node['charge']})"

    # Add the node with attributes
    g.add_node(node_id, label=label, atom=node['atom'], charge=node['charge'])

    # Recurse over substituents
    for sub in node['substituents']:
        bond_order = sub['bond_order'] if sub['bond_order'] else '-'
        child_dict = sub['child']
        # Build the child as well
        g, child_id = build_nx_graph(child_dict, g)
        # Add an edge for this substituent
        g.add_edge(node_id, child_id, bond_order=bond_order, orientation=sub['orientation'])

    return g, node_id


def draw_molecule(g):
    """
    Draw the molecule using NetworkX + Matplotlib.
    We use a force-directed layout (spring_layout) to avoid collisions.
    Double/triple bonds are drawn thicker for demonstration.
    """
    # Compute 2D coordinates (you can also try kamada_kawai_layout, etc.)
    pos = nx.spring_layout(g, k=0.7, iterations=100, seed=42)

    # Separate edges by bond order for different styling
    single_edges = []
    double_edges = []
    triple_edges = []
    for (u, v, edata) in g.edges(data=True):
        bo = edata.get('bond_order', '-')
        if bo == '-':
            single_edges.append((u, v))
        elif bo == '=':
            double_edges.append((u, v))
        elif bo == '#':
            triple_edges.append((u, v))
        else:
            # fallback
            single_edges.append((u, v))

    # Draw edges
    plt.figure(figsize=(8, 6))
    nx.draw_networkx_edges(g, pos, edgelist=single_edges, width=1.0, edge_color='k')
    nx.draw_networkx_edges(g, pos, edgelist=double_edges, width=2.0, edge_color='k')
    nx.draw_networkx_edges(g, pos, edgelist=triple_edges, width=3.0, edge_color='k')

    # Draw nodes
    node_labels = nx.get_node_attributes(g, 'label')
    nx.draw_networkx_nodes(g, pos, node_color='lightblue', node_size=600)
    nx.draw_networkx_labels(g, pos, labels=node_labels, font_size=10)

    plt.title("Molecule via Force-Directed Layout")
    plt.axis('off')
    plt.show()


def main():
    # Diethylamine example
    diethylamine_text = """
    N [ 
      bottom -> lp,
      left   -> -H,
      top    -> -C [
          left  -> -H,
          right -> -H,
          top   -> -C [
              left  -> -H,
              right -> -H,
              top   -> -H
          ]
      ],
      right  -> -C [
          left  -> -H,
          right -> -H,
          top   -> -C [
              left  -> -H,
              right -> -H,
              top   -> -H
          ]
      ]
    ]
    """

    # 1) Parse
    tree = parse_molecule(diethylamine_text)

    # 2) Build the NetworkX graph
    g, root_id = build_nx_graph(tree)

    # 3) Draw the molecule with a standard 2D layout
    draw_molecule(g)


if __name__ == "__main__":
    main()
