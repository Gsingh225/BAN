#!/usr/bin/env python3

import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

###############################################################################
# 1) PARSE THE BRACKET-AND-ARROW FORMAT
###############################################################################

def parse_molecule(text):
    """
    Parses bracket-and-arrow text (e.g., diethylamine).
    Returns a nested dict structure:
        {
          'atom': 'N',
          'charge': None,
          'substituents': [
            {
              'orientation': 'bottom',
              'bond_order': '-',
              'child': {...}  # child dict
            },
            ...
          ]
        }
    """

    for tok in ["[", "]", ",", "->"]:
        text = text.replace(tok, f" {tok} ")
    tokens = text.split()

    def consume_atom_and_block(tokens_list):
        if not tokens_list:
            raise ValueError("Ran out of tokens while expecting an atom token.")
        atom_token = tokens_list.pop(0)
        # e.g. "N", "N(+1)", "C", "lp", ...
        match_atom = re.match(r'^([A-Za-z]+)(\([+-]?\d+\))?$', atom_token)
        if not match_atom:
            raise ValueError(f"Unexpected token for atom: {atom_token}")

        base_atom = match_atom.group(1)
        maybe_charge = match_atom.group(2)
        charge = maybe_charge.strip("()") if maybe_charge else None

        node = {
            'atom': base_atom,
            'charge': charge,
            'substituents': []
        }

        # If next token is '[', parse bracket
        if tokens_list and tokens_list[0] == '[':
            tokens_list.pop(0)  # consume '['
            while tokens_list and tokens_list[0] != ']':
                if tokens_list[0] == ',':
                    tokens_list.pop(0)
                if not tokens_list:
                    raise ValueError("Unexpected end inside bracket block.")

                orientation = tokens_list.pop(0)
                if tokens_list and tokens_list[0] == '->':
                    tokens_list.pop(0)

                bond_order = None
                if tokens_list:
                    peek = tokens_list[0]
                    bond_match = re.match(r'^([-=#])(\S.*)$', peek)
                    if bond_match:
                        bond_order = bond_match.group(1)
                        remainder = bond_match.group(2)
                        tokens_list[0] = remainder

                child = consume_atom_and_block(tokens_list)
                node['substituents'].append({
                    'orientation': orientation,
                    'bond_order': bond_order,
                    'child': child
                })

            if not tokens_list or tokens_list[0] != ']':
                raise ValueError("Missing ']' after bracket block.")
            tokens_list.pop(0)  # consume ']'

        return node

    top_node = consume_atom_and_block(tokens)
    if tokens:
        raise ValueError(f"Extra tokens left after parse: {tokens}")
    return top_node

###############################################################################
# 2) BUILD AN RDKit MOL (SKIPPING LONE PAIRS)
###############################################################################

def build_rdkit_mol(parsed_node):
    """
    Build an RDKit RWMol from the parsed structure.
    - "lp" nodes (lone pairs) are skipped, as RDKit doesn't treat them as separate atoms.
    - Bond orders: '-', '=', '#'.
    - Formal charges: e.g. N(+1).

    Returns an RDKit Mol object.
    """
    rwmol = Chem.RWMol()

    # We'll assign each *non-lp* node a unique integer ID, stored in id_map[id(node_dict)].
    id_map = {}
    node_id_counter = [0]  # store in a list to mutate in nested func

    def traverse_assign_ids(node):
        """
        Recursively assigns integer IDs to each node with a real atom
        (skip 'lp'), storing in id_map[id(node)] = some_id.
        """
        # If atom is 'lp', skip giving an ID, but still traverse children (rare)
        if node['atom'].lower() == 'lp':
            for sub in node['substituents']:
                traverse_assign_ids(sub['child'])
            return

        # Give this node an ID
        this_id = node_id_counter[0]
        node_id_counter[0] += 1
        id_map[id(node)] = this_id

        # Recurse children
        for sub in node['substituents']:
            traverse_assign_ids(sub['child'])

    traverse_assign_ids(parsed_node)

    # Make sure we have enough placeholder atoms
    # The highest ID we assigned is node_id_counter[0], so we need that many atoms
    total_atoms = node_id_counter[0]
    for _ in range(total_atoms):
        rwmol.AddAtom(Chem.Atom(0))  # atomic number=0 placeholder

    def get_atom_info(atom_label, charge_str):
        """
        Convert e.g. 'N', charge_str='-1' => (atomic_num=7, formal_charge=-1).
        If unknown, default to carbon (atomic_num=6).
        If 'lp', return None => skip.
        """
        if atom_label.lower() == 'lp':
            return None, None
        # minimal periodic table
        symbol_map = {
            'H': 1, 'He': 2,
            'Li': 3, 'Be': 4,
            'B': 5, 'C': 6,
            'N': 7, 'O': 8,
            'F': 9, 'P': 15,
            'S': 16, 'Cl': 17,
            'Br': 35, 'I': 53
        }
        atomic_num = symbol_map.get(atom_label, 6)  # default carbon
        formal_charge = 0
        if charge_str is not None:
            formal_charge = int(charge_str)
        return atomic_num, formal_charge

    def add_atoms_and_bonds(node, parent_node=None, bond_order=None):
        """
        Recursively build the graph in rwmol.
        If node is 'lp', skip. If parent != None, connect them via bond_order.
        """
        if node['atom'].lower() == 'lp':
            # skip building an atom, but continue children
            for s in node['substituents']:
                add_atoms_and_bonds(s['child'], None, s['bond_order'])
            return

        # get ID for this node
        this_idx = id_map[id(node)]
        # fill in correct atomic info
        atomic_num, formal_charge = get_atom_info(node['atom'], node['charge'])
        atom = rwmol.GetAtomWithIdx(this_idx)
        atom.SetAtomicNum(atomic_num)
        atom.SetFormalCharge(formal_charge)

        if parent_node is not None:
            parent_idx = id_map[id(parent_node)]
            rd_bond = Chem.BondType.SINGLE
            if bond_order == '=':
                rd_bond = Chem.BondType.DOUBLE
            elif bond_order == '#':
                rd_bond = Chem.BondType.TRIPLE

            # add bond if not present
            if rwmol.GetBondBetweenAtoms(parent_idx, this_idx) is None:
                rwmol.AddBond(parent_idx, this_idx, rd_bond)

        # Recurse children
        for s in node['substituents']:
            add_atoms_and_bonds(s['child'], node, s['bond_order'])

    add_atoms_and_bonds(parsed_node, None, None)

    # Remove any leftover dummy atoms (atomicNum=0)
    remove_list = []
    for a in rwmol.GetAtoms():
        if a.GetAtomicNum() == 0:
            remove_list.append(a.GetIdx())
    # Must remove from highest to lowest
    for idx in sorted(remove_list, reverse=True):
        rwmol.RemoveAtom(idx)

    mol = rwmol.GetMol()
    Chem.SanitizeMol(mol)
    return mol

###############################################################################
# 3) DEMO
###############################################################################

def main():
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
    # 2) Build RDKit molecule (skipping lone pairs)
    mol = build_rdkit_mol(tree)
    # 3) Compute 2D coords
    AllChem.Compute2DCoords(mol)
    # 4) Draw
    img = Draw.MolToImage(mol, size=(300, 300))
    img.show()  # or img.save("diethylamine.png")

if __name__ == "__main__":
    main()
