#!/usr/bin/env python3

import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

###############################################################################
# 1) PARSE THE BRACKET-AND-ARROW FORMAT
###############################################################################

def parse_molecule(text):
    """
    Parses bracket-and-arrow text (e.g., diethylamine) into a nested dict:
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
    # Insert spacing around special tokens
    for tok in ["[", "]", ",", "->"]:
        text = text.replace(tok, f" {tok} ")
    tokens = text.split()

    def consume_atom_and_block(tokens_list):
        if not tokens_list:
            raise ValueError("Ran out of tokens while expecting an atom token.")

        atom_token = tokens_list.pop(0)
        # Example: "N", "N(+1)", "C", "lp", ...
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

        # If next token is '[', parse bracket content
        if tokens_list and tokens_list[0] == '[':
            tokens_list.pop(0)  # consume '['
            while tokens_list and tokens_list[0] != ']':
                if tokens_list[0] == ',':
                    tokens_list.pop(0)  # skip comma

                if not tokens_list:
                    raise ValueError("Unexpected end inside bracket block.")

                orientation = tokens_list.pop(0)  # e.g. 'bottom','top','left','right'
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
# 2) BUILD AN RDKit MOL, INCLUDING "lp" AS DUMMY ATOMS
###############################################################################

def build_rdkit_mol(parsed_node):
    """
    Build an RDKit RWMol from the parsed structure.
    - Real atoms (C, H, N, O, etc.) get normal atomic numbers.
    - 'lp' is a dummy atom (atomic_num=0) labeled "lp", with a single bond to its parent.
    - Bond orders: '-', '=', '#'.
    - Formal charges: e.g. N(+1).
    
    Because dummy atoms with real bonds break normal valence rules, we
    skip standard sanitization (use SANITIZE_NONE). So this is purely
    a hack for visualization.
    """

    rwmol = Chem.RWMol()

    # We'll assign each node an integer ID, stored in id_map[id(node_dict)].
    id_map = {}
    node_id_counter = [0]

    def traverse_assign_ids(node):
        # Give every node, including 'lp', an ID
        this_id = node_id_counter[0]
        node_id_counter[0] += 1
        id_map[id(node)] = this_id

        # Recurse
        for sub in node['substituents']:
            traverse_assign_ids(sub['child'])

    traverse_assign_ids(parsed_node)
    total_atoms = node_id_counter[0]

    # Create placeholder atoms
    for _ in range(total_atoms):
        rwmol.AddAtom(Chem.Atom(0))  # atomicNum=0 => dummy

    # Helper for real atoms
    def get_atom_info(atom_label, charge_str):
        """
        Convert e.g. 'N', charge_str='-1' => (atomic_num=7, formal_charge=-1).
        If unknown, default to carbon (6).
        If 'lp', we keep it dummy => (0, 0).
        """
        if atom_label.lower() == 'lp':
            # We'll keep atomic_num=0 => dummy
            return 0, 0
        # minimal map
        symbol_map = {
            'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
            'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53
        }
        atomic_num = symbol_map.get(atom_label, 6)  # default carbon
        formal_charge = 0
        if charge_str is not None:
            formal_charge = int(charge_str)
        return atomic_num, formal_charge

    def add_atoms_and_bonds(node, parent_node=None, bond_order=None):
        this_idx = id_map[id(node)]
        atom_label = node['atom']
        atom_charge = node['charge']

        # Set up the placeholder atom
        atomic_num, fcharge = get_atom_info(atom_label, atom_charge)
        atom = rwmol.GetAtomWithIdx(this_idx)
        atom.SetAtomicNum(atomic_num)
        atom.SetFormalCharge(fcharge)

        # If it's an lp dummy, attach a property so we can rename it
        if atom_label.lower() == 'lp':
            # We can store a special property so the RDKit drawer
            # will show "lp" instead of "*"
            # "molFileAlias" is recognized by RDKit if we skip normal sanitize
            atom.SetProp("molFileAlias", "lp")

        if parent_node is not None:
            parent_idx = id_map[id(parent_node)]
            rd_bond = Chem.BondType.SINGLE
            if bond_order == '=':
                rd_bond = Chem.BondType.DOUBLE
            elif bond_order == '#':
                rd_bond = Chem.BondType.TRIPLE

            # Check if there's already a bond
            if rwmol.GetBondBetweenAtoms(parent_idx, this_idx) is None:
                rwmol.AddBond(parent_idx, this_idx, rd_bond)

        # Recurse for children
        for sub in node['substituents']:
            add_atoms_and_bonds(sub['child'], node, sub['bond_order'])

    add_atoms_and_bonds(parsed_node, None, None)

    # Convert to normal Mol, skipping standard sanitization
    mol = rwmol.GetMol()
    # We skip full sanitization because dummy atoms with real bonds can fail valence checks
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
    # This partial sanitize ensures we at least get ring info, etc., if needed
    Chem.rdmolops.SetAromaticity(mol)

    return mol

###############################################################################
# 3) DEMO: PARSE, BUILD, DRAW
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

    # 2) Build RDKit molecule (with "lp" as dummy atoms)
    mol = build_rdkit_mol(tree)

    # 3) Compute 2D coords
    # RDKit will lay out the entire graph, including dummy "lp" atoms.
    AllChem.Compute2DCoords(mol)

    # 4) Draw to an image
    img = Draw.MolToImage(mol, size=(400, 400))
    img.show()  # pop-up if environment supports it

    # Optionally save
    # img.save("diethylamine_with_lp.png")

if __name__ == "__main__":
    main()
