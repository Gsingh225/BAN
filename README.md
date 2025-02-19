# Bracket-and-Arrow Molecular Notation (BAN)

## Overview

**Bracket-and-Arrow Notation (BAN)** is a **text-based** molecular notation system designed for **human readability**. It allows representation of organic and inorganic molecules, including:

- **Atoms** (explicitly including hydrogens if desired).
- **Lone pairs** (optional, with positional orientation).
- **Bond orders** (single, double, triple).
- **Formal or literal charges** on atoms.
- **Relative orientation** (e.g., top, bottom, left, right).

BAN is **not** a validation tool—it allows users to write structures freely, even if chemically nonsensical.

## Syntax

### 1. Basic Atom Block

```
AtomSymbol (optionalCharge) [ listOfSubstituents ]
```

#### Atom Symbol
- Standard element symbols: `C`, `N`, `O`, `H`, `S`, etc.
- `lp` for lone pairs.
- Charge notation in parentheses: `N(+1)`, `O(-2)`, etc.

#### Substituents Block
- Enclosed in brackets `[ ... ]`.
- Each substituent is written as:
  ```
  orientation -> bondSpec
  ```
  or simply:
  ```
  orientation -> child
  ```
- **Orientation** can be `top`, `bottom`, `left`, `right`, or any custom label.
- **Bond specification**:
  - `-` = single bond (default if omitted)
  - `=` = double bond
  - `#` = triple bond
- **Child**: Either another **atom block** or a single atom symbol (`H`, `lp`, `Cl(-1)`, etc.).

## Examples

### 1. Methane (CH₄)

```
C [ -H, -H, -H, -H ]
```

or explicitly showing orientation:

```
C [
  top -> -H,
  left -> -H,
  right -> -H,
  bottom -> -H
]
```

### 2. Water (H₂O)
```
O [
  top -> lp,
  left -> -H,
  right -> -H
]
```

### 3. Ammonium Cation (NH₄⁺)
```
N(+1) [ -H, -H, -H, -H ]
```

### 4. Diethylamine (C₄H₁₁N)
```
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
```

## Quick Reference

| Feature         | Syntax Example |
|----------------|---------------|
| **Atom**       | `C`, `N(-1)`, `lp` |
| **Substituents List** | `[ ... ]` |
| **Substituent** | `Orientation -> BondOrder?AtomBlock` |
| **Orientation** | `top`, `bottom`, `left`, `right`, `customLabel` |
| **Bond Order** | `-` (single), `=` (double), `#` (triple) |
| **Charges** | `(+1)`, `(-2)`, `(0)` |
| **Lone Pairs** | `lp` |

### Example
```
C (0) [ top -> -H, bottom -> =O, left -> -H, right -> -OH ]
```

## Features & Design Goals
- **Human-readable** and **structured**.
- **Flexible**—allows multiple lone pairs, custom orientation labels.
- **No validation**—you can freely represent molecular structures.
- **Extensible** for advanced notation or custom labeling.

---

## License
This project is licensed under the MIT License.
