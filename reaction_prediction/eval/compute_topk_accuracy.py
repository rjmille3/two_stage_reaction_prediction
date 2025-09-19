#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import numpy as np
import rdkit
from rdkit import Chem


from rpCHEM.Common.Util import smi_to_unique_smi_fast, smi_to_unique_smi_map, exact_mass  # noqa: F401

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else ''

def load_targets(path, skip_header=False):
    """Load target reaction strings (first column per row). Returns a list[str]."""
    targets = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for idx, row in enumerate(reader):
            if skip_header and idx == 0:
                continue
            if not row:
                # keep alignment: treat empty row as empty string
                targets.append("")
                continue
            targets.append(row[0].strip())
    return targets

def load_preds_csv(path):
    """Load predictions CSV. Returns a list[list[str]] (empty list for blank lines)."""
    rows = []
    with open(path, "r", encoding="utf-8") as pf:
        reader = csv.reader(pf)
        for row in reader:
            # row is [] for a blank line
            rows.append(row)
    return rows

def compute_topk(targets, preds_rows, kmax=10, verbose=False):
    """
    targets: list[str] (reaction lines)
    preds_rows: list[list[str]] (each row: prediction strings; [] means blank)
    Returns: top_k_counts (len kmax), total_used
    """
    n_t = len(targets)
    n_p = len(preds_rows)
    n = min(n_t, n_p)
    if n_t != n_p:
        print(f"[WARN] Length mismatch (targets={n_t}, preds={n_p}). Using first {n} rows.")

    top_k_counts = [0] * kmax

    for i in range(n):
        tgt_line = targets[i]
        # If target row is blank, count as miss (keeps alignment)
        if not tgt_line:
            continue

        # The reaction string may have arrow codes after a space: "<rxnSmiles> <arrowCodes>"
        # Keep only reaction SMILES before the first space:
        rxn_smiles = tgt_line.split(" ")[0].strip()

        # Extract and canonicalize target products
        try:
            tgt_prods = rxn_smiles.split(">>")[1].strip()
        except IndexError:
            # ill-formed line → count as miss
            if verbose:
                print(f"[WARN] Bad target line @ {i}: {tgt_line!r}")
            continue
        tgt_prods = smi_to_unique_smi_fast(tgt_prods)

        row = preds_rows[i]
        if not row:
            # blank predictions line → miss
            continue

        # Normalize predictions: each cell looks like "<rxnSmiles> <arrowCodes>"
        prods_list = []
        for cell in row:
            cell = cell.strip().strip('"').strip("'")
            if not cell:
                continue
            try:
                pred_rxn_smiles = cell.split(" ")[0].strip()
                pred_prods = pred_rxn_smiles.split(">>")[1].strip()
                prods_list.append(smi_to_unique_smi_fast(pred_prods))
            except Exception:
                # Skip malformed prediction cell
                continue

        if not prods_list:
            # No usable predictions → miss
            continue

        # If target products appear at rank r, credit all top-k where k-1 >= r
        try:
            rank = prods_list.index(tgt_prods)  # 0-based
            for kk in range(rank, kmax):
                top_k_counts[kk] += 1
            if verbose:
                print(f"[HIT] idx={i} rank={rank+1} tgt={tgt_prods}")
        except ValueError:
            if verbose:
                print(f"[MISS] idx={i} tgt={tgt_prods}")

    return top_k_counts, n

def parse_args():
    p = argparse.ArgumentParser(description="Compute Top-K accuracy from targets and prediction CSV.")
    p.add_argument("--targets", required=True,
                   help="Path to target reactions file (CSV/TXT; reaction in first column).")
    p.add_argument("--preds", required=True,
                   help="Path to predictions CSV (each row: up to K prediction strings; blank line = miss).")
    p.add_argument("--kmax", type=int, default=10, help="Compute Top-K up to this K (default: 10).")
    p.add_argument("--skip_header", action="store_true",
                   help="Skip the first line of the targets file (if it's a header).")
    p.add_argument("--verbose", action="store_true", help="Print per-row hits/misses.")
    return p.parse_args()

def main():
    args = parse_args()

    targets = load_targets(args.targets, skip_header=args.skip_header)
    preds_rows = load_preds_csv(args.preds)

    top_k_counts, total = compute_topk(targets, preds_rows, kmax=args.kmax, verbose=args.verbose)

    # Summary
    for k in range(args.kmax):
        acc = top_k_counts[k] / total if total else 0.0
        print(f"top-{k+1} accuracy: {acc:.6f}")

    print(f"\nTotal evaluated rows: {total}")

if __name__ == "__main__":
    main()
