#!/usr/bin/env python3

import os
import argparse
import pandas as pd


def load_vaf(filepath):
    df = pd.read_csv(filepath, sep="\t", header=None, dtype=str)

    # build variant key
    df["key"] = df[1] + ":" + df[2] + ":" + df[3] + ":" + df[4] + ":" + df[5]

    # VAF column (last column = index 9)
    df["VAF"] = pd.to_numeric(df[9], errors="coerce")

    return df[["key", "VAF"]].dropna()


def load_all_samples(input_dir):
    sample_data = {}

    for sample in sorted(os.listdir(input_dir)):
        sample_dir = os.path.join(input_dir, sample)
        if not os.path.isdir(sample_dir):
            continue

        t_file = os.path.join(sample_dir, f"{sample}.T.glrc.vaf")
        n_file = os.path.join(sample_dir, f"{sample}.N.glrc.vaf")

        if os.path.exists(t_file):
            print(f"Loading {sample}.T")
            sample_data[f"{sample}.T"] = load_vaf(t_file)

        if os.path.exists(n_file):
            print(f"Loading {sample}.N")
            sample_data[f"{sample}.N"] = load_vaf(n_file)

    return sample_data


def compute_pairwise_corr(sample_data):
    names = list(sample_data.keys())
    corr_matrix = pd.DataFrame(index=names, columns=names, dtype=float)

    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            if j < i:
                continue

            print(f"Computing: {name1} vs {name2}")

            df1 = sample_data[name1].rename(columns={"VAF": "v1"})
            df2 = sample_data[name2].rename(columns={"VAF": "v2"})

            merged = df1.merge(df2, on="key", how="inner")

            if merged.shape[0] < 2:
                corr = None
            else:
                corr = merged["v1"].corr(merged["v2"])

            corr_matrix.loc[name1, name2] = corr
            corr_matrix.loc[name2, name1] = corr  # symmetric

    return corr_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    sample_data = load_all_samples(args.input_dir)

    print("\nCalculating pairwise correlations...\n")
    corr_matrix = compute_pairwise_corr(sample_data)

    corr_matrix.index.name = "Name"
    corr_matrix.to_csv(args.output, sep="\t")

    print(f"\nDone. Output written to: {args.output}")


if __name__ == "__main__":
    main()