#!/usr/bin/env python3

import os
import argparse
import math
import csv


def load_vaf(filepath):
    data = {}

    with open(filepath) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            try:
                key = f"{row[1]}:{row[2]}:{row[3]}:{row[4]}:{row[5]}"
                vaf = float(row[9])
                data[key] = vaf
            except:
                continue

    return data


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


def pearson_corr(x, y):
    n = len(x)
    if n < 2:
        return None

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    den_x = sum((xi - mean_x) ** 2 for xi in x)
    den_y = sum((yi - mean_y) ** 2 for yi in y)

    if den_x == 0 or den_y == 0:
        return None

    return num / math.sqrt(den_x * den_y)


def compute_pairwise_corr(sample_data):
    names = list(sample_data.keys())
    corr_matrix = {n: {} for n in names}

    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            if j < i:
                continue

            print(f"Computing: {name1} vs {name2}")

            data1 = sample_data[name1]
            data2 = sample_data[name2]

            shared_keys = set(data1.keys()) & set(data2.keys())

            if len(shared_keys) < 2:
                corr = None
            else:
                v1 = [data1[k] for k in shared_keys]
                v2 = [data2[k] for k in shared_keys]
                corr = pearson_corr(v1, v2)

            corr_matrix[name1][name2] = corr
            corr_matrix[name2][name1] = corr

    return names, corr_matrix


def write_matrix(names, matrix, outfile):
    with open(outfile, "w") as f:
        f.write("Name\t" + "\t".join(names) + "\n")
        for n1 in names:
            row = [n1]
            for n2 in names:
                val = matrix[n1].get(n2)
                row.append("" if val is None else str(val))
            f.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    sample_data = load_all_samples(args.input_dir)

    print("\nCalculating pairwise correlations...\n")
    names, corr_matrix = compute_pairwise_corr(sample_data)

    write_matrix(names, corr_matrix, args.output)

    print(f"\nDone. Output written to: {args.output}")


if __name__ == "__main__":
    main()