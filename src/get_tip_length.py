#!/usr/bin/env python3

from ete3 import Tree
import argparse


def main():
    tree = Tree(treefile)
    # print the leaf_name and terminal_branch_length to csv
    print("leaf_name,terminal_branch_length")
    for leaf in tree:
        print(f"{leaf.name},{leaf.dist}")


# dev
# iqtree_output_tarfile = Path(
#     "output/090_rerun-full-dataset/030_iqtree/NUC.NT.muscle_align.gap0.8.cov0.8.wscore0.4/iqtree_files.tar"
# )

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract terminal branch lengths from a Newick tree in a tarfile."
    )
    parser.add_argument(
        "treefile", type=str, help="Path to the tarfile containing the Newick tree."
    )

    args = parser.parse_args()
    treefile = args.treefile

    main()
