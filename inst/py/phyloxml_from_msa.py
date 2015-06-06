"""Generate a phylo.xml from a MUSCLE MSA.fasta"""


import argparse


##### PARSE ARGUMENTS #####
argparser = argparse.ArgumentParser()
argparser.add_argument("msa", help="path to the MUSCLE MSA.fasta")
argparser.add_argument("phyloxml", help="path to an output phylo.xml")
args = argparser.parse_args()
args = vars(args)  # Give access to arguments using a dict: e.g. args["msa"]
##### END PARSE ARGUMENTS #####


def phyloxml_from_msa(msa, phyloxml):
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    from Bio import Phylo
    ms_alignment = AlignIO.read(msa, "fasta")
    calculator = DistanceCalculator("ident")
    dist_matrix = calculator.get_distance(ms_alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dist_matrix)
    Phylo.write(tree, phyloxml, "phyloxml")


if __name__ == "__main__":
    msa = args["msa"]
    phyloxml = args["phyloxml"]
    phyloxml_from_msa(msa, phyloxml)

