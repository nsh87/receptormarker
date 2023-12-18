"""Generate a phylo.xml from a MUSCLE MSA.fasta

For usage use 'python phyloxml_from_msa.py -h'
"""


import argparse


##### PARSE ARGUMENTS #####
argparser = argparse.ArgumentParser()
argparser.add_argument("--msa", dest="msa", required=True, metavar="MSA_FILE",
                       help="path to a MUSCLE_msa.fasta to create phylo from")
argparser.add_argument("--dest", dest="dest", required=True,
                       help="path to an output phylo.xml")
argparser.add_argument("--distance", dest="distance", required=True,
                       help="distance calculator to use on MSA")
args = argparser.parse_args()
##### END PARSE ARGUMENTS #####


def phyloxml_from_msa(msa, phyloxml, distance):
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    from Bio import Phylo
    ms_alignment = AlignIO.parse(msa, "fasta")
    # For some reason, we're getting weird chars appending to the IDs now
    for record in ms_alignment:
        if '.' in record.id:
            record.id = record.id.split('.')[0]
    calculator = DistanceCalculator(distance)
    dist_matrix = calculator.get_distance(ms_alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dist_matrix)
    Phylo.write(tree, phyloxml, "phyloxml")


if __name__ == "__main__":
    msa = args.msa
    phyloxml = args.dest
    distance = args.distance
    phyloxml_from_msa(msa, phyloxml, distance)

