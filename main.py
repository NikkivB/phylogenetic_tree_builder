"""
Generating a phylogenetic tree with questionable measures.
The strains will be made equal by adding '-' at the end of them.
If you give contigs or reads, then they will be put after each other to create
one sequence.
Date:
    - 09 november 2026
Authors:
    - Nikki van Braam
    - Aldiana Kučević
    - Timo Klaassen
"""
from Bio import AlignIO, SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.PhyloXML import Phylogeny
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt


def combine_fasta(input_fastas: list[str], output_fasta: str) -> None:
    """
    Combines contigs/reads in a bit of a whacky way.
    Don't expect perfection. This is bruteforced.
    :param input_fastas: A list of the fasta files to combine.
    :param output_fasta: the name of the output fasta file.
    :return: None
    """
    i = 0
    records = []
    ids = ["ATCC_1015",
           "ATCC_13496",
           "CBS_101883",
           "N402",
           "SH2"]

    open(output_fasta, "w").close() # Erase data from file

    for file in input_fastas:

        with open(file, "r") as current_file:
            seqs = [str(record.seq) for record in
                    SeqIO.parse(current_file, "fasta")]
            combined_seq = Seq("".join(seqs))
            new_record = SeqRecord(combined_seq, id="",
                                   description="")
            records.append(new_record)

    # Make the sequences equal length
    max_len = max(len(r) for r in records)
    padded_sequences = [str(s.seq).ljust(max_len, "-") for s in
                        records]

    records_padded = []
    temp = 0
    for item in padded_sequences:
        records_padded.append(SeqRecord(Seq(item), id=str(ids[temp]), description=""))
        temp += 1

    SeqIO.write(records_padded, output_fasta, "fasta")


def align(input_file: str):
    """
    Create a MSA based on the given sequences with the Biopython package.
    :param input_file: The file in which the sequences are. All sequences must
    be in one file.
    :return: alignment, msa_results
    """
    with open("data/output/msa_results.phy", "w") as msa_results:
        alignment = AlignIO.read(input_file, "fasta")
        print("Alignment length %i" % len(alignment))
        print(alignment)

        AlignIO.write(alignment, "data/output/msa_results.phy", "phylip")

        for record in alignment:
            print(record.seq + " " + record.id)

        print(type(alignment))
        msa_results.close()

    return alignment, msa_results


def tree_generating(alignment):
    """

    :param msa_results:
    :return:
    """

    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator, "nj")
    tree = constructor.build_tree(alignment)
    print(tree)
    Phylo.draw_ascii(tree)
    plt.title("Phylogenetic tree of A.niger strains")
    plt.savefig("data/output/phylogenetic_tree.png")



if __name__ == "__main__":
    in_fasta: list[str] = [
        "data/1015/ncbi_dataset/data/GCA_000230395.2/GCA_000230395.2_ASPNI_v3.0_genomic.fna",
        "data/13496/ncbi_dataset/data/GCA_003344705.1/GCA_003344705.1_Aspni_bvT_1_genomic.fna",
        "data/cbs_101883/ncbi_dataset/data/GCA_003184595.1/GCA_003184595.1_Asplac1_genomic.fna",
        "data/N402/ncbi_dataset/data/GCA_900248155.1/GCA_900248155.1_Aniger_ATCC_64974_N402_genomic.fna",
        "data/sh2/ncbi_dataset/data/GCA_000633045.1/GCA_000633045.1_AspAngSh2_genomic.fna"
    ]

    intermediate = "data/output/combined_fasta_aspergillus_niger.fasta"

    combine_fasta(in_fasta, intermediate)
    alignment, msa_results = align(intermediate)
    tree_generating(alignment)