"""
    Nah,
"""
import Bio
from Bio import AlignIO, SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.PhyloXML import Phylogeny


def combine_fasta(input_fastas: list[str], output_fasta: str) -> None:
    """
    Combines contigs/reads in a bit of a whacky way.
    Don't expect perfection. This is bruteforced.
    :param input_fastas: A list of the fasta files to combine.
    :param output_fasta: the name of the output fasta file.
    :return: None
    """
    temp = ""
    i = 0
    records = []

    open(output_fasta, "w").close() # Erase data from file

    for file in input_fastas:
        with open(file, "r") as current_file:
            i += 1
            seqs = [str(record.seq) for record in
                    SeqIO.parse(current_file, "fasta")]
            combined_seq = Seq("".join(seqs))
            new_record = SeqRecord(combined_seq, id=file,
                                   description="")
            records.append(new_record)

            max_len = max(len(r) for r in records)
            padded_sequences = [str(s.seq).ljust(max_len, "-") for s in
                                records]
            records_padded = (SeqRecord(Seq(s)) for s in padded_sequences)

        SeqIO.write(records_padded, output_fasta, "fasta")


# todo: MSA uitvoeren
def align(input_file: str):
    """

    :param input_file:
    :return:
    """
    alignment = AlignIO.read(input_file, "fasta")
    msa_results = AlignIO.write(alignment, "data/output/msa_results.fasta", "fasta")
    print("Alignment length %i" % len(alignment))

    for record in alignment:
        print(record.seq + " " + record.id)

    print(type(alignment))



    return alignment, msa_results


# todo: MSA naar ander newwick

# todo: Fyloboom maken

def tree_generating(msa_results):

    tree = Phylo.read(msa_results, "newick")
    tree = Phylogeny.from_tree(tree)
    Phylo.draw(tree, do_show=False)

if __name__ == "__main__":
    input_fastas: list[str] = [
        "data/1015/ncbi_dataset/data/GCA_000230395.2/GCA_000230395.2_ASPNI_v3.0_genomic.fna",
        "data/13496/ncbi_dataset/data/GCA_003344705.1/GCA_003344705.1_Aspni_bvT_1_genomic.fna",
        "data/cbs_101883/ncbi_dataset/data/GCA_003184595.1/GCA_003184595.1_Asplac1_genomic.fna",
        "data/N402/ncbi_dataset/data/GCA_900248155.1/GCA_900248155.1_Aniger_ATCC_64974_N402_genomic.fna",
        "data/sh2/ncbi_dataset/data/GCA_000633045.1/GCA_000633045.1_AspAngSh2_genomic.fna"
    ]
    output_fasta = "data/output/combined_fasta_aspergillus_niger.fasta"

    combine_fasta(input_fastas, output_fasta)
    alignment, msa_results = align(output_fasta)
    tree_generating(msa_results)
