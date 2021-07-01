import click
import pandas as pd
import numpy as np
from Bio import SeqIO
import swifter


def most_common(lst):
    lst = [l for l in lst if l not in "-NRYKMSWBDHV"]
    return max(set(lst), key=lst.count)


@click.command()
@click.option("-aln",
              help="Alignment file",
              type=str,
              default="/home/devil/Downloads/All_seqs_aln.fasta",
              show_default=True)
@click.option(
    "-pheno",
    help="Phenotype File expecting names same as aln file",
    type=str,
    default=None,
    show_default=True)
def run(aln, pheno):
    """"""
    sequences = {}
    for rec in SeqIO.parse(aln, "fasta"):
        sequences[rec.id] = list(rec.seq.upper())
    sequences = pd.DataFrame(sequences)
    pheno = pd.read_csv(pheno)
    common_seq = list(set(pheno["samp"]) & set(sequences.columns))
    sequences = sequences[common_seq]
    pheno = pheno[pheno.samp.isin(common_seq)]
    pheno = pheno.set_index("samp")
    pheno = pheno.loc[common_seq, ]

    sequences = sequences[
        sequences.swifter.apply(
            lambda values: len(set(values.values) - set("-NRYKMSWBDHV")), axis=1
        )
        > 1
    ]
    # print(list(sequences.apply(lambda x: most_common(x.values), axis=1)))
    most_common_base = list(sequences.apply(
        lambda x: most_common(x.values), axis=1))
    all_bases = sequences.swifter.apply(
        lambda values: set(values.values) - set("-NRYKMSWBDHV"), axis=1
    )
    all_frames = []
    for idx, mc, all_b in zip(sequences.index, most_common_base, all_bases):
        values = np.array(sequences.loc[idx].values)
        remaining_bases = all_b - set([mc])
        for other_base in remaining_bases:
            init = [f"s_{idx}_{mc}_{other_base}", other_base, mc]
            tvalues = values.copy()
            tvalues[tvalues == mc] = 0
            tvalues[tvalues == other_base] = 1
            tvalues[~((tvalues == 0) | (tvalues == 1))] = "NA"
            all_frames.append(init+list(tvalues))
            # print(tvalues)

    # print(all_frames)
    all_frames = pd.DataFrame(
        all_frames, columns=["SNP", "minor", "major"]+common_seq)
    all_frames.to_csv("gemma.geno", index=False, header=False)
    pheno.to_csv("gemma.pheno", index=False, header=False)
    # print(all_frames)

    pass
#         print(idx, mc, all_b, remaining_bases)
#         pass
#
#     sequences.insert(0, "Y", list(sequences.apply(
#         lambda x: most_common(x.values), axis=1)))
#
#     sequences.insert(0, "SNP", [f"s_{pos}" for pos in sequences.index])
#
#     print(sequences.head())


if __name__ == '__main__':
    run()
