import click
from Bio import SeqIO
from os import path, makedirs, getcwd, symlink, system
from shutil import rmtree
from subprocess import call
import pandas as pd
from glob import glob
from multiprocessing import Pool, cpu_count
import os
import subprocess

#todo Arrange the genome contigs for minimum clash

def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull,
                         stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True


def blat(blat_paratmeters):
    pair, blatpath, outfold = blat_paratmeters

    sequences = {}
    for rec in SeqIO.parse(pair[1], 'fasta'):
        sequences[rec.id] = rec.seq

    outpsl = "tmp_psl/"+pair[1].split("/")[-1].split(".")[0]+'.psl'

    call([blatpath, "-noHead", pair[0], pair[1], outpsl])
    try:
        data = pd.read_table(outpsl, header=None)
    except:
        return
    data = data[data[0] > (0.3 * data[10])]
    data = data.sort_values(by=[9, 0], ascending=[1, 0])
    already_selected = []
    for i, row in data.iterrows():
        if row[9] in already_selected:
            continue
        already_selected.append(row[9])
        if row[8] == '-':
            sequences[row[9]] = sequences[row[9]].reverse_complement()

    with open("%s/%s" % (outfold, pair[1].split('/')[-1]), "a") as fout:
        for id_ in set(data[9]):
            fout.write(">%s\n%s\n" % (id_, sequences[id_]))
    with open(pair[1], 'w') as fout:
        for id_ in set(sequences.keys()) - set(data[9]):
            fout.write(">%s\n%s\n" % (id_, sequences[id_]))
    return

@click.command( )
@click.option("--confold", help="Contig folder path", type=str, default="fas")
@click.option("--outfold", help="Output folder path", type=str,
              default='./output')
@click.option("--ncor", help="Number of Processors", type=int, default=1)
def run(confold, outfold, ncor):
    assert path.isdir(confold), "Input folder %s desn't exist" % confold
    assert not path.isdir(outfold), ("Folder %s already exists. Deleting"
                                     % outfold)
    rmtree(outfold)

    if ncor > cpu_count():
        click.echo("Using %d core" % cpu_count())
        ncor = cpu_count()
    pool = Pool(ncor)

    if path.isdir("tmp_con"):
        rmtree("tmp_con")
    makedirs("tmp_con")

    if path.isdir("tmp_psl"):
        rmtree("tmp_psl")
    makedirs("tmp_psl")

    makedirs(outfold)
    for source in glob("%s/*" % confold):
        symlink("%s/%s" % (getcwd(), source),
                "%s/tmp_con/%s" % (getcwd(), source.split("/")[-1]))

    files = glob("tmp_con/*")
    # order the files in increase number of contigs for higher confidence
    # you can use pandas to do that
    file_contigg_count = {'files': [], 'contig_count': []}
    for fl in files:
        count = 0
        for rec in SeqIO.parse(fl, 'fasta'):
            count += 1
        file_contigg_count['files'].append(fl)
        file_contigg_count['contig_count'].append(count)
    dt = pd.DataFrame.from_dict(file_contigg_count)
    dt = dt.sort_values(by=['contig_count'], ascending=[True])
    files = list(dt['files'])
    del file_contigg_count, dt

    while len(files) > 1:
        for i, f1 in enumerate(files):
            pairs = []
            for j, f2 in enumerate(files):
                if i >= j:
                    continue
                pairs.append([f1, f2])
            pool.map(blat,
                     [[pair, blatpath, outfold] for pair in pairs])

            with open("%s/%s" % (outfold, f1.split('/')[-1]), "a") as fout:
                fout.write(open(f1).read())
            system("rm %s" % f1)
            files = files[1:]
            print(files)
            if len(files) == 1:
                with open("%s/%s" % (outfold, files[0].split('/')[-1]),
                          "a") as fout:
                    fout.write(open(files[0]).read())
                break
            if len(pairs) == 1:
                break

    rmtree('tmp_psl')
    rmtree('tmp_con')


if __name__ == '__main__':
    run()
