import click
from multiprocessing import Pool
from subprocess import PIPE, Popen
from os import system, path
from glob import glob
import pandas as pd
import pysam as psm
from tempfile import mktemp

@click.command()
@click.option("-ref", help="CPS locus reference sequences", type=str,
              default=None, show_default=True)
@click.option("-fqf", help="Fastq files folder", type=str, default=None,
              show_default=True)
@click.option("-rt", help="Read type", type=click.Choice(["paired", "single"]),
              default='paired', show_default=True)
@click.option("-fw", help="Forward reads file symbol", type=str, default='_1',
              show_default=True)
@click.option("-rv", help="Reverse reads file symbol", type=str, default='_2',
              show_default=True)
@click.option("-fext", help="File extention", type=str, default='fastq.gz',
              show_default=True)
@click.option("-ncor", help="Number of CPU used", type=int, default=1,
              show_default=True)
@click.option("-ofile", help="Output csv file", type=str,
              default="Serotypes.csv", show_default=True)
@click.option("-best", help="Best hit based. Else All", type=bool,
              default=True, show_default=True)
## Add other option to generate differet kinds of file
def run(ref, fqf, rt, fw, rv, fext, ncor, ofile, best):
    '''in silico serotyping using short reads and CPS reference.'''
    # Generate index
    # Expect gzipped files
    # process = Popen(["bowtie2-build", "--large-index", ref, ref],
    #                 stdout=PIPE, stderr=PIPE)
    # stdout, stderr = process.communicate()
    dfs = []
    tmp_file = mktemp()#dir="/tmp"
    for fl in glob("%s/*%s.%s" % (fqf, fw, fext)):
        pathbs = fl.split("%s.%s" % (fw, fext))[0]
        revpath = "%s%s.%s" % (pathbs, rv, fext)
        strain = path.split(fl)[1].split("%s.%s" % (fw, fext))[0] # Replace this with somthing better
        # print(fl, revpath)
        if not best:
            process = Popen(["bowtie2", "--end-to-end", "--very-sensitive",
                             "-a", "-N", "1", "-L", "4", "-x", ref, "-1", fl,
                              "-2", revpath,
                              '-p', str(ncor), "-S", "%s.sam" % tmp_file],
                            stdout=PIPE, stderr=PIPE)
        else:
            process = Popen(["bowtie2", "--end-to-end", "--very-sensitive",
                              "-N", "1", "-L", "4", "-x", ref, "-1", fl,
                              "-2", revpath,
                              '-p', str(ncor), "-S", "%s.sam" % tmp_file],
                            stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        process = Popen(["samtools", "view","-bS", "-o", "%s.bam" % tmp_file,
                         "%s.sam" % tmp_file], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        process = Popen(["samtools", "sort", "%s.bam" % tmp_file, #"-o",
                         "%s_sorted" % tmp_file], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        process = Popen(["samtools", "index", "%s_sorted.bam" % tmp_file],
                        stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()

        sam = psm.Samfile("%s_sorted.bam" % tmp_file, "rb")

        ref_sizes = {}
        for refr, sz in zip(sam.references, sam.lengths):
            ref_sizes[refr] = sz


        read_count = {'ref':[], 'depth':[]}
        idtt = 1.
        for refr in sam.references:
            # print(ref)
            rcount = 0
            for read in sam.fetch(refr):
                if read.alen == None:
                    continue
                if read.alen >= read.qlen * idtt:
                    rcount += 1
            read_count['ref'].append(refr)
            read_count['depth'].append(rcount * 1./ref_sizes[refr])
        process = Popen(["rm", "-rf",
                         "%s*" % tmp_file], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        rd_count = pd.DataFrame.from_dict(read_count)
        rd_count["samp"] = strain
        # rd_count.sort_values("depth", ascending=False)
        df = rd_count.pivot_table(index="samp", columns="ref",
                                  values="depth").reset_index()
        dfs.append(df)
    process = Popen(["rm", "-rf", "%s*" % tmp_file],
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    dfs = pd.concat(dfs)
    columns = list(dfs.columns)
    columns.remove("samp")
    dfs['max'] = dfs[columns].idxmax(axis=1)
    dfs = dfs[["samp", "max"] + columns]
    dfs.to_csv(ofile, index=False, sep=",")


if __name__ == '__main__':
    run()
