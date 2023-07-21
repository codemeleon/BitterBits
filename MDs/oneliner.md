### Counting reads in `fastq.gz` files in a folder.

- `parallel "echo -n '{},' && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: *.fastq.gz`

## bamstats

- `parallel "echo -n '{},' && bam_stat.py  -i {} > {}.bam_stat" ::: *.bam`
- `parallel "clipping_profile.py -i {} -s 'PE' -o  cpo{}" ::: *.bam`
- `parallel "deletion_profile.py -i {} -l 101 -o dp{}" ::: *.bam`
- `parallel "deletion_profile.py -i {} -l 101 -o dp{}" ::: *.bam`
  > > `bam_stat.py`, `clipping_profile.py`, `deletion_profile.py` are part of `rseqc`

ribodetector_cpu -l 150 -i ../trimmed/N1_R1_001.fastq.gz ../trimmed/N1_R2_001.fastq.gz -o N1_R1_001.fastq.gz N1_R2_001.fastq.gz -r N1_R1_rrna.fastq.gz N1_R2_rrna.fastq.gz -t 5 --chunk_size 2000 && telegram-send 'rRNA removal done'
