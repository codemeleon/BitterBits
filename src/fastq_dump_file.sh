#!/usr/bin/env bash

while read p; do
	fastq-dump --split-3 --gzip $p
done
