#!/usr/bin/env python

import subprocess
import time
from subprocess import DEVNULL

with open(str(snakemake.log), "w") as l:
    s = time.time()
    l.write("*** Getting command ***\n")

    command_client = [
        "diamond-aligner",
        "blastp",
        "-d",
        str(snakemake.input.db),
        "-q",
        str(snakemake.input.fasta),
        "-o",
        str(snakemake.output),
        "-e",
        "1e-5",
        "--sensitive",
        "-f",
        "6",
        "qseqid",
        "qlen",
        "qstart",
        "qend",
        "sseqid",
        "slen",
        "sstart",
        "send",
        "length",
        "pident",
        "ppos",
        "score",
        "evalue",
        "bitscore",
    ]

    command_cluster = [
        "diamond",
        "blastp",
        "-d",
        str(snakemake.input.db),
        "-q",
        str(snakemake.input.fasta),
        "-o",
        str(snakemake.output),
        "-e",
        "1e-5",
        "--sensitive",
        "-f",
        "6",
        "qseqid",
        "qlen",
        "qstart",
        "qend",
        "sseqid",
        "slen",
        "sstart",
        "send",
        "length",
        "pident",
        "ppos",
        "score",
        "evalue",
        "bitscore",
    ]
    l.write("*** Alignment of sequences ***\n")
    try:
        subprocess.call(command_client, stdout=DEVNULL)

    except:
        subprocess.call(command_cluster, stdout=DEVNULL)
    e = time.time()
    l.write(f"Alignment done in {round(e - s,2)} seconds")

