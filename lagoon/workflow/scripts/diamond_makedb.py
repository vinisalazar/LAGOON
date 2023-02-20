#!/usr/bin/env python

import subprocess
import time
from subprocess import DEVNULL

with open(str(snakemake.log), "w") as l:
    s = time.time()
    l.write("*** Getting command ***\n")
    command_client = [
        "diamond-aligner",
        "makedb",
        "--in",
        str(snakemake.input),
        "--db",
        str(snakemake.output),
    ]
    command_cluster = [
        "diamond",
        "makedb",
        "--in",
        str(snakemake.input),
        "--db",
        str(snakemake.output),
    ]
    l.write("*** Creation of Diamond database ***\n")
    try:
        subprocess.call(command_client, stdout=DEVNULL)

    except:
        subprocess.call(command_cluster, stdout=DEVNULL)
    e = time.time()
    l.write(f"Database created in {round(e - s,2)} seconds")