# Snakemake workflow file
# Brendan Furneaux
# April 2021

# This workflow identifies raw PacBio sequencing reads and performs command-line based
# operations like circular consensus calling, demultiplexing, and primer trimming.

import os.path
from glob import glob
import re
import subprocess
from math import gcd
from snakemake.io import glob_wildcards

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = int(subprocess.check_output("nproc").decode())

# find the PacBio RSII movie files
rs2_runs = [r.split(os.path.sep)[1] for r in glob("rawdata/pb_*")]
movies = {r:[os.path.basename(m)[:-7] for m in glob(f"rawdata/{r}/**/*.bas.h5", recursive = True)]
            for r in rs2_runs}


# figure out how to distribute the movie files equally over the cores we have available.
moviejobs = gcd(maxthreads, sum([len(m) for m in movies]))
moviethreads = maxthreads/moviejobs

# create the log directory if needed
if not os.path.exists("logs"):
    os.mkdir("logs")

rule all:
    input: expand("process/{seqrun}.ccs.trimmed.fastq.gz", seqrun = rs2_runs)

# convert a raw RSII-format (.h5) movie to the Sequel format (.bam)
# these files are pretty large, so they are marked as temporary.
rule bax2bam:
    output:
        temp("process/{seqrun}_{movie}.subreads.bam"),
        temp("process/{seqrun}_{movie}.subreads.bam.pbi"),
        temp("process/{seqrun}_{movie}.scraps.bam"),
        temp("process/{seqrun}_{movie}.scraps.bam.pbi")
    input:
        lambda wildcards: glob(f"rawdata/{wildcards.seqrun}/**/{wildcards.movie}.*.bax.h5", recursive = True)
    shadow: "shallow"
    params:
        prefix="process/{seqrun}_{movie}"
    resources:
        walltime=20
    threads: 2
    log: "logs/bax2bam_{seqrun}_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1" # no bax2bam in newer versions
    shell: "bax2bam {input} -o {params.prefix} &> {log}"

# endpoint target: convert all pacbio movies to Sequel format
rule convertmovies:
    input:
        [expand("process/{seqrun}_{movie}.{type}.bam",
               seqrun = s,
               movie = m,
               type = ['subreads', 'scraps'])
         for (s,m) in movies.items()]

# merge all the movies which belong to the same plate
rule mergebam:
    output: temp("process/{seqrun}.subreads.bam")
    input:
        lambda wildcards: expand("process/{m}.subreads.bam", m = movies[wildcards.seqrun])
    shadow: "shallow"
    resources:
        walltime=20
    log: "logs/mergebam_{seqrun}.log"
    conda: "conda/samtools.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools"
    shell: "samtools merge {output} {input} &>{log}"

# a seqrun identifier should start with pb (RSII) or ps (Sequel), then 3 digits,
# and optionally 3 more if multiple plates were sequenced
wildcard_constraints:
    seqrun = "p[bs]_\d{3}(_\d{3})?",
    movie = "m\d{6}_.*"

# generate circular consensus sequences from subreads
rule ccs:
    output: "process/{seqrun}_{movie}.ccs.bam"
    input: "process/{seqrun}_{movie}.subreads.bam"
    resources:
        walltime=120
    shadow: "shallow"
    threads: moviethreads
    log: "logs/ccs_{seqrun}_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1" # ccs from newer versions doesn't accept RSII data
    shell: "ccs --numThreads {threads} --richQVs {input} {output} &>{log}"

# convert a ccs BAM to a fastq
# this loses a lot of PacBio-specific information, but it is useful for other software.
rule bam2fastq:
    output: temp("process/{name}.fastq.gz")
    input:
        bam = "process/{name}.bam",
        pbi = "process/{name}.bam.pbi"
    resources:
             walltime=10
    threads: 1
    log: "logs/bam2fastq_{name}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "bam2fastq -o process/{wildcards.name} {input.bam} &>{log}"

# lima doesn't store any information about orientation when run on subreads,
# so orient using the primers
rule orient:
    output:
        orient = "process/{seqrun}_{movie}.ccs.orient.fastq.gz",
        noprimer = "process/{seqrun}_{movie}.ccs.noprimer.fastq.gz"
    input:
        ccs = "process/{seqrun}_{movie}.ccs.fastq.gz"
    threads: moviethreads
    log: "logs/orient_{seqrun}_{movie}.log"
    conda: "conda/orient.yaml"
    shell:
        """
        cutadapt\\
            -a "TCCGTAGGTGAACCTGC;e=0.15...CGAAGTTTCCCTCAGGA;required;e=0.15"\\
            --action=none\\
            --revcomp\\
            -o {output.orient}\\
            --untrimmed-output {output.noprimer}\\
            -j {threads}\\
            {input.ccs}
        """

# quality filter the ccs and dereplicate
# allow up to 1% expected errors, minimum length 1000, maximum length 2000
# the fasta output has only one entry for each unique sequence
# the label is the ZMW name of one of the first appearance, followed by ";size=n"
# "n" gives the number of times the sequence appears.
rule derep:
    output:
        fasta="process/{seqrun}.ccs.derep.fasta",
        fastq="process/{seqrun}.ccs.orient.fastq.gz",
        trimmed="process/{seqrun}.ccs.trimmed.fastq.gz",
        tooshort="process/{seqrun}.ccs.orient.tooshort.fastq.gz",
        toolong="process/{seqrun}.ccs.orient.toolong.fastq.gz",
        toopoor="process/{seqrun}.ccs.orient.toopoor.fastq.gz",
        uc="process/{seqrun}.ccs.derep.uc"
    input:
        lambda wildcards: expand("process/{seqrun}_{movie}.ccs.orient.fastq.gz",
                                 seqrun = wildcards.seqrun,
                                 movie = movies[wildcards.seqrun])
    resources:
        walltime=10
    shadow: "shallow"
    threads: 2
    log: "logs/derep_{seqrun}.log"
    conda: "conda/orient.yaml"
    shell:
        """
         fastq=$(mktemp --suffix .fastq) &&
         tooshort=$(mktemp --suffix .fastq) &&
         toolong=$(mktemp --suffix .fastq) &&
         toopoor=$(mktemp --suffix .fastq) &&
         trap 'rm ${{fastq}} ${{tooshort}} ${{toolong}} ${{toopoor}}' EXIT &&
         cat {input} > {output.fastq} &&
         cutadapt \\
            -a "TCCGTAGGTGAACCTGC;e=0.15...CGAAGTTTCCCTCAGGA;required;e=0.15"\\
            -o -\\
            -j 1\\
            {output.fastq} |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_maxee_rate 0.01 \\
            --fastq_qmax 93 \\
            --fastqout_discarded ${{toopoor}} \\
            --fastqout - |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_qmax 93 \\
            --fastq_minlen 1000 \\
            --fastqout_discarded ${{tooshort}} \\
            --fastqout - |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_qmax 93 \\
            --fastq_maxlen 2000 \\
            --fastqout_discarded ${{toolong}} \\
            --fastqout ${{fastq}}\\
            --fastaout - |
         vsearch --derep_fulllength - \\
            --threads 1 \\
            --sizeout \\
            --fasta_width 0\\
            --output {output.fasta} \\
            --uc {output.uc} &&
         gzip -c ${{fastq}} >{output.trimmed} &&
         gzip -c ${{tooshort}} >{output.tooshort} &&
         gzip -c ${{toolong}} >{output.toolong} &&
         gzip -c ${{toopoor}} >{output.toopoor}
        """