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
nmovies = sum([len(m) for m in movies.values()])
moviejobs = gcd(maxthreads, nmovies)
moviethreads = maxthreads//moviejobs

# create the log directory if needed
if not os.path.exists("logs"):
    os.mkdir("logs")

rule all:
    input:
        expand("process/{seqrun}/trimmed/{seqrun}.fastq.gz", seqrun = rs2_runs),
        "_targets.R",
        glob("scripts/*.R"),
        makescript = "make.R"
        planscript = "_targets.R"
    output: touch(".finished")
    conda: "conda/soil-fungi.yaml"
    script: "make.R"

# convert a raw RSII-format (.h5) movie to the Sequel format (.bam)
# these files are pretty large, so they are marked as temporary.
rule bax2bam:
    output:
        temp("process/{seqrun}/rawmovie/{movie}.subreads.bam"),
        temp("process/{seqrun}/rawmovie/{movie}.subreads.bam.pbi"),
        temp("process/{seqrun}/rawmovie/{movie}.scraps.bam"),
        temp("process/{seqrun}/rawmovie/{movie}.scraps.bam.pbi")
    input:
        lambda wildcards: glob(f"rawdata/{wildcards.seqrun}/**/{wildcards.movie}.*.bax.h5", recursive = True)
    shadow: "shallow"
    params:
        prefix="process/{seqrun}/rawmovie/{movie}",
        dir = "process/{seqrun}/rawmovie/"
    resources:
        walltime=20
    threads: 2
    log: "logs/bax2bam_{seqrun}_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1" # no bax2bam in newer versions
    shell:
        """
        [ -d {params.dir} ] || mkdir -p {params.dir}
        bax2bam {input} -o {params.prefix} &> {log}
        """

# endpoint target: convert all pacbio movies to Sequel format
rule convertmovies:
    input:
        [expand("process/{seqrun}/rawmovie/{movie}.{type}.bam",
               seqrun = s,
               movie = m,
               type = ['subreads', 'scraps'])
         for (s,m) in movies.items()]

# merge all the movies which belong to the same plate
rule mergebam:
    output: temp("process/{seqrun}/rawmovie/{seqrun}.subreads.bam")
    input:
        lambda wildcards: expand("process/{{seqrun}}/rawmovie/{m}.subreads.bam", m = movies[wildcards.seqrun])
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
    output:
        bam = "process/{seqrun}/ccs/{movie}.bam",
        index = "process/{seqrun}/ccs/{movie}.bam.pbi"
    input: "process/{seqrun}/rawmovie/{movie}.subreads.bam"
    params: dir = "process/{seqrun}/ccs"
    resources:
        walltime=120
    shadow: "shallow"
    threads: moviethreads
    log: "logs/ccs_{seqrun}_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1" # ccs from newer versions doesn't accept RSII data
    shell:
        """
        [ -d {params.dir} ] || mkdir -p {params.dir}
        ccs --numThreads {threads} --richQVs {input} {output.bam} &>{log}
        """

# convert a ccs BAM to a fastq
# this loses a lot of PacBio-specific information, but it is useful for other software.
rule bam2fastq:
    output: temp("process/{seqrun}/ccs/{movie}.fastq.gz")
    input:
        bam = "process/{seqrun}/ccs/{movie}.bam",
        pbi = "process/{seqrun}/ccs/{movie}.bam.pbi"
    params:
        basename = "process/{seqrun}/ccs/{movie}"
    resources:
             walltime=10
    threads: 1
    log: "logs/bam2fastq_{seqrun}_{movie}.ccs.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "bam2fastq -o {params.basename} {input.bam} &>{log}"

# lima doesn't store any information about orientation when run on subreads,
# so orient using the primers
rule orient:
    output:
        orient = "process/{seqrun}/orient/{movie}.fastq.gz",
        noprimer = "process/{seqrun}/lost/noprimer/{movie}.fastq.gz"
    input:
        ccs = "process/{seqrun}/ccs/{movie}.fastq.gz"
    params:
        orientdir = "process/{seqrun}/orient",
        noprimerdir = "process/{seqrun}/lost/noprimer"
    threads: moviethreads
    log: "logs/orient_{seqrun}_{movie}.log"
    conda: "conda/orient.yaml"
    shell:
        """
        [ -d {params.orientdir} ] || mkdir -p {params.orientdir}
        [ -d {params.noprimerdir} ] || mkdir -p {params.noprimerdir}
        cutadapt\\
            -a "TCCGTAGGTGAACCTGC;e=0.15...CGAAGTTTCCCTCAGGA;required;e=0.15"\\
            --action=none\\
            --revcomp\\
            -o {output.orient}\\
            --untrimmed-output {output.noprimer}\\
            -j {threads}\\
            {input.ccs} &>{log}
        """

# quality filter the ccs and dereplicate
# allow up to 1% expected errors, minimum length 1000, maximum length 2000
# the fasta output has only one entry for each unique sequence
# the label is the ZMW name of one of the first appearance, followed by ";size=n"
# "n" gives the number of times the sequence appears.
rule derep:
    output:
        derep="process/{seqrun}/derep/{seqrun}.fasta",
        orient="process/{seqrun}/orient/{seqrun}.fastq.gz",
        trimmed="process/{seqrun}/trimmed/{seqrun}.fastq.gz",
        filtered="process/{seqrun}/filtered/{seqrun}.fastq.gz",
        tooshort="process/{seqrun}/lost/tooshort/{seqrun}.fastq.gz",
        toolong="process/{seqrun}/lost/toolong/{seqrun}.fastq.gz",
        toopoor="process/{seqrun}/lost/toopoor/{seqrun}.fastq.gz",
        uc="process/{seqrun}/derep/{seqrun}.uc"
    input:
        lambda wildcards: expand("process/{seqrun}/ccs/{movie}.fastq.gz",
                                 seqrun = wildcards.seqrun,
                                 movie = movies[wildcards.seqrun])
    params:
        derepdir = "process/{seqrun}/derep/",
        trimmeddir = "process/{seqrun}/trimmed/",
        filtereddir = "process/{seqrun}/filtered/",
        tooshortdir= "process/{seqrun}/lost/tooshort/",
        toolongdir= "process/{seqrun}/lost/toolong/",
        toopoordir= "process/{seqrun}/lost/toopoor/"
    resources:
        walltime=10
    shadow: "shallow"
    threads: 2
    log:
        derep = "logs/{seqrun}/derep.log",
        trim = "logs/{seqrun}/trim.log",
        qualfilter = "logs/{seqrun}/qualfilter.log",
        maxlenfilter = "logs/{seqrun}/maxlenfilter.log",
        minlenfilter = "logs/{seqrun}/minlenfilter.log"
    conda: "conda/orient.yaml"
    shell:
        """
         [ -d {params.derepdir} ] || mkdir -p {params.derepdir}
         [ -d {params.trimmeddir} ] || mkdir -p {params.trimmeddir}
         [ -d {params.filtereddir} ] || mkdir -p {params.filtereddir}
         [ -d {params.tooshortdir} ] || mkdir -p {params.tooshortdir}
         [ -d {params.toolongdir} ] || mkdir -p {params.toolongdir}
         [ -d {params.toopoordir} ] || mkdir -p {params.toopoordir}
         trimmed=$(mktemp --suffix .fastq) &&
         filtered=$(mktemp --suffix .fastq) &&
         tooshort=$(mktemp --suffix .fastq) &&
         toolong=$(mktemp --suffix .fastq) &&
         toopoor=$(mktemp --suffix .fastq) &&
         trap 'rm ${{trimmed}} ${{filtered}} ${{tooshort}} ${{toolong}} ${{toopoor}}' EXIT &&
         cat {input} > {output.orient} &&
         cutadapt \\
            -a "TCCGTAGGTGAACCTGC;e=0.15...CGAAGTTTCCCTCAGGA;required;e=0.15"\\
            -o -\\
            -j 1\\
            {output.orient}\\
            2>{log.trim} |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_maxee_rate 0.01 \\
            --fastq_qmax 93 \\
            --fastqout_discarded ${{toopoor}} \\
            --fastqout - \\
            2>{log.qualfilter} |
         tee ${{trimmed}} |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_qmax 93 \\
            --fastq_minlen 1000 \\
            --fastqout_discarded ${{tooshort}} \\
            --fastqout - \\
            2>{log.minlenfilter} |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_qmax 93 \\
            --fastq_maxlen 2000 \\
            --fastqout_discarded ${{toolong}} \\
            --fastqout ${{filtered}}\\
            --fastaout - \\
            2>{log.maxlenfilter} |
         vsearch --derep_fulllength - \\
            --threads 1 \\
            --sizeout \\
            --fasta_width 0\\
            --output {output.derep} \\
            --uc {output.uc} \\
            &>{log.derep} &&
         gzip -c ${{trimmed}} >{output.trimmed} &&
         gzip -c ${{filtered}} >{output.filtered} &&
         gzip -c ${{tooshort}} >{output.tooshort} &&
         gzip -c ${{toolong}} >{output.toolong} &&
         gzip -c ${{toopoor}} >{output.toopoor}
        """
