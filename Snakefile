import numpy as np
from Bio import SeqIO


THREADS = 32
FASTQ_IDS = [
  "LG1-016252-18horse", "LG2-053057-19horse",
  "LG3-012998-20horse", "LG5-002173-18horse"
]
DOWNSAMPLE = .2


rule all:
  input:
    expand("output/{run}/{reference}/quasirecomb/quasispecies.fasta", run=FASTQ_IDS, reference=['spike']),
    expand("output/{run}/{reference}/variants.vcf", run=FASTQ_IDS, reference=['spike']),
    expand("output/{run}/{reference}/qualimapReport.html", run=FASTQ_IDS, reference=['spike'])

rule downsample:
  input:
    reads="input/{run}_1.fq",
    mates="input/{run}_2.fq"
  output:
    reads="output/{run}/reads_1.fastq",
    mates="output/{run}/reads_2.fastq",
  run:
    reads = SeqIO.parse(input.reads, 'fastq')
    mates = SeqIO.parse(input.mates, 'fastq')
    total_all_reads = 0
    for read, mate in zip(reads, mates):
      assert read.name == mate.name, "Aborting due to mismatch in mate pairs"
      total_all_reads += 1
    print('Reads/mates are harmonized, proceeding to downsample')

    np.random.seed(1)
    reads = SeqIO.parse(input.reads, 'fastq')
    mates = SeqIO.parse(input.mates, 'fastq')
    total_kept_reads = np.ceil(DOWNSAMPLE*total_all_reads).astype(int)
    indices_to_keep = np.random.choice(total_all_reads, size=total_kept_reads, replace=False)
    indices_to_keep.sort()
    # append dummy index to the end to properly handle termination in loop below
    indices_to_keep = np.append(indices_to_keep, total_all_reads)

    keeping_index = 0
    fastq_handle = open(output.reads, 'w')
    for read_index, read in enumerate(reads):
      if read_index == indices_to_keep[keeping_index]:
        keeping_index += 1
        SeqIO.write(read, fastq_handle, 'fastq')

    keeping_index = 0
    fastq_handle = open(output.mates, 'w')
    for mate_index, mate in enumerate(mates):
      if mate_index == indices_to_keep[keeping_index]:
        keeping_index += 1
        SeqIO.write(mate, fastq_handle, 'fastq')

rule fastp:
  input:
    reads=rules.downsample.output.reads,
    mates=rules.downsample.output.mates
  output:
    reads="output/{run}/qc_1.fastq",
    mates="output/{run}/qc_2.fastq",
    html="output/{run}/fastp.html",
    json="output/{run}/fastp.json"
  shell:
    "fastp -i {input.reads} -I {input.mates} -o {output.reads} -O {output.mates} -h {output.html} -j {output.json} --thread %d" % THREADS

rule bowtie2_index:
  input:
    "input/{reference}.fasta"
  output:
    bt1="output/{reference}/{reference}.1.bt2",
    bt2="output/{reference}/{reference}.2.bt2",
    bt3="output/{reference}/{reference}.3.bt2",
    bt4="output/{reference}/{reference}.4.bt2",
    revbt1="output/{reference}/{reference}.rev.1.bt2",
    revbt2="output/{reference}/{reference}.rev.2.bt2"
  shell:
    "bowtie2-build {input} {wildcards.reference}/{wildcards.reference}"

rule reference_index:
  input:
    "input/{reference}.fasta"
  output:
    "input/{reference}.fasta.fai"
  shell:
    "lofreq faidx {input}"

rule bowtie2_alignment:
  input:
    rules.bowtie2_index.output,
    reads=rules.fastp.output.reads,
    mates=rules.fastp.output.mates
  params:
    index=lambda w:'output/%s/%s' % (w.reference, w.reference)
  output:
    temp("output/{run}/{reference}/mapped.sam")
  shell:
    "bowtie2 -x {params.index} -1 {input.reads} -2 {input.mates} -S {output} -p %d" % THREADS

rule sort_and_index:
  input:
    rules.bowtie2_alignment.output[0],
  output:
    bam="output/{run}/{reference}/sorted.bam",
    bai="output/{run}/{reference}/sorted.bam.bai"
  shell:
    """
      samtools sort -@ %d {input} > {output.bam}
      samtools index {output.bam}
    """ % THREADS

rule qualimap:
  input:
    rules.sort_and_index.output.bam
  output:
    "output/{run}/{reference}/qualimapReport.html"
  params:
    "output/{run}/{reference}"
  shell:
    "qualimap bamqc -bam {input} -outdir {params}"

rule lofreq:
  input:
    rules.reference_index.output[0],
    bam=rules.sort_and_index.output.bam,
    reference=rules.bowtie2_index.input[0]
  output:
    "output/{run}/{reference}/variants.vcf"
  shell:
    "lofreq call-parallel -f {input.reference} -o {output} --pp-threads %d {input.bam}" % THREADS

rule quasirecomb:
  input:
    rules.sort_and_index.output.bam
  output:
    "output/{run}/{reference}/quasirecomb/quasispecies.fasta"
  params:
    basedir="{run}/{reference}/quasirecomb"
  shell:
    """
      quasirecomb -conservative -o {params.basedir} -i {input}
    """ 
