#!/usr/bin/py
import pandas as pd
import os

# trim -> align -> removeDup -> index -> shift reads -> peak call -> bigWig -> filter peaks -> motif analysis

configfile: "config.json"
localrules: all, mkdir

df = pd.read_csv(config["meta_file"], sep='\t', header=0, index_col=0)
sample_ids = list(df.index)
df.index = sample_ids

def get_pair_gz(sample_id):
    dir = config["raw_fastq_gz_dir"]
    return tuple(os.path.join(dir, df.loc[str(sample_id), x]) for x in ('ForwardFastqGZ', 'ReverseFastqGZ'))

def get_forward_primer(sample_id):
    return df.loc[sample_id]["Adapter_1"]

def get_reverse_primer(sample_id):
    return df.loc[sample_id]["Adapter_2"]

rule all:
    input:expand("{dir}/{sample.id}.cpm.norm.bw", dir=config["dir_names"]["bigwigs"],sample_id=sample_ids)
    run:
        for sample in sample_ids:
            print("Wrapping up pipeline")
                   
rule mkdir:
    output: touch(config["file_names"]["mkdir_done"])
    params: dirs = list(config["dir_names"].values())
    shell: "mkdir -p {params.dirs}"

rule trim:
    input: 
        rules.mkdir.output,
        all_read1 = lambda wildcards: get_pair_gz(wildcards.sample_id)[0],
        all_read2 = lambda wildcards: get_pair_gz(wildcards.sample_id)[1]
    output: 
        trimmed_read1 = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R1.fastq.gz",
        trimmed_read2 = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R2.fastq.gz",
        trimmed_stats = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.stats"

    version: config["tool_version"]["cutadapt"]
    params:
        adapter1=lambda wildcards: get_forward_primer(wildcards.sample_id),
        adapter2=lambda wildcards: get_reverse_primer(wildcards.sample_id)
    shell: "cutadapt -m 15 -a {params.adapter1} -A {params.adapter2} -n 2 -o {output.trimmed_read1} -p {output.trimmed_read2} {input.all_read1} {input.all_read2} >{output.trimmed_stats}"

rule map:
    input:
        p1 = rules.trim.output.trimmed_read1,
        p2 = rules.trim.output.trimmed_read2
    output:
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam",
        stats = config["dir_names"]["mapped_dir"] + "/{sample_id}.stats",
    version: config["tool_version"]["bwa"]
    params:
        threads = config["params"]["bowtie2"]["threads"],
        map_all = config["params"]["bowtie2"]["all"],
        reference = config["params"]["bowtie2"]["bwa_reference"]
    shell:
        """
        #bwa mem -t {params.threads} {params.reference} {input.p1} {input.p2} | samtools view -F 4 -Sb |  samtools sort -T {output.mapped_bam_file}.align -o {output.mapped_bam_file}
        bowtie2 --threads {params.threads} -1 {input.p1} -2 {input.p2} -x {params.reference} 2> {output.stats} | samtools sort -T ${sample.id}.tmp -@ 6 -O bam -o {output.mapped_bam_file}
        """

rule remove_dup:
    input:
        sorted_bam = rules.map.output.mapped_bam_file,
    output:
        marked_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.marked.bam",
        dedup_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.dedup.bam",
        dedup_stats = config["dir_names"]["mapped_dir"] + "/{sample_id}.dedup.metrics.txt"
    params:
        picard_path = config["picard_path"]
    shell:
        """
        java -jar {params.picard_path} MarkDuplicates INPUT={input.sorted_bam} OUTPUT={output.marked_bam_file} METRICS_FILE={output.dedup_stats};
        samtools view -bh -F 1024 {output.marked_bam_file} > {output.dedup_bam_file}
        """

rule index_bam:
    input:
        dedup_bam_file = rules.remove_dup.output.dedup_bam_file
    output:
        dedup_bam_index = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam.bai"
    shell:
        """
        samtools index {input.filtered_bam} {output.dedup_bam_index}
        """

rule shift_bam:
    input:
        dedup_bam_file = rules.remove_dup.output.dedup_bam_file,
        dedup_bam_index = rules.index_bam.output.dedup_bam_index
    output:
        shifted_bam = config["dir_names"]["mapped_dir"] + "/{sample.id}.shifted.bam"
    params:
        cores = config["params"]["deeptools"]["cores"]
    shell:
	"""
        alignmentSieve --numberOfProcessors {params.cores} -b {input.dedup_bam} -o {output.shifted_bam} --ATACshift  
	"""

rule sort_shifted_bam:
    input:
        shifted_bam = rules.shift_bam.output.shifted_bam
    output:
        sorted_shifted_bam = config["dir_names"]["mapped_dir"] + "/{sample.id}.shifted.sorted.bam"
        sorted_shifted_index = config["dir_names"]["mapped_dir"] + "/{sample.id}.shifted.sorted.bam.bai"
    params:
        cores = config["params"]["deeptools"]["cores"]
    shell:
        """
        samtools sort -@ {params.cores} -T {sample.id}.tmp {input.shifted_bam} > {output.sorted_shifted_bam};
        samtools index {output.sorted_shifted_bam} {output_shifted_index}
        """

rule call_peaks:
    input:
        sorted_shifted_bam = rules.sort_shifted_bam.output.sorted_shifted_bam
    output:
        narrowPeak = config["dir_names"]["peak_calling"] + "/{sample.id}.narrowPeak",
        broadPeak = config["dir_names"]["peak_calling"] + "/{sample.id}.broadPeak"
    params:
        genome = config["params"]["macs2"]["genome"],
        output_directory = config["dire_names"]["peak_calling"]
    shell:
        """
        macs2 callpeak -t {input.sorted_shifted_bam} -g {params.genome} -f BAMPE -n {sample.id} --outdir {params.output_directory} -q 0.01 -B --SPMR --keep-dup all;
        macs2 callpeak -t {input.sorted_shifted_bam} -g {params.genome} -f BAMPE -n {sample.id} --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all
        """

rule normalized_bigwig:
    input:
        sorted_shifted_bam = rules.sort_shifted_bam.output.sorted_shifted_bam
    output:
        normalized_bigwig = config["dir_names"]["bigwigs"] + "/{sample.id}.cpm.norm.bw"
    params:
        cores = config["params"]["deeptools"]["cores"],
        size = config["params"]["deeptools"]["size"]
    shell:
        """
        bamCoverage --bam {input.sorted_shifted_bam} -o {params.normalized_bigwigs} \
            --binSize 10 \
            --normalizeUsing CPM \
            --effectiveGenomeSize {params.size} \
            --numberOfProcessors {params.cores}
        """


