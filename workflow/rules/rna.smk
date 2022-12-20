build_ = config["build"]
scratch_dir = config['scratch_dir']

rule get_fastq_pe_gz:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "results/rna/pe/raw/{accession}_1.fastq.gz",
        "results/rna/pe/raw/{accession}_2.fastq.gz",
    log:
        "logs/pe/{accession}.gz.log"
    params:
        extra="--skip-technical --temp " + scratch_dir
    threads: 6  # defaults to 6
    resources:
        memory="64GB",
        cpu=6
    wrapper:
        "v1.19.1/bio/sra-tools/fasterq-dump"

rule trimmomatic_pe_rna:
    input:
        r1="results/rna/pe/raw/{rna_sample}_1.fastq.gz",
        r2="results/rna/pe/raw/{rna_sample}_2.fastq.gz"
    output:
        r1="results/rna/pe/trimmed/{rna_sample}_1.fastq.gz",
        r2="results/rna/pe/trimmed/{rna_sample}_2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/{rna_sample}_1.unpaired.fastq.gz",
        r2_unpaired="trimmed/{rna_sample}_2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{rna_sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        16
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        memory="1GB",
        cpu=16
    wrapper:
        "v1.19.1/bio/trimmomatic/pe"

rule bowtie2_rna:
    input:
        sample=[rules.trimmomatic_pe_rna.output.r1, rules.trimmomatic_pe_rna.output.r2],
        idx=multiext(
            f"resources/ref_genomes/{build_}/Bowtie2/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "results/rna/pe/mapped/{rna_sample}.bam",
    log:
        "logs/bowtie2/rna/{rna_sample}.log",
    params:
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    resources:
        memory="16GB",
        cpu=8
    wrapper:
        "v1.19.1/bio/bowtie2/align"

rule bam2bed_pe_rna:
    input: rules.bowtie2_rna.output,
    output: "results/rna/pe/mapped/{rna_sample}.bed"
    log:
        "logs/bam2bed/{rna_sample}.log",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (echo "`date -R`: Processing bam file..." && 
        bedtools bamtobed -bedpe -i {input} | grep -v "^\." | awk '{{ print $1"\\t"$2"\\t"$6"\\t"$7"\\t\\.\\t+" }}' | sort -k1,1 -k2,2n -k3,3n > {output} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bed2geneCounts_rna:
    input:
        reads=rules.bam2bed_pe_rna.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output: "results/rna/pe/mapped/{rna_sample}.tsv",
    params:
        fields = lambda w: getSampleFields(config["rna_meta"][w.rna_sample])
    log: "logs/bed2geneCounts_rna/{rna_sample}.log",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Counting reads..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -F 0.5 \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output} &&
        echo "`date -R`: Success! {output} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        """
        )

rule mergeTSNTScounts_rna:
    input:
        lambda w: rna_input4merge(config["rna_sample"]),
    output:
        out=f"results/rna/readCounts.tsv",
    log: f"logs/rna/readCounts.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge counts from samples..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """
