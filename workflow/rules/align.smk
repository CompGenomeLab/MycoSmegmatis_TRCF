build_ = config["build"]
project_ = config["project"]

rule bowtie2_se:
    input:
        sample = rules.cutadapt_se.output.fastq,
        bowtie2=f"resources/ref_genomes/{build_}/Bowtie2/genome.1.bt2",
    output:
        sam=temp("results/{sample}/{sample}.sam"),
        bam="results/{sample}/{sample}.bam",
    params:
        ref_genome=rules.genome_build.params.base,
        extra="--seed 1",
        sample_id="{sample}"
    log:
        "logs/rule/analysis/{sample}/{sample}_bowtie2.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bowtie2.benchmark.txt",
    resources:
        memory="64GB",
        cpu=16
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq file..." &&
        bowtie2 \
        --threads {resources.cpu} \
        {params.extra} \
        -x {params.ref_genome} \
        --rg-id {params.sample_id} \
        --rg SM:{params.sample_id} \
        -U {input.sample} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh {output.sam} | samtools sort -o {output.bam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bam2bed_se:
    input:
        bam= rules.bowtie2_se.output.bam,
        fasta="resources/ref_genomes/" + config["build"] + "/singletons.fa"
    output:
        bed="results/{sample}/{sample}.bed",
        bam="results/{sample}/{sample}_sorted.bam",
        sam="results/{sample}/{sample}_sorted.sam",
        idx="results/{sample}/{sample}_sorted.bam.bai",
    params:
        q_trim=config["samtools_q_trim_se"], 
    log:
        "logs/rule/analysis/{sample}/{sample}_bam2bed.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bam2bed.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools view -h {input.bam} | python3 workflow/scripts/filterSam.py -f {input.fasta} -o {output.sam} && 
        samtools view -S {output.sam} -b | samtools sort > {output.bam} &&
        echo "`date -R`: Success! Bam file is filtered and sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} -b {output.bam} |&
        bedtools bamtobed | sort -k1,1 -k2,2n -k3,3n > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1


        """

rule fixTo12:
    input: rules.bam2bed_se.output.bed,
    output: "results/{sample}/random/{sample}_fix12.bed"
    log:
        "logs/rule/analysis/{sample}/{sample}_fixTo12.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_fixTo12.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:
        """
            (echo "`date -R`: Filtering bed file to clip last nucleotide of 13 mer" && 
            awk '{{if($3 - $2 == 12 || $3 - $2 == 13) print}}' {input} |  awk '{{if ($3 - $2 == 13 && $6 == "+") $3 = int($3) - 1; else if ($3 - $2 == 13 && $6 == "-") $2 = int($2) + 1; print}}' > {output} &&
            echo "`date -R`: Success! We clipped last nuc of 13mer." || 
            {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule randombed:
    input: 
        bed=rules.fixTo12.output,
        singleton="results/" + config["project"] + "/mappableReads.bed",
        # singleton="resources/ref_genomes/" + config["build"] + "/singletons.fa",
    output: "results/{sample}/random/{sample}_randomMappable.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_randombed.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_randombed.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Processing bed file to get random headers..." && 
        python3 workflow/scripts/random_reads.py -i {input.singleton} -b {input.bed} | sort -k2n > {output} &&
        echo "`date -R`: Success! We obtained the random bed file." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
    """


# rule uniqbed:
#     input: 
    
#     output: "results/{sample}/uniq/{sample}.bed",
#     log:
#         "logs/rule/analysis/{sample}/{sample}_uniqbed.log",
#     benchmark:
#         "logs/rule/analysis/{sample}/{sample}_uniqbed.benchmark.txt",
#     resources:
#         memory="16GB",
#         cpu=1
#     conda:
#         "../envs/bedtools.yaml"
#     shell:  
#         """
#         (echo "`date -R`: Processing bed file to get unique locations..." && 
#         awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t.\\t\"$5\"\\t\"$6}}' {input} | uniq > {output} &&
#         echo "`date -R`: Success! We obtained the unique bed file." || 
#         {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
#     """


rule sep_strands:
    input:
        rules.fixTo12.output,
    output:
        plus="results/{sample}/{sample}_plus.bed",
        minus="results/{sample}/{sample}_minus.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_sepStrands.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_sepStrands.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """


rule sep_byLength:
    input:
        "results/{sample}/{sample}_{strand}.bed",
    output:
        "results/{sample}/lengthSeparated/{sample}_{strand}_{readLength}.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_{readLength}_sepByLength_{strand}.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_{readLength}_sepByLength_{strand}.benchmark.txt",
    resources:
        memory="4GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Separating reads by their lengths..." &&
        awk '($3 - $2)=={wildcards.readLength}{{print}}' {input} > {output} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  )  > {log} 2>&1
        """

    
rule length_dist_eachLength:
    input:rules.fixTo12.output,
    output:
        "results/{sample}/{sample}_bedLength.txt",
    log:
        "logs/rule/analysis/{sample}/{sample}_bedLength.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bedLength.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Calculating the read length distribution..." &&
        awk '{{print $3-$2}}' {input} |&
        sort -k1,1n |& 
        uniq -c |& 
        sed 's/\s\s*/ /g' |&
        awk '{{print $2"\\t"$1}}' > {output} &&
        echo "`date -R`: Success! Length distribution is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

# use rule sep_byLength_plus as sep_byLength_minus with:
#     input:
#         rules.sep_strands.output.minus,
#     output:
#         "results/{sample}/{sample}_minus_{readLength}.bed",
#     log:
#         "logs/rule/analysis/{sample}/{sample}_{readLength}_sepByLength_minus.log",
#     benchmark:
#         "logs/rule/analysis/{sample}/{sample}_{readLength}_sepByLength_minus.benchmark.txt",


use rule sep_strands as sep_strands_sim with:
    input: 
        "results/{sample}/simulation/{sample}.bed"
    output:
        plus="results/{sample}/simulation/{sample}_plus.bed",
        minus="results/{sample}/simulation/{sample}_minus.bed",
    log:
        "logs/rule/analysis/{sample}/simulation/{sample}_sepStrands.log",
    benchmark:
        "logs/rule/analysis/{sample}/simulation/{sample}_sepStrands.benchmark.txt",

rule genomecov:
    input:
        bed="results/{sample}/{sample}_{strand}.bed", 
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg=temp("results/{sample}/{sample}_{strand}.bdg"),
        bw="results/{sample}/{sample}_{strand}.bw",
    params:
        read=lambda w, input: mappedReads(input['bed']),
    log:
        "logs/rule/analysis/{sample}/{sample}_{strand}_genomeCov.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_{strand}_genomeCov.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell: 
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Calculating genome coverage of {input.bed}..." &&
        bedtools genomecov \
        -i {input.bed} \
        -g {input.ref_genome_index} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.bdg} \
        && \
        bedGraphToBigWig {output.bdg} {input.ref_genome_index} {output.bw} &&
        echo "`date -R`: Success! Genome coverage for located damaged is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
        )


rule mergeReplicates:
    input: lambda w: input4mergeReplicates(config, w.strain, w.strand),
    output: "results/mergedReplicates/{strain}_{strand}.bed"
    log:
        "logs/rule/analysis/mergeReplicates_{strain}_{strand}.log",
    benchmark:
        "logs/rule/analysis/mergeReplicates_{strain}_{strand}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell: 
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Merging replicates {input}..." &&
        cat {input} | sort -k2,2n > {output} &&
        echo "`date -R`: Success! Genome coverage for located damaged is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
        )

rule bigwigWithMergedBedReplicates:
    input: 
        bed=rules.mergeReplicates.output,
        ref_genome_index=rules.samtools_index.output,
    output: 
        bdg=temp("results/mergedReplicates/{strain}_{strand}.bdg"),
        bw=report("results/mergedReplicates/{strain}_{strand}.bw", category="bigwig")
    params:
        read=lambda w, input: mappedReads(input['bed']),
    log:
        "logs/rule/analysis/bigwigWithMergeBedReplicates_{strain}_{strand}.log",
    benchmark:
        "logs/rule/analysis/bigwigWithMergeBedReplicates_{strain}_{strand}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell: 
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Calculating genome coverage of {input.bed}..." &&
        bedtools genomecov \
        -i {input.bed} \
        -g {input.ref_genome_index} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.bdg} \
        && \
        bedGraphToBigWig {output.bdg} {input.ref_genome_index} {output.bw} &&
        echo "`date -R`: Success! Genome coverage for located damaged is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
        )

use rule genomecov as genomecov_readlength with:
    input:
        bed="results/{sample}/lengthSeparated/{sample}_{strand}_{readLength}.bed", 
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg=temp("results/{sample}/lengthSeparated/{sample}_{strand}_{readLength}.bdg"),
        bw="results/{sample}/lengthSeparated/{sample}_{strand}_{readLength}.bw",
    params:
        read=lambda w, input: mappedReads(input['bed']),
    log:
        "logs/rule/analysis/{sample}/{sample}_{strand}_{readLength}_genomeCov.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_{strand}_{readLength}_genomeCov.benchmark.txt",

rule mappableReads:
    input: input4mappableReads(config['sample'])
    output: f'results/{project_}/mappableReads.bed'
    log:
        "logs/rule/mappableReads.log",
    benchmark:
        "logs/rule/mappableReads.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell: 
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Retrieving the mappable reads from {input}..." &&
        cat {input} | sort -k1,1 -k2,2n -k3,3n -k6 -u > {output} &&
        echo "`date -R`: Mappable reads are generated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
        )

# use rule genomecov_plus as genomecov_minus with:
#     input:
#         bed=rules.sep_strands.output.minus, 
#         ref_genome_index=rules.samtools_index.output,
#     output:
#         bdg=temp("results/{sample}/{sample}_neg.bdg"),
#         bw="results/{sample}/{sample}_neg.bw",
#     log:
#         "logs/rule/analysis/{sample}/{sample}_genomeCovN.log",
#     benchmark:
#         "logs/rule/analysis/{sample}/{sample}_genomeCovN.benchmark.txt",

# use rule genomecov_plus as genomecov_plus_sim with:
#     input:
#         bed=rules.sep_strands.output.plus, 
#         ref_genome_index=rules.samtools_index.output,
#     output:
#         bdg=temp("results/{sample}/simulation/{sample}_pos.bdg"),
#         bw=report("results/{sample}/simulation/{sample}_pos.bw"),
#     params:
#         read=lambda w, input: mappedReads(input['bed']),
#     log:
#         "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovP.log",
#     benchmark:
#         "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovP.benchmark.txt",


# use rule genomecov_minus as genomecov_minus_sim with:
#     input:
#         bed=rules.sep_strands.output.minus, 
#         ref_genome_index=rules.samtools_index.output,
#     output:
#         bdg=temp("results/{sample}/simulation/{sample}_neg.bdg"),
#         bw=report("results/{sample}/simulation/{sample}_neg.bw"),
#     params:
#         read=lambda w, input: mappedReads(input['bed']),
#     log:
#         "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovN.log",
#     benchmark:
#         "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovN.benchmark.txt"

