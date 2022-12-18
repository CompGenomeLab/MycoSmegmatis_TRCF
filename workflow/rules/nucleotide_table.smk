project_ = config['project']

rule fastq2fa:
    input: rules.cutadapt_se.output.fastq
    output: "results/{sample}/{sample}_cut.fa"
    log: "logs/rule/analysis/{sample}/{sample}_fastq2fa.log",
    benchmark: "logs/rule/analysis/{sample}/{sample}_fastq2fa.benchmark.txt",
    resources:
        memory="4GB",
        cpu=1
    shell:
        """
        (echo "`date -R`: Converting fastq to fasta..." &&
        python3 workflow/scripts/fastq2fa.py \
        -i {input} \
        > {output} &&
        echo "`date -R`: Success! Fastq to fasta conversion" ||
        {{ echo "`date -R`: Fastq to fasta conversion cannot be done..."; rm {output}; exit 1; }}  ) \
        > {log} 2>&1
        """

rule fa2filteredFa:
    input: 
        fasta=rules.fastq2fa.output,
        length_dist=rules.fastq2length.output.nozero
    output: "results/{sample}/{sample}_cutN.fa"
    log: "logs/rule/analysis/{sample}/{sample}fa2filteredFa.log",
    benchmark: "logs/rule/analysis/{sample}/{sample}fa2filteredFa.benchmark.txt",
    resources:
        memory="4GB",
        cpu=1
    shell:
        """
        (echo "`date -R`: Filtering fasta by sequence length..." &&
        length=$(sort -k2,2n {input.length_dist} | tail -1 | cut -f 1) &&
        python3 workflow/scripts/fa2lengthFilteredFa.py \
        -i {input.fasta} \
        -l $length \
        > {output} &&
        echo "`date -R`: Success! Filtering fasta by sequence length" ||
        {{ echo "`date -R`: Filtering fasta by sequence length cannot be done..."; rm {output}; exit 1; }}  ) \
        > {log} 2>&1
        """

rule nucleotide_table:
    input:
        fasta=rules.fa2filteredFa.output,
        length_dist=rules.fastq2length.output.nozero
    output:
        nuc=temp("results/{sample}/{sample}_nuc.txt"),
        dinuc=temp("results/{sample}/{sample}_dinuc.txt"),
    log:
        "logs/rule/analysis/{sample}/{sample}_nucleotide.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_nucleotide.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        length=$(sort -k2,2n {input.length_dist} | tail -1 | cut -f 1) &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input.fasta} \
        -k 2 \
        -n 13 \
        -l $length \
        -o {output.dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Dinucleotide abundace table cannot be calculated..."; rm {output.dinuc}; exit 1; }}  ) \
        > {log} 2>&1

        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        length=$(sort -k2,2n {input.length_dist} | tail -1 | cut -f 1) &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input.fasta} \
        -k 1 \
        -n 13 \
        -l $length \
        -o {output.nuc} &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Nucleotide abundace table cannot be calculated..."; rm {output.nuc}; exit 1; }}  ) \
        > {log} 2>&1
        """


rule nucleotide_table_forEachLength:
    input:
        fasta=rules.fastq2fa.output,
    output:
        nuc="results/{sample}/lengthSeparated/{sample}_{readLength}.txt",
        dinuc="results/{sample}/lengthSeparated/{sample}_{readLength}_dinuc.txt",
        organized_nuc="results/{sample}/lengthSeparated/{sample}_{readLength}_organized.txt",
        organized_dinuc="results/{sample}/lengthSeparated/{sample}_{readLength}_dinuc_organized.txt",
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample])
    log:
        "logs/rule/analysis/{sample}/{sample}_nucleotide_{readLength}.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_nucleotide_{readLength}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input.fasta} \
        -k 2 \
        -l {wildcards.readLength} \
        -o {output.dinuc} && \
        Rscript workflow/scripts/nucleotideTableChangeFormat.r \
        -i {output.dinuc} \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output.organized_dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Dinucleotide abundace table cannot be calculated..."; rm {output.dinuc}; exit 1; }}  ) \
        > {log} 2>&1

        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input.fasta} \
        -k 1 \
        -l {wildcards.readLength} \
        -o {output.nuc} && \
        Rscript workflow/scripts/nucleotideTableChangeFormat.r \
        -i {output.nuc} \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output.organized_nuc} &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Nucleotide abundace table cannot be calculated..."; rm {output.nuc}; exit 1; }}  ) \
        > {log} 2>&1
        """

rule merge_nucleotide_table_forEachLength:
    input:
        lambda w: input4mergeNucleotideContent(config),
    output:
        out=f"results/{project_}/nucleotide_content.tsv",
    log: f"logs/rule/analysis/{project_}/nucleotide_content.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge nucleotide content files..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """



rule plot_nuc:
    input:
        nuc=rules.nucleotide_table.output.nuc,
        dinuc=rules.nucleotide_table.output.dinuc,
    output:
        nuc=report("results/{sample}/{sample}_nuc.pdf", 
                    category="Nucleotide Content"),
        dinuc=report("results/{sample}/{sample}_dinuc.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.sample, config["meta"][w.sample]["product"]),
        name=lambda w: returnItself(config["meta"][w.sample]["title"]) + "",
    log:
        "logs/rule/figs/{sample}/{sample}_plot_nuc.log",
    benchmark:
        "logs/rule/figs/{sample}/{sample}_plot_nuc.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.nuc} \
        -k 1 \
        -s {params.name} \
        -o {output.nuc} \
        -l {log}

        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.dinuc} \
        -k 2 \
        -s {params.name} \
        -f {params.motif} \
        -o {output.dinuc} \
        -l {log}
        """



rule plot_nuc_sim:
    input:
        dinuc="results/{samples}/simulation/{samples}_dinucleotideTable.txt",
        nuc="results/{samples}/simulation/{samples}_nucleotideTable.txt",
    output:
        dinuc=report("results/{samples}/simulation/{samples}_dinucleotideTable.pdf", 
                    category="Nucleotide Content"),
        nuc=report("results/{samples}/simulation/{samples}_nucleotideTable.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.samples, config["meta"][w.samples]["product"]),
        name=lambda w: returnItself(config["meta"][w.samples]["title"]) + "_simulation",
    log:
        "logs/rule/figs/{samples}/simulation/{samples}_plot_nuc.log",
    benchmark:
        "logs/rule/figs/{samples}/simulation/{samples}_plot_nuc.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.dinuc} \
        -k 2 \
        -s {params.name} \
        -f {params.motif} \
        -o {output.dinuc} \
        -l {log}

        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.nuc} \
        -k 1 \
        -s {params.name} \
        -o {output.nuc} \
        -l {log}
        """