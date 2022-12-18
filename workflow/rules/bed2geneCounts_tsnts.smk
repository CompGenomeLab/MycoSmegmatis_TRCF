import re
build_ = config["build"]
project_ = config["project"]

rule bed2geneCounts_tsnts:
    input:
        reads=rules.fixTo12.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output:
        TS=temp("results/{sample}/{sample}_TS.tsv"),
        NTS=temp("results/{sample}/{sample}_NTS.tsv"),
        mergedTSNTS="results/{sample}/{sample}_TSNTS.tsv",
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2geneTSNTScounts.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bed2geneTSNTScounts.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample])
    conda:
        "../envs/bedtools.yaml"
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Counting TS reads {output.TS}..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -S \
        > {output.TS} &&
        echo "`date -R`: Success! {output.TS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Counting NTS reads {output.NTS}..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -s \
        > {output.NTS} &&
        echo "`date -R`: Success! {output.NTS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Merging TS NTS files..." &&
        paste \
        {output.TS} \
        <(cut -f 7 {output.NTS}) \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output.mergedTSNTS} &&
        echo "`date -R`: Success! {output.mergedTSNTS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        """
        )

use rule bed2geneCounts_tsnts as bed2geneCounts_tsnts_random with:
    input:
        reads=rules.randombed.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output:
        TS=temp("results/{sample}/random/{sample}_TS.tsv"),
        NTS=temp("results/{sample}/random/{sample}_NTS.tsv"),
        mergedTSNTS="results/{sample}/random/{sample}_TSNTS.tsv",
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2geneTSNTScounts_random.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bed2geneTSNTScounts_random.benchmark.txt",
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample])


rule bed2geneCounts_tsnts_mappable:
    input:
        reads=rules.mappableReads.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output:
        TS=temp(f"results/{project_}/mappable_TS.tsv"),
        NTS=temp(f"results/{project_}/mappable_NTS.tsv"),
        mergedTSNTS=f"results/{project_}/mappable_TSNTS.tsv",
    log:
        f"logs/rule/analysis/{project_}/mappable_bed2geneTSNTScounts_random.log",
    benchmark:
        f"logs/rule/analysis/{project_}/mappable_bed2geneTSNTScounts_random.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    conda:
        "../envs/bedtools.yaml"
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Counting TS reads {output.TS}..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -S \
        > {output.TS} &&
        echo "`date -R`: Success! {output.TS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Counting NTS reads {output.NTS}..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -s \
        > {output.NTS} &&
        echo "`date -R`: Success! {output.NTS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Merging TS NTS files..." &&
        paste \
        {output.TS} \
        <(cut -f 7 {output.NTS}) \
        > {output.mergedTSNTS} &&
        echo "`date -R`: Success! {output.mergedTSNTS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        """
        )


rule mergeTSNTScounts:
    input:
        lambda w: input4mergeTSNTScounts(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"]),
    output:
        out=f"results/{project_}/readCountsTSNTS.tsv",
    log: f"logs/rule/analysis/{project_}/readCountsTSNTS.log",
    benchmark: f"logs/rule/analysis/{project_}/readCountsTSNTS.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge TSNTS counts from true sample..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """

rule mergeTSNTScounts_random:
    input:
        lambda w: input4mergeTSNTScounts_random(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"]),
    output:
        out=f"results/{project_}/random/readCountsTSNTS.tsv",
    log: f"logs/rule/analysis/{project_}/random/readCountsTSNTS.log",
    benchmark: f"logs/rule/analysis/{project_}/random/readCountsTSNTS.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge TSNTS counts from true sample..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """

use rule bed2geneCounts_tsnts as bed2geneCounts_tsnts_sim with:
    input:
        reads="results/{sample}/simulation/{sample}.bed",
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output:
        TS=temp("results/{sample}/simulation/{sample}_TS.tsv"),
        NTS=temp("results/{sample}/simulation/{sample}_NTS.tsv"),
        mergedTSNTS="results/{sample}/simulation/{sample}_TSNTS.tsv",
    log:
        "logs/rule/analysis/{sample}/simulation/{sample}_bed2geneCounts_tsnts.log",
    benchmark:
        "logs/rule/analysis/{sample}/simulation/{sample}_bed2geneCounts_tsnts.benchmark.txt",


use rule mergeTSNTScounts as mergeTSNTScounts_sim with:
    input:
        lambda w: input4mergeTSNTScounts_sim(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"]),
    output:
        out="results/simulation/readCountsTSNTS.tsv",
    log:
        "logs/rule/analysis/simulation/readCountsTSNTS.log",
    benchmark:
        "logs/rule/analysis/simulation/readCountsTSNTS.benchmark.txt",