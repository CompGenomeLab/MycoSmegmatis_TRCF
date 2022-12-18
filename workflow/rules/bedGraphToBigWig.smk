
rule bedGraphToBigWig:
    input:
        bdg="results/{sample}/{sample}_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome.fa.fai",
    output:
        report("results/{sample}/{sample}_{strand}.bw", 
                category="BigWig"),
    log:
        "logs/rule/analysis/{sample}/{sample}_bedGraphToBigWig_{strand}.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bedGraphToBigWig_{strand}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

use rule bedGraphToBigWig as bedGraphToBigWig_readLength with:
    input:
        bdg="results/{sample}/{sample}_{strand}_{readLength}.bdg",
        index="resources/ref_genomes/{build}/genome.fa.fai",
    output: "results/{sample}/{sample}_{strand}_{readLength}.bw", 
        # report("results/{sample}/{sample}_{strand}_{readLength}.bw", 
                # category="BigWig"),
    log:
        "logs/rule/analysis/{sample}/{sample}_bedGraphToBigWig_{strand}_{readLength}.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bedGraphToBigWig_{strand}_{readLength}.benchmark.txt",
  

rule bedGraphToBigWig_sim:
    input:
        bdg="results/{samples}/simulation/{samples}_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome.fa.fai",
    output:
        report("results/{samples}/simulation/{samples}_{strand}.bw", 
                category="BigWig"),
    log:
        "logs/rule/analysis/{samples}/simulation/{samples}_bedGraphToBigWig_{strand}.log",
    benchmark:
        "logs/rule/analysis/{samples}/simulation/{samples}_bedGraphToBigWig_{strand}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """