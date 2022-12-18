rule bgzip_variation:
    input:
        "resources/ref_genomes/{build}/variation.vcf.gz",
    output:
        fin="resources/ref_genomes/{build}/variation_bgzip.vcf.gz",
        decomp=temp("resources/ref_genomes/{build}/variation.vcf"),
    log:
        "logs/rule/analysis/{build}/log/bgzip_variation.log",
    benchmark:
        "logs/rule/analysis/{build}/log/bgzip_variation.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bgzip_variation.yaml"
    shell:  
        """
        (echo "`date -R`: Decompressing vcf file..." &&
        gunzip -c {input} > {output.decomp} &&
        echo "`date -R`: Success! vcf file is decompressed." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Bgzip vcf file..." &&
        bgzip -c {output.decomp} > {output.fin} &&
        echo "`date -R`: Success! vcf file is compressed." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """