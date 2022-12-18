rule prepareSingletons:
    input:
        "resources/ref_genomes/{build}/genome.fa",
    output:
        theoretical="resources/ref_genomes/{build}/theoreticalReads.fa",
        singletons="resources/ref_genomes/{build}/singletons.fa",
        bedsingletons=report("resources/ref_genomes/{build}/singletons.bed", category="genome"),
    log:
        "logs/rule/analysis/{build}/log/prepareSingletons.log",
    benchmark:
        "logs/rule/analysis/{build}/log/prepareSingletons.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    params:
        kmerlist=config['readLength']
    shell:  
        """
        (echo "`date -R`: Preparing theoretical reads..." &&
        python3 workflow/scripts/fa2theoreticalReads.py -i {input} -k {params.kmerlist} > {output.theoretical} &&
        echo "`date -R`: Success! Theoretical reads file is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Preparing singletons (reads mapped onto genome only for once)..." &&
        python3 workflow/scripts/fa2singletons.py -i {output.theoretical} -b {output.bedsingletons} > {output.singletons} &&
        echo "`date -R`: Success! Singletons file is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """