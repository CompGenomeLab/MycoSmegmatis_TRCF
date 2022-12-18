rule check_presence:
    output: 
        gtf = report(f"resources/ref_genomes/{build}/genome.gtf", category="genome")