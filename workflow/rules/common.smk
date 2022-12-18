#!/bin/env python

import warnings
import os
import subprocess
import re
# configfile: "config/config_XR_initial.yaml" 

################### Helper Functions ###########################################

def getSRR(sample, srrList, sampleList):

    try:
        idx = sampleList.index(sample)
    except:
       raise(ValueError("Designated wildcard cannot be found in sample list."))
        
    return srrList[idx]

def isSingle(sample, sampleList, srrEnabled, srrList, sample_dir):

    if srrEnabled:

        mySRR = getSRR(sample, srrList, sampleList)

        if mySRR == "NA":

            single = f"{sample_dir}{sample}.fastq.gz"
            pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
            paired1 = f"{sample_dir}{sample}_1.fastq.gz"
            
            if os.path.isfile(pairedR1) or os.path.isfile(paired1):
                return False
            elif os.path.isfile(single):
                return True
            else:
                raise(ValueError(f"{paired1}, {pairedR1}, or {single} not found..."))

        if ":" in mySRR:
            mySRR = mySRR.split(":")[0]

        shellCommand = f'fastq-dump -X 1 -Z --split-spot {mySRR} | wc -l'
        #print(shellCommand)
        p=subprocess.getoutput(shellCommand)
        #print(p)
        lineNum = int(p.split("\n")[2])
        #print(lineNum)

        if lineNum == 4:
            return True
        else:
            return False

    else:

        single = f"{sample_dir}{sample}.fastq.gz"
        pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
        paired1 = f"{sample_dir}{sample}_1.fastq.gz"
        
        if os.path.isfile(pairedR1) or os.path.isfile(paired1):
            return False
        elif os.path.isfile(single):
            return True
        else:
            raise(ValueError(f"{paired1}, {pairedR1}, or {single} not found..."))

def getPaired(sample, sampleList, read, sample_dir):

    pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
    paired1 = f"{sample_dir}{sample}_1.fastq.gz"
    
    if os.path.isfile(pairedR1) and read == "forward":
        return f"{sample_dir}{sample}_R1.fastq.gz"

    elif os.path.isfile(pairedR1) and read == "reverse":
        return f"{sample_dir}{sample}_R2.fastq.gz"

    elif os.path.isfile(paired1) and read == "forward":
        return f"{sample_dir}{sample}_1.fastq.gz"

    elif os.path.isfile(paired1) and read == "reverse":
        return f"{sample_dir}{sample}_2.fastq.gz"

def input4filter(wildcards, sampleList, srrEnabled, srrList):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList, "resources/samples/"):
        return "results/{sample}/{sample}_se.bed"
    else:    
        return "results/{sample}/{sample}_pe.bed"

def input4filter_sim(wildcards, sampleList, srrEnabled, srrList):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList, "resources/samples/"):
        return "results/simulation/{sample}/{sample}.bed"
    else:    
        return "results/simulation/{sample}/{sample}.bed"

def input4fasta(wildcards, sampleList, srrEnabled, srrList):

    if isSingle(wildcards.samples, "input", sampleList, srrEnabled, srrList, "resources/input/"):
        return "results/input/{sample}/{sample}_se.bed"
    else:    
        return "results/input/{sample}/{sample}_pe.bed"

def input4PCA(sampleList, srrEnabled, srrList, build, duplicate):

    inputList = []
    for sample in sampleList:

        if isSingle(sample, sampleList, srrEnabled, srrList, "resources/samples/"):
            inputList.append(f"results/{sample}/{sample}_se_sortedbyCoordinates.bam")
        else:    
            inputList.append(f"results/{sample}/{sample}_pe_sortedbyCoordinates.bam")

    return inputList

def input4computeMatrix():
    inputList = []
    for strand in ["plus", "minus"]:
        inputList.append(f"results/{{samples}}/{{samples}}_{{build}}_{{duplicate}}_sorted_located_{strand}.bw")        
    return inputList

def input4computeMatrix_sim():
    inputList = []
    for strand in ["plus", "minus"]:
        inputList.append(f"results/{{samples}}/simulation/{{samples}}_{{build}}_{{duplicate}}_sorted_located_{strand}.bw")        
    return inputList

def input4mergeTSNTScounts(sampleList, srrEnabled, srrList):
    inputList = []
    for sample in sampleList:
        inputList.append(f"results/{sample}/{sample}_TSNTS.tsv")
    return inputList

def input4mergeNucleotideContent(config):
    inputList = []
    for sample in config['sample']:
        for readLength in config['readLengthForNucleotide']:
            inputList.append(f"results/{sample}/lengthSeparated/{sample}_{readLength}_organized.txt")
    return inputList

def input4mergeLength(config):
    inputList = []
    for sample in config['sample']:
        for readLength in config['readLengthForNucleotide']:
            inputList.append(f"results/{sample}/{sample}_length.txt")
    return inputList

def input4mergeReplicates(config, strain, strand):
    d = {}    
    for sample in config['sample']:
        strain_ = config['meta'][sample]['strain']
        if strain_ in d.keys():
            d[strain_].append(f"results/{sample}/{sample}_{strand}.bed")
        else:
            d[strain_] = [f"results/{sample}/{sample}_{strand}.bed"]
    return d[strain]

def input4mergeTSNTScounts_sim(sampleList, srrEnabled, srrList):
    inputList = []
    for sample in sampleList:
        inputList.append(f"results/{sample}/simulation/{sample}_TSNTS.tsv")
    return inputList

def input4mergeTSNTScounts_random(sampleList, srrEnabled, srrList):
    inputList = []
    for sample in sampleList:
        inputList.append(f"results/{sample}/random/{sample}_TSNTS.tsv")
    return inputList

def input4mappableReads(sampleList):
    inputList = []
    for sample in sampleList:
        inputList.append(f"results/{sample}/{sample}.bed")
    return inputList

def input4nucTable():
    return "results/{sample}/{sample}_lengthMode.fa"


def getMotif(sample, product):

    if product.lower() in ["oxaliplatin", "cisplatin", "bpdedg"]: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif product.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'.{4}(c|t|C|T){2}.{4}'"

    elif product.lower() == "na":
        return "'.{10}'"

def getDinuc(sample, product):

    if product.lower() in ["oxaliplatin", "cisplatin"]: 
        return "'GG'"
    
    elif product.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'CC','CT','TC','TT'"

def returnItself(x):
    return x

def getSampleFields(sample):
    l = []
    for key in sorted(sample.keys()):
        l.append(str(sample[key]))
    return '\t'.join(l)

def getInput(sample, inputExist, inputList, inputIdx, sampleList, build):

    if inputExist:
        inpDict={}
        for inp_idx in range(len(inputIdx)):
            idx_split = inputIdx[inp_idx].strip().split(",")
            indexList=[]
            for sample_idx in idx_split:
                sample_idx = sample_idx.strip() 
                if "-" in sample_idx:
                    for range_idx in range(int(sample_idx.split("-")[0]), int(sample_idx.split("-")[1])+1):
                        indexList.append(int(range_idx)) 
                else:
                    indexList.append(int(sample_idx)) 
                
            for sample_idx in indexList:
                if inputList[inp_idx] not in inpDict:
                    inpDict[inputList[inp_idx]] = [sampleList[sample_idx]]
                else:
                    inpDict[inputList[inp_idx]].append(sampleList[sample_idx])
        for k,v in inpDict.items():
        
            if sample in v:
            
                return f"results/input/{k}/{k}.fasta"
                
    else:
        return f"resources/ref_genomes/{build}/genome.ron" 

def lineNum(file):
    
    linenum = 0
    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

    warnMessage = (f"\n{file} file is either empty or does not exists!\n" + 
        "It is expected if this is a dry-run. The file will be produced " + 
        "after the execution.")

    if linenum == 0:
        warnings.warn(warnMessage)

    return linenum

def mappedReads(*files):

    lineNumber = 0
    for file in files:
        lineNumber += lineNum(str(file))

    return lineNumber

def allInput(config):

    build = config["build"]
    sampleList = config["sample"]
    project = config['project']
    strands = config['strand']
    readLengths = config['readLength']
    readLengthForNucleotide = config['readLengthForNucleotide']


    strains = set()
    for sample in sampleList:
        strains.add(config['meta'][sample]['strain'])

    inputList = []

    # gtf = report(f"resources/ref_genomes/{build}/genome.gtf", category="genome")

    # inputList.append(f"resources/ref_genomes/{build}/operons.bed") 
    inputList.append(f"resources/ref_genomes/{build}/singletons.bed") 
    inputList.append(f"resources/ref_genomes/{build}/genome.gtf") 
    inputList.append(f"results/{project}/mappableReads.bed") 

    duplications = ['nodedup', 'dedup']

    for sample in sampleList:
        sample_dir = f"results/{sample}/"

        
        inputList.append(f"{sample_dir}{sample}_removedDup.fastq")
        inputList.append(f"{sample_dir}{sample}_cut.fastq")
        inputList.append(f"{sample_dir}{sample}_length_distribution.pdf")
        inputList.append(f"{sample_dir}{sample}_dinuc.pdf")
        inputList.append(f"{sample_dir}{sample}_sorted.bam")
        inputList.append(f"{sample_dir}{sample}_TSNTS.tsv")
        inputList.append(f"{sample_dir}{sample}_report.txt")
        # inputList.append(f"{sample_dir}uniq/{sample}.bed")
        inputList.append(f"{sample_dir}{sample}_bedLength.txt")
        inputList.append(f"{sample_dir}random/{sample}_random.bed")
        inputList.append(f"results/{project}/readCountsTSNTS.tsv")
        inputList.append(f"results/{project}/random/readCountsTSNTS.tsv")
        inputList.append(f"results/{project}/mappable_TSNTS.tsv")
        inputList.append(f"results/{project}/nucleotide_content.tsv")
        inputList.append(f"results/{project}/length.tsv")
        inputList.append(f"{sample_dir}{sample}_TSNTS.tsv")
        for strand in strands:
            inputList.append(f"{sample_dir}{sample}_{strand}.bw")
            for strain in list(strains):
                inputList.append(f"results/mergedReplicates/{strain}_{strand}.bw")
            
            # for readLength in readLengths:
            #     inputList.append(f"{sample_dir}lengthSeparated/{sample}_{strand}_{readLength}.bw")
        # for readLength in readLengthForNucleotide:
        #     inputList.append(f"{sample_dir}lengthSeparated/{sample}_{readLength}_organized.txt")
        #     inputList.append(f"{sample_dir}lengthSeparated/{sample}_{readLength}_dinuc_organized.txt")

        # for duplication in duplications:
            # inputList.append(f"{sample_dir}{sample}_{duplication}_sorted_nucleotideTable.pdf")
            # inputList.append(f"{sample_dir}{sample}_{duplication}_sorted_dinucleotideTable.pdf")
            # inputList.append(f"{sample_dir}simulation/{sample}_{duplication}_nucleotideTable.pdf")
            # inputList.append(f"{sample_dir}simulation/{sample}_{duplication}_dinucleotideTable.pdf")
            # inputList.append(f"{sample_dir}{sample}_{duplication}_length_distribution.pdf")
            # inputList.append(f"{sample_dir}{sample}_{duplication}_fslp.bw")
            # inputList.append(f"{sample_dir}{sample}_{duplication}_fslm.bw")
            # inputList.append(f"{sample_dir}simulation/{sample}_{duplication}_fslp.bw")
            # inputList.append(f"{sample_dir}simulation/{sample}_{duplication}_fslm.bw")
            # inputList.append(f"results/readCountsTSNTS_{duplication}.tsv")
            # inputList.append(f"results/plots/{sample}_{duplication}_TSSprofile.png")
            # inputList.append(f"results/plots/plot_heatmap/{sample}_{duplication}_heatmap.png")
            
            # inputList.append(f"results/plots/simulation/{sample}_{duplication}_TSSprofile.png")

            
    #         for strand in ['plus', 'minus']:
    #             inputList.append(f"{sample_dir}{sample}_{duplication}_sorted_{strand}.bed") 
    #             inputList.append(f"{sample_dir}{sample}_{duplication}_sorted_located_{strand}.bed")
    #             inputList.append(f"{sample_dir}{sample}_{duplication}_sorted_located_{strand}_matrix.gz")

            # inputList.append(f"results/scatterplot_PearsonCorr_bigwigScores_{duplication}.pdf")
            # inputList.append(f"results/PearsonCorr_bigwigScores_{duplication}.tab")
            # inputList.append(f"results/heatmap_SpearmanCorr_readCounts_{duplication}.pdf")
            # inputList.append(f"results/SpearmanCorr_readCounts_{duplication}.tab")
            # inputList.append(f"results/PCA_readCounts_{duplication}.pdf")
            # inputList.append(f"results/readCountsTSNTS_{duplication}.tsv")
            # inputList.append(f"results/simulation/readCountsTSNTS_{duplication}.tsv")

    # for i in inputList:
    #     print(i)    

    return inputList


def chromosome():
    keyword = config['build']
    initialTwo = keyword[0:2].upper()
    chromosome = initialTwo + '_' + keyword[2:] + '.1'
    return chromosome




wildcard_constraints:
    duplicate='dedup|nodedup'

# rule all:
#     input:
#         lambda w: allInput(config["build"], config["sample"], 
#             config["srr"]["enabled"], config["srr"]["codes"]),
 
include: "fastqc.smk"

include: "genome_build.smk"
    
include: "prepareSingletons.smk"

include: "faidx.smk"
include: "genome_idx2ron.smk"
include: "fastq-sort.smk"
include: "removeDuplicatesAtFastq.smk"
include: "adaptor_handling.smk"

include: "length_distribution.smk"
include: "plot_length.smk"

include: "nucleotide_table.smk"


include: "align.smk"
include: "bed2geneCounts_tsnts.smk"
include: "report_stats.smk"
include: "check_presence.smk"
# include: "bedGraphToBigWig.smk"
# include: "simulation.smk"
# include: "nucleotide_table_sim.smk"
# include: "bam_correlation.smk"
# include: "bam_corr_graphs.smk"
include: "prepareOperonBed.smk"
include: "computeMatrix.smk"
# include: "plotProfile.smk"




################################################################################
