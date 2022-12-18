# The Mfd Protein is the Transcription-Repair Coupling Factor (TRCF) in *Mycobacterium Smegmatis*
In vitro and in vivo experiments with *Escherichia coli* have shown that the Mfd translocase is responsible for transcription-coupled repair which is defined as the faster rate of repair of the transcribed strand than the non-transcribed strand by nucleotide excision repair. Even though the mfd gene is conserved in all bacterial lineages, there is only limited information on whether it performs the same function in other bacterial species. Here, by genome scale analysis of repair of UV-induced cyclobutene dimers we find that the Mfd protein is the Transcription-Repair Coupling Factor (TRCF) in *Mycobacterium smegmatis*. This finding, combined with the inverted strandedness of UV-induced mutations in wild-type and mfd- *Escherichia coli* and *Bacillus subtilis* indicate that the Mfd protein is the universal TRCF in bacteria.

This repository contains the data analysis workflow.

![Workflow](/results/figs/workflow.png "Workflow")

<br>

## Installation

- This workflow is prepared using 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management 
system and [conda](https://docs.conda.io/en/latest/)

- To run the workflow, you should have conda installed for environment 
management. All the other packages including Snakemake and their dependencies 
can be obtained automatically through environments prepared for each step of 
the workflow. You can follow the installation steps from 
[the link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).

- Initially, you should clone the repository and navigate into the directory: 

    ```
    git clone https://github.com/CompGenomeLab/MycoSmegmatis_TRCF.git
        
    cd MycoSmegmatis_TRCF
    ```

- Next, you should create a conda environment with the defined packages. 
Install [mamba](https://mamba.readthedocs.io/en/latest/) 
and create the environment using mamba:

    ```
    conda install -c conda-forge mamba

    mamba create -c bioconda -c conda-forge -c r -n repair snakemake=6.3.0 python=3.8 rust=1.50 sra-tools=2.11.0

    conda activate repair
    ```

<br>

## Directory Structure

This workflow is prepared according to the 
[structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 
recommended by Snakemake: 

- `config/`: contains the configuration files.

- `logs/`: contains the log files of each step. 
This folder will automatically appear when you run the workflow.

- `report/`: contains the description files of figures,
which will be used in reports.

- `resources/`: contains `samples/` where the raw XR-seq and Damage-seq data  are stored, 
`input/` where the input files are stored, 
and `ref_genomes/` where the reference genome files are stored. 
Reference genome files can be automatically produced by the workflows, 
if they are properly defined in the config files.  

- `results/`: contains the generated files and figures. 
This folder will automatically appear when you run the workflow.

- `workflow/`: contains `envs/` where the environments are stored, 
`rules/` where the Snakemake rules are stored, and 
`scripts/` where the scripts used inside the rules are stored. 

<br>

## Configuration file

The configuration file with "_initial_" prefix shouldn't be modified by the user since they are containing configuration settings that are common for all XR-seq experiments. 
For more detail about these configuration files, check out the readme file in `config/` directory. 
The parameters for "config.yaml" as below:

- `sample`: The name of the sample file w/o the extension. 
Multiple sample names can be given in the below format:

    ```
    sample: 
        - "SAMPLE_1"
        - "SAMPLE_2"
        - "SAMPLE_3"
    ```

    - Using the given sample name, the workflow will look for 
    `{SAMPLE}.fastq.gz` as raw data. 
    Therefore, the fastq file must be gzipped before running the workflow.


- `damage_type`: Damage type of each sample should be provided here in the 
same order of the samples:

    ```
    damage_type: 
        - "64"
        - "CPD"
        - "oxaliplatin"
    ```

    - Currently damages below are available can be provided as (case-insensitive):

        - (6-4)PP: `64`, `64pp`, `(6-4)pp`, `6-4pp`;
        - CPD: `CPD`;
        - Cisplatin: `cisplatin`;
        - Oxaliplatin: `oxaliplatin`.

<br>

## Usage

After adjusting the configuration file, you can run the workflow from this directory.

This workflow runs on Slurm Workload Manager](https://slurm.schedmd.com/srun.html)

    ```
    source al.sh
    rna 
    ```

After rna pipeline is completed you should run

    ```
    source al.sh
    smm 
    ```


<br>

To generate detailed HTML report files, 
the code below should be run after workflow:

```
snakemake --report report.zip
```