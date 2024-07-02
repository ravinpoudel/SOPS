This SOP is a part of [HLA-Evolutionary-Divergence ~ HED project](https://gitlab.com/neoantigen1/hla-evolutionary-divergence), in which we aimed to study the contribution of Microbiome Derived Epitopes (MDEs) on HLA divergence. Our pipeline can be divided into  two main parts:
- Analyses of shotgun metagenomics data and selection of genomes based on abundance (i.e. top 10). `snakefile_metaphlan.py` is the main script for this step.

- Predicting the MDEs for given HLA allotypes per patient. HLA allotypes were predicted using whole-exome sequence using a tool called OptiType. `snakefile_netmhcpan.py` is the main script for this step.

The full pipeline has the following 2 main steps:
1. NGS -> Composition -> Genome Selection :white_check_mark:
   - Metaphlan3
   - ncbi-genome-download
2. Selected Genome -> Epitopes Binding :white_check_mark:
   - netmhcpan

# Methods to access HED SOP:
Scripts and snakemake file associated with HED SOP can be downloaded from gitlab repo.


# 1. NGS -> Composition -> Genome Selection :white_check_mark:
``` bash

# clone the repo with code and env files

git clone https://gitlab.com/neoantigen1/hed_sop.git
cd hed_sop

# Create a conda environment and install the base tools. This has snakemake install.
conda env create --file envs/hed_base_env.yml -n hed

# Activate conda env
conda activate hed

```
_**Note: If needed, rule specific conda env file (env.yml) is used. `envs/` folder has all the `.yml` files used in the pipeline.**_

Preparation steps:
- Copy NGS data (Forward and Reverse) reads into a folder: `/data/ngs_data/`
- Prepare the metadata according to the provided template. Three metadata are needed:
- [ ] `sequence_info.csv`: **Required column**- First column is `SeqID`, Second column is `R1`, and Third column is `R2`
- [ ] `sample_df.csv`: **Required column**- `SeqID` and `WES_ID`
- [ ] `hla_metadata.csv` **Required column**- `WES_ID` and columns with HLA allotypes per WES_ID. This file is needed to link the bacterial compositional information(taxon) with HLA types of a patient.



Before we run the metaphlan3, we need to prepare the config file, which is in the json format. This can be created using `prepare_config_for_metaphlan.py` script. The script requires multiple inputs:
 - --fastq_path : is the folder where input files are located. Users have to put these files inside a folder: /data/ngs_data, both R1 (forward) and R2 (reverse) files are needed.
- --out_dir: this is the output folder, where all the outputs are stored. E.g. `metaphlan_out`

- --hla_metadata: A CSV file with information about the sample HLA type/allotypes. Columns: **WES_ID,A1,A2,B1,B2,C1,C2**. **Required Columns:** WES_ID and one of the columns for allotypes. WES_ID refers to the ID for Whole Exome Sequencing. WES_ID could remain same for different SeqID (Sequence ID), especially for longitudinal data.

- --sample_metadata: A CSV file with metadata. **Required Columns:** SeqID,and WES_ID. WES_ID refers to the ID for Whole Exome Sequencing. WES_ID could remain same for different SeqID (Sequence ID), especially for longitudinal data.

- --sequence_info: CSV file with a list of SeqID. The first column is SeqID. The second column is the file name of R1. The third column is the file name of R2.

- --metaphlanDB_path: Path to the metaphlan databases i.e chocophlan data (use: /share/projects/Neoantigen/Bioinformatics_tools/data/mpa_latest/ in Chinese HPC)

- --n_genomes: Number of genomes to download based on the composition analyses. For example, if 10 then top 10 genomes based on the abundance are selected

- --config_json: This is the name of the config file. E.g. `metaphlan_config.json`


**Script to run:**

```python
python script/prepare_config_for_metaphlan.py \
--fastq_path data/ngs_data \
--out_dir metaphlan_out \
--config_json metaphlan_config.json \
--hla_metadata data/metadata/hla_metadata.csv \
--sample_metadata data/metadata/sample_df.csv \
--sequence_info data/metadata/sequence_info.csv \
--relative_abundance_threshold 1 \
--updated_metadata_with_filepath updated_sample_df.csv \
--metaphlanDB_path /fsx/Bioinformatics_tools/data/mpa_latest

```



- Output: ` metaphlan_config.json`

```json
{
    "fastq_path": "data/ngs_data",
    "out_dir": "metaphlan_out",
    "hla_metadata": "data/metadata/hla_metadata.csv",
    "sample_metadata": "data/metadata/sample_df.csv",
    "sequence_info": "data/metadata/sequence_info.csv",
    "metaphlanDB_path": "/fsx/Bioinformatics_tools/data/mpa_latest",
    "relative_abundance_threshold": "1",
    "updated_metadata_with_filepath": "updated_sample_df.csv",
    "seqids": [
        "temp_715",
        "temp_1593"
    ],
    "temp_715": {
        "R1": "temp_715.R1.fastq.gz",
        "R2": "temp_715.R2.fastq.gz"
    },
    "temp_1593": {
        "R1": "temp_1593.R1.fastq.gz",
        "R2": "temp_1593.R2.fastq.gz"
    }
}

```

Now, we can run `snakefile_metaphlan.py` file.

```python


snakemake -s snakefile_metaphlan.py \
--configfile metaphlan_config.json \
--cores 8 \
--conda-frontend conda \
--use-conda


```

Output:
- `metaphlan_out/top_n_with_hla.csv` -- file with top n genomes per patient/SeqID and HLA metadata
- `metaphlan_out/acc_filepath_df.csv` -- filewith metadata and genome-path for the downloaded genomes
- `refseq` - folder with downloaded genomes

**NOTE:** While running snakemake, parameters defined in the config file can be over written from terminal. For instance, if we want to over write `n_genomes`:

```python

snakemake -s snakefile_metaphlan.py \
--configfile metaphlan_config.json \
--cores 8 \
--conda-frontend conda \
--use-conda \
--config n_genomes=5

```


# 2. Selected Genome -> Epitopes Binding :white_check_mark:

Now, we have finished running metaphlan and selected the genomes based on abundance, out next step in HED analyses is to predict the binding affinity of MDEs for a given patient HLA type. For this, as in the first step we need to prepare the config file (in json format). We will use the script called `prepare_config_for_netmhcpan.py`. This script requires:
 
Required inputs for the script are:
- [ ] -df == `acc_filepath_df.csv`: This is an output from the first step. Includes information about the HLA type and genome location.
- [ ] -idx == Single or multiple columns that we want to use as an index. Records in the index column are unique. For instance, in this case we want to run netmhcpan at unique combination of taxon and HLA, we are using these two columns as index columns. In the config file, information from these two columns will be used a key.
- [ ] -t == Number of threads to use.
- [ ] --config_json == Name of the config file. E.g. `netmhcpan_config.json`
- [ ] --humanEpitopesDB == Path to the folder with human epitopes. In AWS it it located at:`/fsx/Bioinformatics_tools/data/human_self_epitopes_db/`


**Running the script:**

```python

python script/prepare_config_for_netmhcpan.py \
-df metaphlan_out/updated_sample_df.csv \
-idx HLA taxon \
-t 8 \
--out_dir netmhcpan_out \
--config_json netmhcpan_config.json \
--humanEpitopesDB /fsx/Bioinformatics_tools/data/human_self_epitopes_db/


```


Output:
```json

{
    "out_dir": "netmhcpan_out",
    "humanEpitopesDB": "/fsx/Bioinformatics_tools/data/human_self_epitopes_db/",
    "indexid": [
        "HLA-C06:02_Alistipes_putredinis"
    ],
    "HLA-C06:02_Alistipes_putredinis": {
        "Allele": "A1",
        "rel_abund": 20.90179,
        "SeqID": 1593,
        "taxon": "Alistipes_putredinis",
        "WES_ID": "P100",
        "HLA": "HLA-C06:02",
        "f_tax": "Alistipes_putredinis",
        "file": "GCF_019132025.1_ASM1913202v1_protein.faa.gz",
        "Accession": "GCF_019132025.1",
        "filepath": "/fsx/ravin/neoantigen/hed_sop/refseq/bacteria/GCF_019132025.1/GCF_019132025.1_ASM1913202v1_protein.faa.gz"
    },
}


```

Now, using the config file we can run netmhcpan.

```python 

snakemake -s snakefile_netmhcpan.py \
--configfile netmhcpan_config.json \
--cores 8 \
--conda-frontend conda \
--use-conda

```

Output:
- `HLA-A02:01_Alistipes_putredinis_netmhcpan.log`
- `HLA-A02:01_Alistipes_putredinis_result_simplify.csv`
- `HLA-A02:01_Alistipes_putredinis_result_summary.txt.gz`
- `HLA-A02:01_Alistipes_putredinis_result.txt.gz`





