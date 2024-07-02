## Release: v0.4

This SOP is a part of **Neoantigen project** aims to understand the bacterial immuno-oncology potential (BIOP) using microbiome-derived epitomics (MDE). This pipeline can be divided into  three main parts:

- Analyses of shotgun metagenomics data and selection of genomes based on relative abundance. Based on our benchmarking analyses, we found selecting genomes at 0.01% relative abundance is computationally and biologically optimal. Thus, this SOP uses 0.01% relative abundance threshold to select the genomes. `snakefile_metaphlan_for_pangenome.py` is the main script for this step.

- Constructing pangenome for the selected taxon from compositional analyses. `snakefile_pangenome_with_filter.py` is the script for this step. In this step, we are downloading 10 different accessions, if available, for each taxon to construct a pangenome. If the total number of available accessions is less than 10, then all available accessions are considered. In addition, we also apply a gene prevalence-based filter to select the genes that are present in at least 50% of the selected accessions for constructing pangenomes. This approach will cover all the core genes, plus the accessory genes that are common in at least 50% of the accessions. Users have the option to define those parameters in the config file if wish to change or optimize for their study systems.

- Predicting the MDEs for given HLA allotypes per patient. [HLA allotypes were predicted using a whole-exome sequence using a tool called OptiType. We have a SOP to perform the HLA typing available at our sop WIKI](https://gitlab.xbiome.com/neoantigen1/neoantigen/-/wikis/SOP/Inhouse_CRAgs/Neoantigen-calling-on-Chinese-HPC). For runing netmhcpan, `snakefile_netmhcpan.py` is the main script.

The full pipeline has the following 3 main steps:
1. NGS -> Composition -> Genome Selection :white_check_mark:
   - Metaphlan3
   - ncbi-genome-download
2. Selected Genome -> Pangenomes :white_check_mark:
   - Prokka
   - Roary
3. Pangenomes -> Epitopes Binding :white_check_mark:
	- netmhcpan
	

## Quickstart using the input files that have been prepared for test samples

### Step 0: Setting up directory

- Create the following directory structure in the folder `/fsx/{USER}/biop_pangenome_test`

```plaintext
├── code
├── log
├── input_data
└── output
└── README
```

- `code`: Location where we want to keep all the functions, running codes, as well as envs
- `input_data`: Location where we keep all of our ngs_data and metadata
- output: Location to output files during analyses
- `log`: Folder to save the logs generated during the analyses


### Step 0.1: Cloning code and envs

- In the directory `/fsx/{USER}/biop_pangenome_test/code`, git clone the following:

```plaintext
git clone https://gitlab.xbiome.com/neoantigen1/BIOP_PANGENOME_SOP.git

```
After cloning, we see a folder called `BIOP_PANGENOME_SOP` inside which there are multiple snakemake files as well addition two folders:

```
└── BIOP_PANGENOME_SOP/
    ├── README.md
    ├── envs/
    ├── script/
    ├── snakefile_metaphlan_for_pangenome.py
    ├── snakefile_pangenome_with_filter.py
    ├── snakefile_prokka.py
    └── snakefile_update_metadata.py

```

`envs`: This contains yml files describing various programs needed to run the jobs

`script`: This contains all the scripts that are called by the snakemake file

`snakefile` is the one that we will be using to run the pipeline. 


### Step 0.2. Install a conda environment:

- The conda yml file is under the folder `code/envs` 

```
# Create a conda environment and install the base tools. This has snakemake installation as well.
conda env create --file code/BIOP_PANGENOME_SOP/envs/snakemake.yml -n pangenome

# Activate conda env
conda activate pangenome
```


## Getting started
For quick start testing, please copy over the following sample data to your folder: `/fsx/{USER}/biop_pangenome_test/input_data`
  
```plaintext
aws s3 cp s3://xbiome-us-aws/bkp_cur/ravin/BIOP_PANGENOME_SOP/script/misc/input_data/ input_data/ --recursive

```

Once metadata is copied, you can see the following nested folder structures with `metadata` folder:

```plaintext
/fsx/{USER}/biop_pangenome_test/input_data/
├── metadata
│   ├── hla_metadata.csv
│   ├── sample_df.csv
│   └── sequence_info.csv
└── ngs_data
    └── test_samples
        ├── temp_1593.R1.fastq.gz
        ├── temp_1593.R2.fastq.gz
        ├── temp_715.R1.fastq.gz
        └── temp_715.R2.fastq.gz

```

## Structure of data and metadata
Before we run the actual pipeline, SOP requires to structure data and metadata in a specific format. Inside `input_data` folder there are two sub-folders: `metadata` and `ngs_data` for storing metadata and MGS data respectively. For metadata, we need three files:
- `hla_metadata.csv`
	This file contains information about HLA typing for each WES_ID. `WES_ID` is the required column. Other columns in the file could be A1,A2,B1,B2,C1,C2 for storing HLA typing information.
- `sample_df.csv`
	This file contains metadata for each sample. `SeqID` and `WES_ID` are the required columns. Other columns could represent other metadata.
- `sequence_info.csv`
	This file contains information regarding the nsg data. Required columns are: `SeqID`, `R1` , and `R2`. R1 for forward fastq and R2 for reverse fastq. NOTE: First line starts with `#`.

For storing ngs data, you can add a folder inside within `metadata/data/ngs_data/`. 



# Step 1. NGS -> Composition -> Genome Selection :white_check_mark:
Before we run the metaphlan3, we need to prepare the config file, which is in the JSON format. This can be created using prepare_config_for_metaphlan.py script. The script requires multiple inputs:

`--fastq_path`: path to where the fastq/ngs data

`--out_dir_metaphlan`: path to the output directory for metaphlan

`--ncbi_genomes_download_dir`: path to the directory to download genomes from ncbi

`--hla_metadata`: hla metadata

`--sample_metadata`: sample metadata

`--sequence_info`: file with sequence information: SeqID, R1 , and R2

`--metaphlanDB_path`: path to the metaphlan database. For AWS: /fsx/Bioinformatics_tools/data/mpa_latest/

`--relative_abundance_threshold`: Relative abundance threshold to select genomes. Genomes with RA >= relative_abundance_threshold will be selected per sample,Recommendation: 0.01 percent

`--maximum_number_of_accessions`: Maximum number of accessions to be considered for pangenome if possible

`--log_genome_download`: Log file for tracking genome download

`--updated_metadata_with_filepath`:  Name for the updated_metadata_with_filepath, .csv. This file incorporates the path for the downloaded genomes.

`--config_json`: Name of the json file

Now, Let's create a config file for metaphlan.

```python
python code/BIOP_PANGENOME_SOP/script/prepare_config_for_metaphlan_for_pangenome.py \
--fastq_path input_data/ngs_data/test_samples/ \
--out_dir_metaphlan output/test_metaphlan_out \
--ncbi_genomes_download_dir output/test_ncbi_genome_download \
--hla_metadata input_data/metadata/hla_metadata.csv \
--sample_metadata input_data/metadata/sample_df.csv \
--sequence_info input_data/metadata/sequence_info.csv \
--metaphlanDB_path /fsx/shared/neoantigen_grp/metadata/mpa_latest/ \
--relative_abundance_threshold 0.01 \
--maximum_number_of_accessions 10 \
--log_genome_download log/test_genome_download.log \
--updated_metadata_with_filepath updated_sample_df.csv \
--config_json output/test_sample.json
```

Once we have the config file, we can use the config file to pass all the parameters for running metaphlan3. `snakefile_metaphlan_for_pangenome.py` does:

:one:  Run metaphlan3 for each sample

:two:  Concatenate results from each sample to a single file. Two files are generated as output: **species_known_taxa.csv**, **species_known_abundance.csv**

:three: Select taxon with relative abundance greater than 0.01 perent

:four:  Download accession for the selected taxon. The number of accessions downloaded for each taxon is defined by `maximum_number_of_accessions`


```bash
# Run metaphlan
#snakemake -s code/BIOP_PANGENOME_SOP/snakefile_metaphlan_for_pangenome.py --configfile output/test_sample.json --cores 8 --conda-frontend conda --use-conda -n

rm nohup.out; nohup snakemake -s code/BIOP_PANGENOME_SOP/snakefile_metaphlan_for_pangenome.py --configfile output/test_sample.json --jobs 10 --conda-frontend conda --use-conda --cluster "sbatch -D $PWD/log/ -n 8" 1>log/metaphlan_out.log &
```


# Step 2. Selected Genome -> Pangenomes :white_check_mark:

```bash
## Run prokka
# snakemake -s code/BIOP_PANGENOME_SOP/snakefile_prokka.py --configfile output/test_sample.json --cores 8 --conda-frontend conda --use-conda
rm nohup.out; nohup snakemake -s code/BIOP_PANGENOME_SOP/snakefile_prokka.py --configfile output/test_sample.json --jobs 10 --conda-frontend conda --use-conda --cluster "sbatch -D $PWD/log/ -n 8" 1>log/prokka_out.log &
```

```python
# Prepare config for pangenome
python code/BIOP_PANGENOME_SOP/script/prepare_config_for_pangenome_v2.py \
--pangenome_out_dir test_pandir \
--gene_prevalence_threshold 50 \
--ncbi_genomes_download_dir output/test_ncbi_genome_download \
--out_dir_metaphlan output/test_metaphlan_out \
--metadata_file_before_pangenome output/test_metaphlan_out/updated_sample_df.csv \
--indexid f_tax \
--metadata_file_after_pangenome test_sample_updated_after_pangenome.csv \
--config_json output/test_pangenome_config.json
```
**What happens at prepare_config_for_pangenome?**
In this step, we are preparing a config file for running pangenomes. In our metaphlan3 step, we defined the maximum number of accessions to be considered for each taxon. However, not all the taxon will have a defined number of accessions. In such cases, the maximum possible number of accessions are downloaded. In some cases, we have only one accession per taxon, in which it doesn't make sense to run pangenomes. Any taxon with more than one accession is considered for pangenome constriction. In addition, we are applying `gene_prevalence_threshold` .i.e. prevalence percentage criteria for any genes to be considered in the pangenome. `prepare_config_for_pangenome` prepares config for the genomes to run pangenomes at the defined gene_prevalence_threshold.

```bash
## Run Pangenome steps
rm nohup.out; nohup snakemake -s code/BIOP_PANGENOME_SOP/snakefile_pangenome_with_filter.py --configfile output/test_pangenome_config.json --jobs 10 --conda-frontend conda --use-conda --cluster "sbatch -D $PWD/log/ -n 8" 1> log/pangenome_out.log &
```

**Updating metadata file:** We have been tracking the changes and results at each step in the metadata file. We need to include all the final results in the metadata file so that we will be able to pull the right file using the file path.

```bash
snakemake -s code/BIOP_PANGENOME_SOP/snakefile_update_metadata.py --configfile output/test_pangenome_config.json --cores 8 --conda-frontend conda --use-conda

```


# 3. Selected Genome -> Epitopes Binding :white_check_mark:

Now, we have finished running metaphlan and selected the genomes based on abundance, out next step in BIOP analyses is to predict the binding affinity of MDEs for a given patient HLA type.In addition to the binding prediction of epitopes, we also plan to check the quality of the predicted epitopes by matching with epitopes/kmers generated from four publicly available databases. For this, as in the first step we need to prepare the config file (in json format). We will use the script called `prepare_config_for_netmhcpan_with_Epitopes_quality.py`. This script requires:
 
Required inputs for the script are:
- [ ] -df == `test_sample_updated_after_pangenome.csv `: This is an output from the first step. Includes information about the HLA type and genome location.
- [ ] -idx == Single or multiple columns that we want to use as an index. Records in the index column are unique. For instance, in this case we want to run netmhcpan at unique combination of taxon and HLA, we are using these two columns as index columns. In the config file, information from these two columns will be used a key.
- [ ] -t == Number of threads to use.
- [ ] --config_json == Name of the config file. E.g. `netmhcpan_config.json`
- [ ] --humanEpitopesDB == Path to the folder with human epitopes. In AWS it it located at:`/fsx/Bioinformatics_tools/data/human_self_epitopes_db/`
- [ ] vfdb_all_fas == path to vfdb_all_fasta file. This is a Virulence Factor Databases. http://www.mgc.ac.cn/VFs/
- [ ] --vfdb_curated_fas  == path to vfdb_curated_fasta file
- [ ] --card_fas == path to card file. This is a Comprehensive Antibiotic Resistance Database.https://card.mcmaster.ca
- [ ] --iedb_fas == path to iedb downloaded epitopes file. https://www.iedb.org
- [ ] --kalaora_fas ==path to Kalaora et al. fasta. This dataset was obtained from Kalaora et al. 2021 paper https://www.nature.com/articles/s41586-021-03368-8


**Running the script:**

```
python code/BIOP_PANGENOME_SOP/script/prepare_config_for_netmhcpan_with_Epitopes_quality.py \
-df output/test_metaphlan_out/test_sample_updated_after_pangenome.csv \
-idx HLA taxon \
-t 8 \
--out_dir output/test_netmhcpan_out \
--config_json output/netmhcpan_config.json \
--humanEpitopesDB /fsx/Bioinformatics_tools/data/human_self_epitopes_db \
--vfdb_all_fas /fsx/shared/neoantigen_grp/metadata/vfdb/VFDB_full_setB_pro.fas \
--vfdb_curated_fas /fsx/shared/neoantigen_grp/metadata/vfdb/VFDB_core_setA_pro.fas \
--card_fas /fsx/shared/neoantigen_grp/metadata/card/card_protein_fasta.fasta \
--iedb_fas /fsx/shared/neoantigen_grp/metadata/iedb/epitope_table_export_1640023078.fasta \
--kalaora_fas /fsx/shared/neoantigen_grp/metadata/epitopes_from_literature/Kalaora_et_al_Nat2021_bacteria_HLA_peptides_filled.fasta
```


Now, using the config file we can run netmhcpan.

```python 

#snakemake -s code/BIOP_PANGENOME_SOP/snakefile_netmhcpan_without_selfepitope_with_epitopes_quality.py --configfile output/netmhcpan_config.json --cores 8 --conda-frontend conda --use-conda -n

rm nohup.out; nohup snakemake -s code/BIOP_PANGENOME_SOP/snakefile_netmhcpan_without_selfepitope_with_epitopes_quality.py --configfile output/netmhcpan_config.json --jobs 200 --conda-frontend conda --use-conda --use-conda --cluster "sbatch -D $PWD/log/" 1>log/netmhcapan_out.log &

```

`output/test_netmhcpan_out/`:

- ` HLA-A01:01_Alistipes_putredinis_result_simplify.csv`
- `HLA-A01:01_Alistipes_putredinis_netmhcpan.log`
- `HLA-A01:01_Alistipes_putredinis_result.txt.gz`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers.txt`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_card`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_card.N_summary`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_iedb`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_iedb.N_summary`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_kalaora`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_kalaora.N_summary`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_vfdb_core`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_vfdb_core.N_summary`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_vfdb_full`
- `HLA-A01:01_Alistipes_putredinis_result_simplify_unique_8_mers_finalMDEs_vs_vfdb_full.N_summary`









