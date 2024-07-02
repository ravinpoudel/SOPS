## SOP for predicting MDEs using denovo genome assembly approach

**Date: 08/17/
2023** |  _**Author: Ravin Poudel**_  |   _**Email: ravinp@xbiome.us**_ | **Version**: [**denovo MGS SOP: 05/05/2023- v1.0**](https://gitlab.xbiome.com/neoantigen/denovo_mgs_sop/-/releases/v1.0) | [**GitLab**](https://gitlab.xbiome.com/neoantigen/denovo_mgs_sop)



This SOP is a part of **Neoantigen project** aims to understand the bacterial immuno-oncology potential (BIOP) using microbiome-derived epitomics (MDE). The primary focus of this analyses is to predict the MDEs from denovo assembled genomes, which are obtained from MGS data. This pipeline can be divided into  three main parts:



- Assembled the bacterial genomes using MGS data.`snakefile_denovo_to_binning_optimized.py` is the main script to run this step.

- Binning the assembled genomes, and refining the genome bins. We used `metabat2` and `maxbin2` for binning, while `bin refiner` for refining the bins. `snakefile_BinRefine_checkM_Quant_Taxonomy.py` is the script for running these analyses. Once the we obtained the refined bins, we predict the genes and translate the genes into protein sequence  using `prodigal`.

- Predicting the MDEs for given HLA allotypes per patient. [HLA allotypes were predicted using a whole-exome sequence using a tool called OptiType. We have a SOP to perform the HLA typing available at our sop WIKI](https://gitlab.xbiome.com/neoantigen/neoantigen/-/wikis/SOP/Inhouse_CRAgs/Neoantigen-calling-on-Chinese-HPC). For running netmhcpan, `snakefile_netmhcpan_for_denovo_with_human_background_and_TCR_recognization_potential.py` is the main script.

The full pipeline has the following 3 main steps:
1. MGS -> MAGS :white_check_mark:
   - megahit

2. MAGs -> Refined Bins -> Bins selection -> Bins Quantification -> Bins Taxonomy :white_check_mark:
   - metabat2
   - maxbin2
   - bin refiner
   - prodigal
   - checkM
   - GTDBtk
   
3. MAGs protein sequence -> MDEs -> Remove human background :white_check_mark:
	- netMHCpan
	- PEPMatch
	

## Quickstart using the input files that have been prepared for test samples

### Step 0: Setting up directory

- Create the following directory structure in the folder `/fsx/{USER}/denovo_mgs_sop_test`

```plaintext
├── code
├── log
├── input_data
├──────metadata: Copy from S3 bucket
├──────ngs_data: Directly read from S3 bucket
└── output
└── README

```

- `code`: Location where we want to keep all the functions, running codes, as well as envs

- `input_data`: Location where we keep all of our ngs_data and metadata. We have incorporated direct reading of ngs data from S3 bucket.For this, user need to setup/provide AWS credentials. This can be saved as csv file and path to the file can be used. Following is the format of my credential file/csv. 

```
Access key ID,Secret access key
XXXXXXXXXXXXXX,YYYYYYYYYYYYYYYY
,
us-east-1,
```

- output: Location to output files during analyses
- `log`: Folder to save the logs generated during the analyses


### Step 1: Cloning code and envs

- In the directory `/fsx/{USER}/denovo_mgs_sop_test/code`, git clone the following:

```plaintext
git clone https://gitlab.xbiome.com/neoantigen/denovo_mgs_sop.git

```
After cloning, we see a folder called `denovo_mgs_sop` inside which there are multiple snakemake files as well addition two folders:
- scripts- contains additional code that are called in the snakemake files
- envs - contains yaml files with conda envs information


Now copy metadata to run the test example:

```
aws s3 cp s3://xbiome-us-aws/finalized_projects/ravin/denovo_mgs_sop/input_data/metadata input_data/metadata/
```

### Step 2: Create a config file

```

python code/denovo_mgs_sop/script/prepare_config_for_denovo_optimized_for_tax_quant.py \
--fastq_path_S3 s3://xbiome-us-aws/finalized_projects/ravin/denovo_mgs_sop/input_data/ngs_data/test_samples/ \
--path_ngs_input_local /fsx/ravin/denovo_mgs_sop/input_data/ngs_data/test_samples/ \
--S3_key_file /fsx/ravin/.aws/Ravin_AWS_US_credential.csv \
--out_dir_denovo output/test_denovo \
--outdir_S3 s3://xbiome-us-aws/bkp_cur/ravin/denovo_mgs_sop/output/test_denovo \
--hla_metadata input_data/metadata/hla_metadata.csv \
--sample_metadata input_data/metadata/sample_df.csv \
--sequence_info input_data/metadata/sequence_info.csv \
--checkm_data /fsx/shared/neoantigen_grp/input_data/checkM_data/ \
--GTDBTK_DATA_PATH  /fsx/shared/neoantigen_grp/input_data/gtdbtk_r207/db \
--config_json output/test_sample_config.json


```

`--fastq_path`: path to where the fastq/ngs data. Here we are importing data directly from S3 bucket.

`--S3_key_file`: Credential for S3 and amazon account

`--out_dir_denovo`: Path to the output

`--outdir_S3`: Path to the folder in S3 bucker where we will be exporting the results

`--hla_metadata`: hla metadata

`--sample_metadata`: sample metadata

`--sequence_info`: file with sequence information: SeqID, R1 , and R2

`--checkm_data`: Reference data for running checkM

`--GTDBTK_DATA_PATH` : Reference database for running GTDBtk

`--config_json`: Config file

`Config` file contains all the information that we need to run the analyses. 

### Step 3: Run denovo and binning steps

Here we are assembling the microbial genomes per sampling using Megahit. Megahit is a tool that allow to assemble the genomes without the reference (denovo). Genomes are assembled by constructing De Bruijn graph with kmer sequences of various sizes. The final contigs from the megahit are then binned to multiple bins using metabat2 and maxbin2. Each bins represent a 
Metagenome-assembled genome (MAG), which is equivalent to a species. `snakefile_denovo_to_binning_optimized.py` is the main code to run the denovo assembly, the binning steps.


```
#  snakemake -s code/denovo_mgs_sop/snakefile_denovo_to_binning_optimized.py --configfile output/test_sample_config.json --conda-frontend conda --use-conda -n

nohup snakemake -s code/denovo_mgs_sop/snakefile_denovo_to_binning_optimized.py \
--configfile output/test_sample_config.json \
--conda-frontend conda --use-conda \
--jobs 17 \
--cluster "sbatch -D $PWD/log/ --cpus-per-task={threads}" 1> log/test_denovo_binning_nohup.out &

```

Multiple tools are available for predicting MAGs or bins from the large contigs. In order to improve bin/MAGs prediction, we used two binning tools: metabat2 and maxbin2. Bins predicted b these two tools are further refined using bin refiner tool. Bin refiner improves genome bins through the combination of different binning programs. Selected bins from the bin refiner were further access for genome completeness and contamination using checkM. Following the guide of [Minimum information about metagenome-assembled genome (MIMAG)](https://www.nature.com/articles/nbt.3893) we selected MAGs at `--length 50000 --completeness 50 --contamination 10`. The selected MAGs were quantity by mapping back the reads to the assembled genomes using salmon, and taxonomically annotated using GTDB.Genes and proteins sequences were predicted from the selected MAGs  using `prodigal`. Bin refining, checkM, gene/protein prediction can be made using `snakefile_BinRefine_checkM_Quant_Taxonomy.py` script. 


### Step 4: Bin refining and prodigal

```


nohup snakemake -s code/denovo_mgs_sop/snakefile_BinRefine_checkM_Quant_Taxonomy.py \
--configfile output/test_sample_config.json \
--conda-frontend conda --use-conda \
--jobs 14 \
--cluster "sbatch -D $PWD/log/ --cpus-per-task={threads}" 1> log/BinRefine_checkM_Quant_Taxonomy.out &



```

### Step 5: Prepare config for netmhcpan

Now, we have finished running denovo assembly, binning, and gene/protein prediction , out next step is to predict the binding potential of epitopes generated from MAGs for a given patient HLA type.For this, as in the first step we need to prepare the config file (in json format). We will use the script called `prepare_config_for_netmhcpan.py`. 


**Running the script:**


```

python code/denovo_mgs_sop/script/prepare_config_for_netmhcpan.py \
--metadata output/test_denovo/df_hla_meta_denovo_faa.csv \
--indexid WES_ID SeqID HLA faa_name \
--out_dir output/netmhcpan_output/test_sample \
--pepmatch_human_ref /fsx/shared/neoantigen_grp/Bioinformatics_tools/pepmatch_human_ref/human_GCF_000001405.40_protein.faa \
--iedb_positive_luksza_ref /fsx/shared/neoantigen_grp/Bioinformatics_tools/tcr_face_IEDB_postive_T_Cell_ref_Luksza/IEDB_positive_T-cell_assays.fasta \
--config_json output/test_sample_netmhcpan_config.json



```

This script requires:
 
Required inputs for the script are:
- [ ] --metadata == `df_hla_meta_denovo_faa.csv `: This is an output from the first step. Includes information about the HLA type and faa location.
- [ ] --indexid == Single or multiple columns that we want to use as an index. Records in the index column are unique. For instance, in this case we want to run netmhcpan at unique combination of taxon and HLA, we are using these two columns as index columns. In the config file, information from these two columns will be used a key.
- [ ] --config_json == Name of the config file. E.g. `test_sample_netmhcpan_config.json`
- [ ] --out_dir  == Output directory for netmhcpan
- [ ] --pepmatch_human_ref == Path to the human reference protein sequences.
- [ ] --iedb_positive_luksza_ref == Path to the iedb_positive_luksza_ref file. This is the positive peptide that are recognized by TCR based on experimental data, as provided by Luksza et., al. paper.



Now, using the config file we can run netMHCpan.

### Step 6: Running netMHCpan

```
nohup snakemake -s code/denovo_mgs_sop/snakefile_netmhcpan_for_denovo_with_human_background_and_TCR_recognization_potential.py \
--configfile output/test_sample_netmhcpan_config.json \
--conda-frontend conda --use-conda \
--jobs 10 \
--cluster "sbatch -D $PWD/log/ --cpus-per-task={threads}" 1> log/test_nohup_netmhc_run.out &

```



### Step 7: Copy results to the S3 bucket
In each of the snakeamake files we are flagging the intermediate files as `temp` so that once they are used by the pipeline, will get removed. The final `output` folder should mostly contain only the files that we want to save for the future or as a part of the final deliverable. We need to zip the `result` folder, move to S3, them remove from the local drive. In the config file `test_sample_config.json` we have already defined the location in S3 where we want to store our `result` folder. This can be done by using `snakefile_copy_to_S3_and_clean.py` script. 

```

nohup snakemake -s code/denovo_mgs_sop/snakefile_copy_to_S3_and_clean.py \
--configfile output/test_sample_config.json \
--conda-frontend conda --use-conda \
--jobs 1 \
--cluster "sbatch -D $PWD/log/ --cpus-per-task={threads}" 1> log/test_nohup_s3copy_clean_run.out &
```
Lets check the files by running following

```
aws s3 ls xbiome-us-aws/bkp_cur/ravin/denovo_mgs_sop/output/

```

Output:
- netmhcpan_output/
- test_denovo/
- test_sample_config.json
- test_sample_netmhcpan_config.json




