## BIOP-16S SOP
This SOP aims to study the contribution of Microbiome Derived Epitopes (MDEs) on HLA divergence using 16S rRNA gene sequence data. Our pipeline can be divided into  two main parts:
- Selection of genomes based on abundance (i.e. top 10). `snakefile_16S.py` is the main script for this step.

- Predicting the MDEs for given HLA allotypes per patient. HLA allotypes were predicted using whole-exome sequence using a tool called OptiType. `snakefile_netmhcpan.py` is the main script for this step.

The full pipeline has the following 2 main steps:
1. ASV (phyloseq Object)-> Genome Selection :white_check_mark:
   - Phyloseq
   - ncbi-genome-download
2. Selected Genome -> Epitopes Binding :white_check_mark:
   - netmhcpan

# Methods to access HED SOP for 16S data:
Scripts and snakemake file associated with HED SOP can be downloaded from gitlab repo.


# 1. ASV(Phyloseq) -> Genome Selection :white_check_mark:
``` bash

# clone the repo with code and env files

git clone https://gitlab.com/neoantigen1/hed_sop_16s.git
cd hed_sop_16s

# Create a conda environment and install the base tools. This has snakemake install.
conda env create --file envs/hed_16s_base_env.yml -n hed_16s

# Activate conda env
conda activate hed_16s

```
_**Note: If needed, rule specific conda env file (env.yml) is used. `envs/` folder has all the `.yml` files used in the pipeline.**_

Preparation steps:
- Copy `.rds_data_with_phyloseqobject` and `hla_metadata.csv`into a folder: `/data/metadata/`. These are the two required files.


**How to prepare `.rds_data_with_phyloseqobject` ?**

 - In addition to the row names, additional column called “SeqID” must be created in the otu_table from phyloseq object. Phyloseq requires sample_id or id as row names. This SOP requires a column named “SeqID” in the otu_table. Provided script does this for you. Basically, it will use the rownames of out_table to create a new column “SeqID”.

- Metadata requires two columns: `SeqID`, and `WES_ID`. `WES_ID` will be used to link `hla_metadata` with `metadata`. If not present in the ` .rds_data_with_phyloseqobject`, these information/columns need to be added in the metadata prior creating a phyloseq object.

```
# Example of adding a new column on metadata in the phyloseq object
library(phyloseq)
phylo_16s <- readRDS("phyloseq.rds")
o_data <- otu_table(phylo_16s)
t_data <- tax_table(phylo_16s)
m_data <- data.frame(sample_data(phylo_16s))

# Export it as csv and add "WES_ID" if needed
write.csv(m_data, "m_data.csv")

## Once WES_ID added, read in and create a new phyloseq object
new_m_data <- sample_data(read.csv("m_data.csv", header = TRUE))


# create a phyloseq object
tiny_phylo <- phyloseq(o_data, t_data, new_m_data)

# Save phyloseq as rds object
write_rds(tiny_phylo, "tiny_asv.rds")


```

**NOTE:** While selecting the genomes, first selection is done at the kingdom level, where Kingdom != "Eukaryota" and Genus != "NA" are selected. In this way, we will get bacteria and archaea genomes, and only the taxa that have information available at the genus level will be selected. Similar approach has been used in the recent paper by [Kalaora S et al., Nature 2021](https://www.nature.com/articles/s41586-021-03368-8). The authors have used 16S rRNA gene sequencing to select the bacterial genomes, then the protein sequence from the the selected "species" was downloaded from UniProt. Same tumor samples were sequenced using WGS, and the results were compared. Results suggest 16S rRNA gene sequencing as cost effective method.

- [ ] `hla_metadata.csv` **Required column**- `WES_ID` and columns with HLA allotypes per WES_ID. 



Before we run the selection steps, we need to prepare the config file, which is in the json format. This can be created using `prepare_config_for_16S.py` script. The script requires multiple inputs:

- --out_dir: this is the output folder, where all the outputs are stored. E.g. `16s_out`

- --hla_metadata: A CSV file with information about the sample HLA type/allotypes. Columns: **WES_ID,A1,A2,B1,B2,C1,C2**. **Required Columns:** WES_ID and one of the columns for allotypes. WES_ID refers to the ID for Whole Exome Sequencing. WES_ID could remain same for different SeqID (Sequence ID), especially for longitudinal data.

- --phyloseq_rds: rds data with phyloseq object

- --n_genomes: Number of genomes to download based on the composition analyses. For example, if 10 then top 10 genomes based on the abundance are selected

- --config_json: This is the name of the config file. E.g. `16s_config.json`


**Script to run:**

```python
python script/prepare_config_for_16S.py \
--out_dir 16s_out \
--hla_metadata data/metadata/hla_metadata.csv \
--phyloseq_rds data/metadata/tiny_asv.rds \
--n_genomes 5 \
--config_json 16s_config.json


```



- Output: ` 16s_config.json`

```json
{
    "out_dir": "16s_out",
    "hla_metadata": "data/metadata/hla_metadata.csv",
    "phyloseq_rds": "data/metadata/tiny_asv.rds",
    "n_genomes": "5"
}

```

Now, we can run `snakefile_16S.py` file.

```python


snakemake -s snakefile_16S.py \
--configfile 16s_config.json \
--cores 8 \
--conda-frontend conda \
--use-conda


```

Output:
- `16s_out/top_n_with_hla.csv` -- file with top n genomes per patient/SeqID and HLA metadata
- `16s_out/acc_filepath_df.csv` -- filewith metadata and genome-path for the downloaded genomes
- `refseq` - folder with downloaded genomes

**NOTE:** While running snakemake, parameters defined in the config file can be over written from terminal. For instance, if we want to over write `n_genomes`:

```python

snakemake -s snakefile_16S.py \
--configfile 16s_config.json \
--cores 8 \
--conda-frontend conda \
--use-conda \
--config n_genomes=10

```


# 2. Selected Genome -> Epitopes Binding :white_check_mark:

Now, we have finished selecting genomes based on abundance, our next step in HED analyses is to predict the binding affinity of MDEs for a given patient HLA type. For this, as in the first step, we need to prepare the config file (in json format). We will use the script called `prepare_config_for_netmhcpan.py`. This script requires:
 
Required inputs for the script are:
- [ ] -df == `acc_filepath_df.csv`: This is an output from the first step. Includes information about the HLA type and genome location.
- [ ] -idx == Single or multiple columns that we want to use as an index. Records in the index column are unique. For instance, in this case we want to run netmhcpan at unique combination of taxon and HLA, we are using these two columns as index columns. In the config file, information from these two columns will be used a key.
- [ ] -t == Number of threads to use.
- [ ] --config_json == Name of the config file. E.g. `netmhcpan_config.json`


**Running the script:**

```python

python script/prepare_config_for_netmhcpan.py \
-df 16s_out/acc_filepath_df.csv \
-idx HLA taxon \
-t 8 \
--out_dir netmhcpan_out \
--config_json netmhcpan_config.json


```


Output:
```json

{
    "out_dir": "netmhcpan_out",
    "indexid": [
       "HLA-A01:01_Bacteroides_plebeius",
    ],
 "HLA-A01:01_Bacteroides_plebeius": {
        "Accession": "GCF_014334015.1",
        "f_tax": "Bacteroides_plebeius",
        "SeqID": "seq_id9143",
        "filepath": "/fsx/ravin/hed_sop_16s/refseq/bacteria/GCF_014334015.1/GCF_014334015.1_ASM1433401v1_protein.faa.gz",
        "WES_ID": "P100",
        "HLA": "HLA-A01:01",
        "taxon": "Bacteroides_plebeius",
        "file": "GCF_014334015.1_ASM1433401v1_protein.faa.gz",
        "Allele": "A1",
        "rel_abund": 704
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




