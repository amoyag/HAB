### Introduction to Functional Analysis in Bioinformatics

Functional analysis in bioinformatics involves understanding the biological function of genes or proteins within an organism. After sequencing or identifying a gene of interest, researchers often need to determine its function and role in cellular processes. Functional analysis provides insights into gene ontology, biological pathways, molecular interactions, and overall cellular function.

One common approach for functional analysis is **overrepresentation analysis**, which aims to identify biological functions that are significantly enriched within a group of genes of interest compared to a reference group. This is done by comparing the annotated functions of the genes in the target group against those in the reference (or background) group.

A typical statistical method used for this purpose is the **Fisher exact test**. For each functional category, the proportion of genes in the target group that are annotated with a given function is compared to the proportion of genes in the background group with the same annotation. If the target group has a significantly higher proportion of genes with that function, the function is considered overrepresented.

The Fisher exact test can be represented using a contingency table. For example, consider a functional category (e.g., cell cycle regulation) and the number of genes in two lists: the target list (genes of interest) and the background list:

|                | In Category | Not in Category | Total |
|----------------|-------------|-----------------|-------|
| Genes of Interest | 3           | 9               | 12    |
| Background       | 4           | 20              | 24    |
| Total            | 7           | 29              | 36    |

Using the Fisher exact test, we calculate the p-value to determine if the proportion of genes of interest in the category is significantly higher than expected by chance. The null hypothesis is that the representation of the genes of interest in a functional category is equal to the representation of the background in that category.

The Fisher exact test formula for calculating the p-value is:


```math
p = \frac{ \binom{a+b}{a} \binom{c+d}{c} }{ \binom{a+b+c+d}{a+c} }
```


where:
- $a$ is the number of genes of interest in the category,
- $b$ is the number of genes of interest not in the category,
- $c$ is the number of background genes in the category,
- $d$ is the number of background genes not in the category.

The p-value represents the probability of observing this distribution or one more extreme under the null hypothesis.

The p-value obtained from the Fisher exact test can be very small, such as **$p$ = 1.13 × 10⁻²⁹**, indicating significant overrepresentation. Since multiple categories are often tested, **p-value correction** (e.g., using the Benjamini-Hochberg method) is applied to control the false discovery rate.

Gene Set Enrichment Analysis (GSEA) is a computational method used to determine whether a predefined set of genes shows statistically significant, concordant differences between two biological states (e.g., healthy vs. diseased). Unlike overrepresentation analysis, which focuses on discrete gene lists, GSEA considers all genes ranked by their expression levels, making it more robust against arbitrary thresholds for differential expression. It is particularly useful for capturing subtle but coordinated changes across a set of related genes, such as those involved in the same pathway or biological process. GSEA uses permutation testing to assess the statistical significance of the observed enrichment, providing insights into the functional context of experimental data.

With the vast amount of genomic data being generated, functional analysis plays a key role in making sense of raw sequences, converting them into biologically meaningful insights. Several Python modules and tools are available to facilitate functional analysis, allowing researchers to automate and simplify their workflows. Below are some Python packages commonly used for this purpose, along with example scripts illustrating how to use each of them.

### 1. **Biopython**

`Biopython` is a comprehensive toolkit for biological computation. It provides functions for sequence analysis, querying biological databases, and functional analysis, including annotations and domain searches. `Biopython` is not typically used for performing functional analysis in the same way as other specialized tools. It is primarily a toolkit for biological computation, providing modules for: sequence analysis (e.g., parsing FASTA/GenBank files, sequence alignment), database querying (e.g., accessing NCBI Entrez, retrieving protein/gene records) and annotation retrieval (e.g., extracting basic information about genes and proteins). Thus, Biopython allows you to **retrieve gene or protein data**, such as functional annotations, GO terms, and pathways. It can be useful as a complementary tool to **fetch data** (e.g., retrieving sequences, annotations, or even related pathways for genes of interest), which can then be analyzed with other specialized software for functional analysis.

```python
from Bio import Entrez
import time

# Use Entrez to search for genes and fetch their functional annotations
Entrez.email = "your_email@example.com"

# Define the list of glycolysis genes
gene_symbols = ["GAPDH", "PKM", "ALDOA", "HK1", "PFKM", "LDHA"]

# Iterate over each gene to fetch annotations
for gene in gene_symbols:
    try:
        # Search for the gene in Entrez
        handle = Entrez.esearch(db="gene", term=f"{gene}[Gene] AND Homo sapiens[Organism]")
        record = Entrez.read(handle)
        handle.close()

        if record['IdList']:
            gene_id = record['IdList'][0]

            # Fetch the gene summary to get functional annotations
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
            gene_record = handle.read()
            handle.close()

            # Print the fetched annotation
            print(f"Annotations for {gene}:\n")
            print(gene_record)
            print("\n" + "="*50 + "\n")
        else:
            print(f"No record found for gene: {gene}")

        # Pause to avoid overwhelming the Entrez servers
        time.sleep(1)

    except Exception as e:
        print(f"An error occurred while fetching annotations for {gene}: {e}")
```

### 2. **GOATOOLS**

`GOATOOLS` is a Python library that helps with Gene Ontology (GO) enrichment analysis. It allows researchers to perform statistical analysis on a set of genes to determine which biological processes, cellular components, or molecular functions are significantly enriched.

```python
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf
import random
import requests
import gzip
import shutil
import os
import mygene

### Download stuff

## Get the GO ontology obo file
# URL of the go-basic.obo file
obo_url = "http://current.geneontology.org/ontology/go-basic.obo"
obo_file_path = "go-basic.obo"

# Download only if the file does not already exist
if not os.path.exists(obo_file_path):
    try:
        response = requests.get(obo_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        with open(obo_file_path, 'wb') as file:
            file.write(response.content)
        print(f"Successfully downloaded 'go-basic.obo' to {obo_file_path}")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while downloading the file: {e}")
else:
    print(f"'{obo_file_path}' already exists. Skipping download.")

## Download the GAF file
# URL of the goa_human.gaf file
gaf_url = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
compressed_gaf_file_path = "goa_human.gaf.gz"
extracted_gaf_file_path = "goa_human.gaf"

# Download and extract only if the file does not already exist
if not os.path.exists(extracted_gaf_file_path):
    try:
        response = requests.get(gaf_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        with open(compressed_gaf_file_path, 'wb') as file:
            file.write(response.content)
        print(f"Successfully downloaded 'goa_human.gaf.gz' to {compressed_gaf_file_path}")

        # Extract the .gz file
        with gzip.open(compressed_gaf_file_path, 'rb') as f_in:
            with open(extracted_gaf_file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        print(f"Successfully extracted '{extracted_gaf_file_path}'")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while downloading the file: {e}")
else:
    print(f"'{extracted_gaf_file_path}' already exists. Skipping download.")

### Convert genes ID to UniProt ID
# Initialize MyGene client
mg = mygene.MyGeneInfo()

# Convert gene symbols to UniProt IDs
gene_symbols = ["GAPDH", "PKM", "ALDOA", "HK1", "PFKM", "LDHA"]
query_result = mg.querymany(gene_symbols, scopes='symbol', fields='uniprot', species='human')

# Extract UniProt IDs for the study genes
study_gene_uniprot_ids = []
for gene in query_result:
    uniprot_id = gene.get('uniprot', {}).get('Swiss-Prot')
    if uniprot_id:
        study_gene_uniprot_ids.append(uniprot_id)

# Print the UniProt IDs for each gene
print("Converted UniProt IDs for the study genes:")
for gene_symbol, uniprot_id in zip(gene_symbols, study_gene_uniprot_ids):
    print(f"Gene Symbol: {gene_symbol}, UniProt ID: {uniprot_id}")

### Set the background list and perform ORA

# Extract all genes from the GAF file
all_genes = set()

try:
    with open(extracted_gaf_file_path, 'r') as gaf_file:
        for line in gaf_file:
            # Skip comments
            if line.startswith("!"):
                continue

            # Split the GAF line by tab
            columns = line.strip().split("\t")

            # Column 2 contains the systematic identifier (e.g., UniProt ID)
            gene_identifier = columns[1]

            # Add the gene identifier to the set of all genes
            all_genes.add(gene_identifier)

    print(f"Successfully extracted {len(all_genes)} gene identifiers from the GAF file.")

except FileNotFoundError:
    print(f"File '{extracted_gaf_file_path}' not found. Please make sure it is available.")
except Exception as e:
    print(f"An error occurred while reading the GAF file: {e}")

# Convert to list for sampling
all_genes_list = list(all_genes)

# Randomly sample 1000 genes from the full list
random.seed(42)  # Set seed for reproducibility
background_gene_list = random.sample(all_genes_list, 1000)

# Ensure that all study genes are included in the background list
for gene in study_gene_uniprot_ids:
    if gene not in background_gene_list:
        background_gene_list.append(gene)

# Perform GO enrichment analysis with GOATOOLS using the updated study list

# Load the GO ontology
go = obo_parser.GODag(obo_file_path)

# Load gene-to-GO associations in GAF format
gene2go = read_gaf(extracted_gaf_file_path)

# Perform GO enrichment analysis
go_enrich = GOEnrichmentStudy(
    background_gene_list,  # Updated background with genes
    gene2go,               # Gene to GO term associations
    go,                    # Gene Ontology DAG
    propagate_counts=False,
    alpha=0.05,            # Significance cutoff
    methods=['fdr_bh']     # Multiple testing correction
)

# Run study with UniProt IDs from the study set
enriched_results = go_enrich.run_study(study_gene_uniprot_ids)

# Print the enriched GO terms
for res in enriched_results:
    if res.p_fdr_bh < 0.05:
        print(f"GO ID: {res.GO}, Description: {res.name}, P-value: {res.p_fdr_bh}")

```


### 3. **gseapy**

`gseapy` is a Python wrapper for Gene Set Enrichment Analysis (GSEA), a common technique used to determine if a predefined set of genes shows statistically significant differences between two biological states (e.g., healthy vs. diseased).

```python
import gseapy as gp
import pandas as pd

# Define the list of glycolysis genes and convert them to uppercase
glycolysis_genes = ["GAPDH", "PKM", "ALDOA", "HK1", "PFKM", "LDHA"]
glycolysis_genes_upper = [gene.upper() for gene in glycolysis_genes]

# Create a DataFrame to mimic expression data for GSEA analysis
# In practice, this would come from experimental results (e.g., RNA-seq)
# Here we create mock expression data for illustration purposes
mock_expression_data = {
    'Gene': glycolysis_genes_upper,
    'Expression': [2.3, 1.8, 2.1, 1.9, 2.5, 2.0]  # Random fold change values for example purposes
}

# Convert the mock data to a DataFrame
expression_df = pd.DataFrame(mock_expression_data)
expression_df.set_index('Gene', inplace=True)

# Print the mock expression data to see what it looks like
print("Mock Expression Data:")
print(expression_df)

# Save expression data to a CSV file to use for GSEA
expression_df.to_csv('expression_data.csv')

# Run GSEA Preranked analysis with gseapy
# Use a predefined pathway gene set like KEGG
try:
    gsea_results = gp.prerank(
        rnk='expression_data.csv',  # Ranked gene list, usually the output of a differential analysis
        gene_sets='KEGG_2016',      # Gene set database; using KEGG pathways for example
        outdir='gsea_output',       # Output directory
        permutation_num=100,        # Number of permutations for significance testing
        min_size=3,                 # Minimum size of the gene set to include
        max_size=5000,              # Maximum size of the gene set to include
    )

    # Print top 5 results
    print("\nTop 5 Enriched Pathways:")
    print(gsea_results.res2d.head())

except Exception as e:
    print(f"An error occurred during GSEA analysis: {e}")

```

### 4. **Enrichr**

`Enrichr` is a powerful web-based tool and API for performing gene set enrichment analysis. It helps researchers determine if a predefined set of genes is significantly overrepresented in a user-provided list, providing insights into biological processes, pathways, and functions related to the genes of interest. Enrichr supports a vast collection of gene set libraries, including KEGG pathways, Gene Ontology (GO), Reactome, and many more, making it highly versatile for functional analysis. The tool integrates advanced scoring methods, such as combined scores, which consider both statistical significance and magnitude of enrichment, making it especially useful for prioritizing biologically relevant results.

```python
import requests
import json

# Define the list of glycolysis genes and convert them to uppercase
glycolysis_genes = ["GAPDH", "PKM", "ALDOA", "HK1", "PFKM", "LDHA"]
glycolysis_genes_upper = [gene.upper() for gene in glycolysis_genes]

# Enrichr API URLs
add_list_url = "https://maayanlab.cloud/Enrichr/addList"
enrich_url = "https://maayanlab.cloud/Enrichr/enrich"

# Add genes to Enrichr
gene_str = "\n".join(glycolysis_genes_upper)  # Enrichr requires newline-separated gene symbols
payload = {
    'list': (None, gene_str),
    'description': (None, 'Glycolysis Genes List')
}

response = requests.post(add_list_url, files=payload)
if response.status_code == 200:
    result = response.json()
    user_list_id = result['userListId']
    print(f"Successfully added genes to Enrichr. User List ID: {user_list_id}")
else:
    print(f"Error adding genes to Enrichr: {response.status_code}")
    exit()

# Perform enrichment analysis
gene_set_library = 'KEGG_2019_Human'  # Example: KEGG pathways
params = {
    'userListId': user_list_id,
    'backgroundType': gene_set_library
}

response = requests.get(enrich_url, params=params)
if response.status_code == 200:
    enrichment_results = response.json()
    if gene_set_library in enrichment_results:
        print("\nEnrichment Analysis Results (Top Results):")
        for result in enrichment_results[gene_set_library][:5]:  # Display top 5 results
            # result is a list, so we access items by index
            term_name = result[1]  # Term name
            pvalue = result[2]     # P-value
            combined_score = result[4]  # Combined score

            print(f"Term: {term_name}, P-value: {pvalue}, Combined Score: {combined_score}")
    else:
        print("The expected gene set library was not found in the enrichment results.")
else:
    print(f"Error retrieving enrichment results: {response.status_code}")

```

### 5. **STRINGdb**

`STRINGdb` is a Python module that interacts with the STRING database, a resource for exploring known and predicted protein-protein interactions. It allows users to retrieve interaction data and perform enrichment analysis for a set of proteins.

```python
import requests
import json

# STRINGdb API URL and method details
string_api_url = "https://version-11-5.string-db.org/api"
output_format = "json"
method = "enrichment"

# Construct the request URL
request_url = "/".join([string_api_url, output_format, method])

# Glycolysis genes to use for functional enrichment analysis
glycolysis_genes = ["GAPDH", "PKM", "ALDOA", "HK1", "PFKM", "LDHA"]

# Define the parameters for the request
params = {
    "identifiers": "%0d".join(glycolysis_genes),  # Protein list formatted with %0d (newline separator)
    "species": 9606,  # Species NCBI identifier for Homo sapiens (human)
    "caller_identity": "test_HAB"  # Replace with your own identifier or email address
}

# Make the POST request to STRINGdb API
response = requests.post(request_url, data=params)

# Parse the JSON response
data = json.loads(response.text)

# Print the significant GO Biological Process (Process) annotations
print("Functional Enrichment Analysis Results (GO Biological Process Only, FDR < 0.01):")
for row in data:
    term = row["term"]
    preferred_names = ",".join(row["preferredNames"])
    fdr = float(row["fdr"])
    description = row["description"]
    category = row["category"]

    # Filter results to only include GO Biological Processes with FDR < 0.01
    if category == "Process" and fdr < 0.01:
        print("\t".join([term, preferred_names, str(fdr), category, description]))
```


