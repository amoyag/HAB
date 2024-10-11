#!/usr/bin/env python
# coding: utf-8

# ## Functional Analysis with STRING-dB

# [STRING](https://string-db.org) has an extensive API. We can perform statistical enrichment analysis as explained [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531103/).
# 
# The STRING enrichment API method allows you to retrieve functional enrichment for any set of input proteins. It will tell you which of your input proteins have an enriched term and the term's description. The API provides the raw p-values, as well as, False Discovery Rate and Bonferroni corrected p-values. 
# See the documentation of the STRING API [here](https://string-db.org/cgi/help.pl?subpage=api%23getting-protein-protein-interaction-enrichment)

# ### Example. Get significantly enriched (FDR < 1%) GO Processes for the given set of proteins

# #### Setup

# In[2]:


import requests ## python -m pip install requests 
import json


string_api_url = "https://version-11-5.string-db.org/api"
output_format = "json" # tsv tab separated values, with a header line; tsv-no-header tab separated values, without header line; json JSON format; xml XML format
method = "enrichment"


# #### Construct the request

# In[3]:


request_url = "/".join([string_api_url, output_format, method])


# #### Set parameters 
# Note that the genes IDs are STRING IDs. You might need to convert them from other IDs to STRING. See the documentation for [mapping identifiers.](https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers)
# 
# Also note that you need to specify the species NCBI identifier. Finally, note how the identifiers are joined to perform the query:
# 
# The Python code `"%0d".join(my_genes)` is using the `join` method to concatenate a list or iterable of strings, `my_genes`, with a specific string separator. However, the separator used here seems to be a bit unusual. 
# 
# Let's break it down:
# 
# - `join`: This is a string method in Python that is used to join the elements of an iterable (like a list or tuple) into a single string. The elements of the iterable must be strings.
# 
# - `"%0d"`: This is a string that contains a percent sign (%) followed by a 0 and a letter d. The percent sign is often used in string formatting, but in this context, it is treated as a regular character, and there is no string formatting happening.
# 
# So, if `my_genes` is a list of strings, this code will concatenate those strings using "%0d" as the separator between the elements. 
# 
# For example, if `my_genes = ["A", "T", "G", "C"]`, the result would be `"A%0dT%0dG%0dC"`.

# In[4]:


my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
            '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']


params = {

    "identifiers" : "%0d".join(my_genes), # your protein
    "species" : 7227, # species NCBI identifier 
    "caller_identity" : "test_HAB" # your app name

}


# #### Call STRING

# In[5]:


response = requests.post(request_url, data=params)


# #### Read and parse the results
# We have selected JSON as the output format, so we will process the enrichment data accordingly. Note that we need to define which category do we need to print out.

# In[8]:


data = json.loads(response.text)

for row in data:

    term = row["term"]
    preferred_names = ",".join(row["preferredNames"])
    fdr = float(row["fdr"])
    description = row["description"]
    category = row["category"]

    # if category == "Process" and fdr < 0.01: # GO BP
    # if category == "KEGG" and fdr < 0.01: # KEGG Pathways
    if fdr < 0.01: # All STRING Categories

        ## print significant GO Process annotations

        print("\t".join([term, preferred_names, str(fdr), category, description]))

