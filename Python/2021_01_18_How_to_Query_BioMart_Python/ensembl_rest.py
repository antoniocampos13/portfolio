import glob
import subprocess
from collections import defaultdict
from pathlib import Path

import dask.dataframe as dd
import pandas as pd
from dask import delayed

from src.files_to_pandas import files_to_pandas
from src.query_biomart import query_biomart

# Set up the file paths:
project_folder = Path().resolve() # project root
data_folder = project_folder / "data"
go_list = data_folder / "go_demo.xlsx"

# Load file
df = pd.read_excel(go_list, engine="openpyxl")

# Export query commands
with open(data_folder / "commands.sh", "w", newline="\n") as commands:
    
    commands.write("#!/bin/bash" + "\n")
    
    for go_id in df["go_ids"]:
        commands.write(query_biomart(filter_names=["go_parent_term"],
        values=[go_id], attribute_names=["external_gene_name"],
        uniqueRows=True) + "\n")

# Download queries results into data folder...
subprocess.call(str(data_folder / "commands.sh"), cwd=data_folder)

#... Then import data into a pandas.DataFrame with the help of dask
file_list = data_folder.glob("*.txt")

dfs = [delayed(files_to_pandas)(filename) for filename in file_list]

ddf = dd.from_delayed(dfs)

# Create a defauldict to create GO id --> genes relation
go_to_genes = defaultdict(list)

for go_id, gene in zip(ddf["filename"], ddf["gene_name"]):
    go_to_genes[go_id].append(gene)

go_genes = pd.DataFrame(go_to_genes.items(), columns=["go_ids", "gene_name"])

# Error checking. Network errors. "Query" is used as an example. Other error messages can be used.
# "Query ERROR: caught BioMart::Exception::Database: Could not connect to mysql database ..."
for go_id, gene in zip(go_genes["go_ids"],go_genes["gene_name"]):
    if "Query" in str(gene):
        print(f"Error in {go_id}")
    else:
        continue

print("Error check done.")

df_annotated = df.set_index("go_ids").join(go_genes.set_index("go_ids"))

df_annotated["gene_string"] = [", ".join(map(str, element)) for element in df_annotated["gene_name"]]

with pd.ExcelWriter(go_list, mode="a") as writer:
    df_annotated.to_excel(writer, sheet_name="go annotated", index=True)
