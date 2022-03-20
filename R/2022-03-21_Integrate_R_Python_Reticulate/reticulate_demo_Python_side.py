import pyranges as pr
import pandas as pd

# Use 'r.' preffix to retrieve objects from R session
grch38_gff = pr.read_gff3(r.gffPath)

# Convert into dataframe
df = grch38_gff.df

# For construction of string pattern to use with pd.str.contains()'s regex
suffix = "$"

search_string = "|".join([gene + suffix for gene in r.genes])

df = df[(df["gene"].str.contains(search_string)) & (df["Feature"] == "exon")]
