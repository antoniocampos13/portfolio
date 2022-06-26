import pandas as pd

TABLE = "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff"
EDITED_TABLE = "human_gene_exons_GRCh38.gff"
CHUNKSIZE = 20000000

column_names = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
]

chunks = pd.read_table(TABLE, names=column_names, chunksize=CHUNKSIZE)

# # Or:
# chunks = pd.read_csv(TABLE, names=column_names, chunksize=CHUNKSIZE, sep="\t")

for chunk in chunks:
    # Step 1
    temp_df = chunk[chunk["type"] == "exon"]

    # Step 2
    temp_df = temp_df.drop(["source", "score", "strand", "phase", "attributes"])

    # Step 3
    temp_df.to_csv(EDITED_TABLE, header=False, mode="a", sep="\t", index=False)
