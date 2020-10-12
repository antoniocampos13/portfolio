import glob
import gzip
import ntpath
import shutil
from pathlib import Path

import dask.dataframe as dd
import pandas as pd
from dask import delayed
from sqlalchemy import create_engine

from settings import *

# 1. Set up PostgreSQL connection object
connection_string = (
    f"{DB_FLAVOR}+{DB_PYTHON_LIBRARY}://{USER}:{PASSWORD}@{DB_HOST}:{PORT}/{DB_NAME}"
)
engine = create_engine(connection_string)

# 2. Set up project paths
project_folder = Path().resolve()
src_folder = project_folder / "src"
data_folder = project_folder / "data"
output_folder = project_folder / "output"

counts_folder = data_folder / "counts"

# 3. Decompress files
compressed_files = counts_folder.glob("*.gz")


def uncompress(path):
    with gzip.open(str(path), "rb") as src, open(str(path).rstrip(".gz"), "wb") as dest:
        shutil.copyfileobj(src, dest)


for file in compressed_files:
    uncompress(file)

# 4. Getting all counts files
file_list = counts_folder.glob("*.counts")

# 5. Create a function ready to return a pandas.DataFrame
def read_and_label_csv(filename):
    """ Converts each file into a pandas.DataFrame """

    tmp = pd.read_csv(filename, delimiter="\t", names=["gene_id", "gene_count"])

    tmp["filename"] = ntpath.basename(filename)+".gz"
    return tmp

# 6. Create a list of commands applying the read_and_label_csv function to all files
dfs = [delayed(read_and_label_csv)(filename) for filename in file_list]

# 7. Using delayed, assemble the pandas.DataFrames into a dask.DataFrame
ddf = dd.from_delayed(dfs)

# 7.5 Optional: uncomment below to export the dask.DataFrame as HUGE CSV file. WARNING: IT WILL USE A LOT OF RAM AND CPU.
# ddf.to_csv(Path(output_folder, "gene_counts.csv"), single_file=True)

# 8. Send dask.DataFrame to database
ddf.to_sql("gene_counts", connection_string, if_exists="fail", index=False, chunksize=10000)
