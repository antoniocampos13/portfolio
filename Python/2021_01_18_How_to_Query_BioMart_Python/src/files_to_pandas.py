import ntpath
import pandas as pd

def files_to_pandas(filename):
    """ Converts each file into a pandas.DataFrame """

    tmp = pd.read_csv(filename, delimiter="\t", names=["gene_name"])

    tmp["filename"] = ntpath.basename(filename).replace("GO", "GO:").replace(".txt","")
    return tmp
    