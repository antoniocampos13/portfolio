import re
from typing import List


def query_biomart(
    filter_names: List[str],
    values: List[str],
    attribute_names: List[str],
    dataset_name: str = "hsapiens_gene_ensembl",
    formatter: str = "TSV",
    header: bool = False,
    uniqueRows: bool = False
) -> str:
    """Function to build XML query string for data-mining Ensembl's BioMart via RESTful access (datasetConfigVersion = 0.6).

    Args:
        filter_names (List[str]): A list of strings that define restrictions on the query.
        values (List[str]): A list of strings containing the values that will be used for the filters. Must have same length of filter_names.
        attribute_names (List[str]): A list of strings that define the data characteristics that must be retrieved for the filtered data.
        dataset_name (str, optional): A string indicating which dataset will be queried. Each species has their respective dataset. Defaults to "hsapiens_gene_ensembl" (humans).
        formatter (str, optional): A string to indicate how output must be formatted. Options: "HTML", "CSV", "TSV" and "XLS". Defaults to "TSV".
        header (bool, optional): A boolean indicating if the output must have a header. Defaults to False.
        uniqueRows (bool, optional): A boolean indicating if the output must have unique rows (deduplicate data). Defaults to False.

    Returns:
        str: A string containig a wget command to get the desired BioMart data stored inside a file named acording to the given values (the file extension will depend on the formatter chosen: .html, .csv, .txt or .xls). The commands can be conveniently saved to a bash script for easy execution.
    """

    formatter_dict = {"HTML": ".html", "CSV": ".csv", "TSV": ".txt", "XLS": ".xls"}

    value_fix = [val.replace(":", "").replace(".","") for val in values]

    FiltValList = zip(filter_names, values)

    FilterList = [
        f'<Filter name = "{filt}" value = "{val}"/>' for filt, val in FiltValList
    ]

    AttrList = list(map(lambda x: f'<Attribute name = "{x}" />', attribute_names))

    script = f"""wget -O {"_".join(str(val) for val in value_fix)}{formatter_dict[formatter]}
    'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default"
    formatter = "{formatter}"
    header = "{int(header)}"
    uniqueRows = "{int(uniqueRows)}"
    count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "{dataset_name}" interface = "default" >
            {"".join(str(filt) for filt in FilterList)}
            {"".join(str(attr) for attr in AttrList)}
        </Dataset>
    </Query>'"""

    return re.sub("\n|( )+", " ", script)
