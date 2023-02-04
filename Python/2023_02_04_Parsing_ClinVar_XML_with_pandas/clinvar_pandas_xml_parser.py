# %% Install pandas=1.5.0 and lxml modules
import pandas as pd

# %% Path to CLinVar's full release huge XML. Updated weekly. 
## Only the release on the first Thursday of the month is archived. 
## Link for the most recent monthly archive: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml
## Link for the most recent weekly version: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarFullRelease_00-latest_weekly.xml.gz
xml_file_path = "ClinVarFullRelease_00-latest.xml"

# %% The key is the parent XML node. The value is a list contanining all child or grandchild nodes, tags, or attributes at any node level.
iterparse_dict = {
    "ClinVarAssertion": [
        "ID",
        "SubmissionName",
        "localKey",
        "submittedAssembly",
        "submitter",
        "submitterDate",
        "Acc",
        "RecordStatus",
        "OrgID",
        "DateCreated",
        "DateUpdated",
        "Version",
        "DateLastEvaluated",
        "Description",
        "ReviewStatus",
        "Comment"
    ]
}

# %%
df = pd.read_xml(xml_file_path, iterparse=iterparse_dict)

# %% Save data frame as a pickled object
df.to_pickle("pandas_parsed.pkl")

# %% Restore pickle. REMEMBER: Loading pickled data received from untrusted sources can be unsafe
df = pd.read_pickle("pandas_parsed.pkl")