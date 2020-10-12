# "Creating SQL database named 'tcga'"
psql -U postgres -c "CREATE DATABASE tcga ENCODING 'UTF-8' LC_COLLATE 'English_United States' LC_CTYPE 'English_United States' TEMPLATE template0"

# Create cases, demographic, follow_up and files tables
psql -U postgres -d tcga -a -f "src/tcga_create_tables.sql"

# Populate tables created above
psql -U postgres -d tcga -c "\COPY allcases FROM 'data/cases.tsv' DELIMITER E'\t' CSV HEADER"
psql -U postgres -d tcga -c "\COPY demographic FROM 'data/demographic.tsv' DELIMITER E'\t' CSV HEADER"
psql -U postgres -d tcga -c "\COPY follow_up FROM 'data/follow_up.tsv' DELIMITER E'\t' CSV HEADER"
psql -U postgres -d tcga -c "\COPY allfiles FROM 'data/files.tsv' DELIMITER E'\t' CSV HEADER"

# Install virtualenv if necessary
python -m pip install --user virtualenv

# Create virtual environment with name 'venv'
python -m venv venv

# Activate virtual environment
venv\Scripts\activate

# Install Python modules
pip install "dask[complete]" psycopg2-binary sqlalchemy

# Creating Python user and role so Python can access the PostgreSQL server
psql -U postgres -d tcga -c "CREATE USER <USER_NAME> with encrypted password '<PASSWORD>'"
psql -U postgres -d tcga -c "GRANT ALL PRIVILEGES ON DATABASE tcga TO <USER_NAME>"

# Use dask to join files into unified dataframes; export to postgres tcga database; uncomment option inside script to export as one huge csv
python src\tcga_processing_counts.py

# Finishing touches
psql -U postgres -d tcga -c "CREATE TABLE gene_counts_cases AS SELECT DISTINCT case_id, gene_id, gene_count FROM gene_counts LEFT JOIN allfiles ON gene_counts.filename = allfiles.file_uuid WHERE gene_id LIKE '%ENSG%'"
