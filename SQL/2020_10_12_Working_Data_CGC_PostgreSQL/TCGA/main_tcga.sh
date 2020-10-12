#!/bin/bash

set -e

echo "Installing PostgreSQL server"
sudo apt-get -y update
sudo apt-get -y install postgresql postgresql-contrib

echo "Starting PostgreSQL server"
sudo service postgresql start

echo "Creating SQL database named 'tcga'"
sudo -u postgres psql -q -c "CREATE DATABASE tcga ENCODING 'UTF-8' TEMPLATE template0"

echo "Creating tables and copying data to tcga database"
sudo -u postgres psql -d tcga -a -f "src/tcga_create_tables.sql"

echo "Populate tables"
psql -U postgres -d tcga -c "\COPY allcases FROM 'data/cases.tsv' DELIMITER E'\t' CSV HEADER"
psql -U postgres -d tcga -c "\COPY demographic FROM 'data/demographic.tsv' DELIMITER E'\t' CSV HEADER"
psql -U postgres -d tcga -c "\COPY follow_up FROM 'data/follow_up.tsv' DELIMITER E'\t' CSV HEADER"
psql -U postgres -d tcga -c "\COPY allfiles FROM 'data/files.tsv' DELIMITER E'\t' CSV HEADER"

echo "Installing virtualenv"
pip install virtualenv

echo "Creating virtual environment"
virtualenv venv

echo "Activating virtual environment"
source venv/bin/activate

echo "Installing Python modules"
pip install "dask[complete]" psycopg2-binary sqlalchemy

echo "Creating Python user and role so Python can access the PostgreSQL server"
sudo -u postgres psql -q -c "CREATE USER <USER_NAME> with encrypted password <PASSWORD>"
sudo -u postgres psql -q -c "GRANT ALL PRIVILEGES ON DATABASE tcga TO <USER_NAME>"

echo "Use dask to join copy number files into unified dataframe; export to postgres tcga database; uncomment option inside script to export as one huge csv"
python src/tcga_processing.py

echo "Finishing touches"
psql -U postgres -d tcga -c "CREATE TABLE gene_counts_cases AS SELECT DISTINCT case_id, gene_id, gene_count FROM gene_counts LEFT JOIN allfiles ON gene_counts.filename = allfiles.file_uuid WHERE gene_id LIKE '%ENSG%'"
