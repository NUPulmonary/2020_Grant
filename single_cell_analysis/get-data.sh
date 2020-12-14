#!/bin/bash

if [[ ! -d data ]]; then
    mkdir data
fi

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155249&format=file" -O data/GSE155249_RAW.tar
cd data
tar xf GSE155249_RAW.tar
cd -
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155249&format=file&file=GSE155249_main-metadata%2Ecsv%2Egz" -O data/main-metadata.csv.gz
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155249&format=file&file=GSE155249_supplement-metadata%2Ecsv%2Egz" -O data/supplement-metadata.csv.gz

gunzip data/main-metadata.csv.gz
gunzip data/supplement-metadata.csv.gz

