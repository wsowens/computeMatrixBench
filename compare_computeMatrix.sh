#!/bin/bash

TEST_BW="input/random.bigWig"
TEST_BED="input/random_genes.bed"

OUT_BG="output/bigtools.bg"
MY_MAT="output/test_matrix.tsv"
OUT_MAT="output/deeptools.mat.gz"
mkdir -p output

echo "Testing deeptools computeMatrix reference-point against Will's ~200 line script"
set -x
date

: "Testing Will's script with this version of bigtools"
bigtools --version
time ./compute_matrix_faster.py "$TEST_BW" "$TEST_BED" > "$MY_MAT"

: "Testing deeptools computeMatrix"
deeptools --version
#time computeMatrix reference-point -S "$TEST_BW" -o "$OUT_MAT" -R "$TEST_BED"

OUT_MAT_BS1="output/$(basename "$OUT_MAT" .mat.gz).bs1.mat.gz"
time computeMatrix reference-point -bs 1 -S "$TEST_BW" -o "$OUT_MAT_BS1" -R "$TEST_BED"

set +x
#gunzip -kf "$OUT_MAT"
gunzip -kf "$OUT_MAT_BS1"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "Confirming that matrices are the same..."
# for now, just ignore the parameter info on the first line looks like this:
# '@{"upstream":[500],"downstream":[1500],"body":[0],"bin size":[1],"
zcat "$OUT_MAT_BS1" | tail -n+2 | diff -q  - "$MY_MAT" 

if [ $? -eq 0 ]
then
    echo -e "${GREEN}SUCCESS! :)${NC}"
else
    echo -e "${RED}Fail :( Matrix values differ.${NC}"
fi

