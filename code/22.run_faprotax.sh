#!bin/bash

code/collapse_table.py \
-i data/intermediate/faprotax.tsv \
-o output/functional_table.tsv \
-g  data/raw/FAPROTAX.txt \
-c "#" \
-d "taxonomy" \
--column_names_are_in first_data_line \
-r data/intermediate/report.txt \
-n columns_after_collapsing \
-v
