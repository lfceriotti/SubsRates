#!/bin/bash

# ---------------------------------------------------------------
# Script: subsTree_r2t.sh
# Description:
#   Parses PAML (codeml) output to extract:
#     - Terminal branch omega (w), dN, and dS values
#     - Newick-formatted dN and dS trees
#   Then calls an R script to calculate root-to-tip distances
#   and generate final summary table and annotated tree plots.
#
# Usage:
#   bash subsTree_r2t.sh codeml.ctl output_file.txt
#
# Arguments:
#   $1 - codeml.ctl config file
#   $2 - PAML output file
# ---------------------------------------------------------------

code="$1"       # codeml control file (optional, unused)
file="$2"       # output file from codeml run

### STEP 0: Run codeml ###
arch -arch i386 /usr/local/paml4.8/bin/codeml $code

### STEP 1: Extract omega values per taxon (terminal branches) ###

# Extract tree with omega (w) values from PAML output
awk '/w ratios as labels for TreeView:/{print; getline; print}' "$file" \
  | perl -pe 's/w\ ratios.*\n//g' \
  | perl -pe 's/\ //g' \
  | perl -pe 's/#/:/g' \
  > "${file%%.*}_w.tree"

# Extract taxa and corresponding omega values to a tab-delimited table
echo -e "Taxa\tw" > "${file%%.*}_w.txt"
perl -pe 's/\(//g' "${file%%.*}_w.tree" \
  | perl -pe 's/\)/ /g' \
  | perl -pe 's/\,/\n/g' \
  | perl -pe 's/^\ //g' \
  | perl -pe 's/\ \ .*//g' \
  | perl -pe 's/\ :.*//g' \
  | perl -pe 's/;//g' \
  | perl -pe 's/:/\t/g' \
  | sort -k1 >> "${file%%.*}_w.txt"

### STEP 2: Extract dS and dN trees from codeml output ###

# Get dS tree, remove unnecessary whitespace and internal branch to outgroup (e.g. Liriodendron)
awk '/dS tree/{print; getline; print}' "$file" \
  | perl -pe 's/dS.tree.*\n//g' \
  | perl -pe 's/\ //g' \
  | perl -pe 's/,Liriodendron.*;/\);/g' \
  > "${file%%.*}_dS.tree"

# Get dN tree, processed similarly
awk '/dN tree/{print; getline; print}' "$file" \
  | perl -pe 's/dN.tree.*\n//g' \
  | perl -pe 's/\ //g' \
  | perl -pe 's/,Liriodendron.*;/\);/g' \
  > "${file%%.*}_dN.tree"

### STEP 3: Extract dN and dS terminal branch lengths per taxon ###

# Initialize header lines
echo -e "Taxa\tdN_terminal_branch" > "${file%%.*}_dN.txt"
echo -e "Taxa\tdS_terminal_branch" > "${file%%.*}_dS.txt"

# Parse the trees into taxon/branch_length tables
perl -pe 's/\(//g' "${file%%.*}_dN.tree" \
  | perl -pe 's/\)/ /g' \
  | perl -pe 's/\,/\n/g' \
  | perl -pe 's/^\ //g' \
  | perl -pe 's/\ \ .*//g' \
  | perl -pe 's/\ :.*//g' \
  | perl -pe 's/;//g' \
  | perl -pe 's/:/\t/g' \
  | sort -k1 >> "${file%%.*}_dN.txt"

perl -pe 's/\(//g' "${file%%.*}_dS.tree" \
  | perl -pe 's/\)/ /g' \
  | perl -pe 's/\,/\n/g' \
  | perl -pe 's/^\ //g' \
  | perl -pe 's/\ \ .*//g' \
  | perl -pe 's/\ :.*//g' \
  | perl -pe 's/;//g' \
  | perl -pe 's/:/\t/g' \
  | sort -k1 >> "${file%%.*}_dS.txt"

### STEP 4: Run R script to compute root-to-tip distances ###

# Call the R script with the extracted trees and tables
# NOTE: Update path to script below as needed
/Users/fede/Documents/Laboratorio/scripts/root-to-tip.R \
  "${file%%.*}_dS.tree" \
  "${file%%.*}_dN.tree" \
  "${file%%.*}_w.txt" \
  "${file%%.*}_dN.txt" \
  "${file%%.*}_dS.txt" \
  "${file%%.*}_r2t.tmp"

### STEP 5: Format and sort final output ###

# Clean up quotes, format columns, and sort output
perl -pe 's/"//g' "${file%%.*}_r2t.tmp" \
  | perl -pe 's/Taxa/1\tTaxa/' \
  | cut -f 2,3,4,5,6,7 \
  > r2t.tmp

# Keep header, then sort remaining lines alphabetically by taxon
(head -n 1 r2t.tmp && tail -n +2 r2t.tmp | sort) > "${file%%.*}_r2t.txt"

### Optional cleanup ###
# rm *.tmp *.tree *_w.txt *_dN.txt *_dS.txt

