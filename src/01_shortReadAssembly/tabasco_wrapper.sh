#!/usr/bin/env bash

#--------------------------------------------------------------------
# Script Name:   tabasco_wrapper.sh
# Description:   Simple wrapper for running TABASCO
# Author:        Brandon Monier, Charlie Hale
# Created:       2021-06-15 at 11:52:44
# Last Modified: 2024-01-08 at 14:27:01
#--------------------------------------------------------------------

FASTA_FILE=$1
TRANSCRIPT_FILE=$2
THRESHOLD=$3
OUTPUT_DIR=$4

perl src/01_shortReadAssembly/tabasco_pipeline.pl \
    $FASTA_FILE \
    $TRANSCRIPT_FILE \
    $THRESHOLD \
    $OUTPUT_DIR
