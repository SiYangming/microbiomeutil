#!/bin/sh

cmd="../amosCmp16Spipeline.py `pwd`/reads.fasta `pwd`/reads.quals `pwd`/reads.pairs `pwd`/HMP_MOCK_16S_references.fasta"

echo $cmd

$cmd


