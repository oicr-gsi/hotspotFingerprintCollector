#!/bin/bash
set -euo pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all vcf and fin files, return their md5sums to std out, list all file types
find . -name *.vcf.gz -xtype f -exec sh -c "zcat {} | grep -v ^# | md5sum" \;
find . -name *.fin -xtype f -exec sh -c "cat {} | md5sum" \;
ls | sed 's/.*\.//' | sort | uniq -c
