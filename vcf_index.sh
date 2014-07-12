#!/bin/bash
bcftools view -O z -o $1.gz $1
bcftools index -m0 $1.gz
