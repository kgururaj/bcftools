#!/bin/bash
bcftools view -O z -o $1.gz $1
bcftools index $1.gz
