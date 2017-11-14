# 3kRG-SNPUniverse

Author: Roven Rommel Fuentes

This utility merges variant calls from multiple VCF files. It applies the threshold QUAL>=30
and assumes that GATK UnifiedGenotypr was ran with -emit-all-sites parameter to detect both SNPs and Indels.
In retrieving variants, it considers issues caused by running GATK UG with SNPs and Indels together. 
For parallelization, this code runs with a set of samples assigned to specific thread.
There is also a feature to perform reduction when detecting variant calls from multiple references; it follows
an order of the references and also recovers new SNP and Indel positions.
The output format is a tab delimited file ready for loading to SNP-Seek database. 
