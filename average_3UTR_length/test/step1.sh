#!/bin/bash
#3UTR_control.bed is extract from Aseq2 pipline 
#chrom PAC_start PAC_stop strand PAS Total_TPM TPM_rep1 TPM_rep2
#extract three_prime_utr from gtf_file 
#1       10000172        10000455        +       AT1G28440
#1       1000232 1000437 +       AT1G03910
#1       1000232 1000629 +       AT1G03910
#1       10004524        10004581        +       AT1G28450
#1       10011868        10012258        +       AT1G28470
#1       10013434        10013633        -       AT1G28480

intersectBed -a 3UTR_control.bed -b .........../three_prime_utr.gtf.bed -wo > control_3UTR.bed

awk '$4==$12' control_3UTR.bed > control_3UTR_1.bed

awk 'BEGIN{OFS="\t"}{if($4=="-"){print $1,$2,$3,$4,$5,$7,$8,$11,$12,$13}}{if($4=="+"){print $1,$2,$3,$4,$5,$7,$8,$10,$12,$13}}' control_3UTR_1.bed > control_3UTR_2.bed

sort control_3UTR_2.bed | uniq > uniq_control_3UTR_2.bed

awk '{a[$10"\t"$9"\t"$8]=a[$10"\t"$9"\t"$8]$5","$6","$7"\t"}END{for(i in a)print i"\t"a[i]}' uniq_control_3UTR_2.bed > control.bed



