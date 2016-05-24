#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l mem=20gb
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

echo "testing"
#grep "Chr1\t" -P /rhome/cjinfeng/BigData/00.RD/Variations/SV_postprocess/Local_Assembly/bin/insertion.gff > insertion.gff
#perl GenomeValid.pl --step 1234 > log 2> log2 &

echo "HuRef"
#perl GenomeValid.pl --step 1234 -fa /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Validation/Genome_Validation/input/hg18.fa --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Validation/Genome_Validation/input/HuRef.fa --bam HuRef.bam.list --gff insertion.HuRef_RelocaTE2.gff --project HuRef_Genome_RelocaTE2
#perl GenomeValid.pl --step 1234 -fa /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Validation/Genome_Validation/input/hg18.fa --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Validation/Genome_Validation/input/HuRef.fa --bam HuRef.bam.list --gff insertion.HuRef_TEMP.gff --project HuRef_Genome_TEMP

echo "IR64"
#perl GenomeValid_new.pl --step 1234 -fa /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Validation/Genome_Validation/input/IR64.RM.fa --gff ALL.all_nonref_insert.gff --project IR64_Genome_RelocaTE2 
perl GenomeValid_new.pl --step 1234 -fa /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Validation/Genome_Validation/input/IR64.RM.fa --gff IR64_TEMP_TE14_ns.gff --project IR64_Genome_TEMP 


echo "Done"
