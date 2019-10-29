# python vcfCompare.py

import pysam
import argparse
import sys
import math

vcfReal = pysam.VariantFile('../HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz','r')

counts = {"match": 0, "inTrue": 0, "inTest": 0, "samePos_noMatch":0}

startIndex = 1000000
endIndex = 1050000
chrom = "chr20"

with open('cigar_output.vcf', "r") as vcfTest:
    header = vcfTest.readline()
    pos = startIndex - 1
    for var in vcfReal.fetch(chrom, start=startIndex, end=endIndex):
        counts["inTrue"] += 1
        rchr = var.chrom
        rref = var.ref
        ralt = var.alts
        rpos = var.pos
        print(rchr, rref, ralt[0], rpos)
        while pos < rpos:
            line = vcfTest.readline()
            if line:
                line = line.rstrip().split("\t")
                chr = line[0]
                ref = line[3]
                alt = line[4]
                pos = int(line[1])
                counts["inTest"] += 1

            else:
                break

        if rpos == pos:
            for a in ralt:
                if rchr == chr and a == alt and rref == ref:
                    # print("good", chr, alt, ref, pos)
                    counts["match"] += 1
                else:
                    counts["samePos_noMatch"] += 1

    print(counts)
