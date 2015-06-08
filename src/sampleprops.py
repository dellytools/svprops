#! /usr/bin/env python

from __future__ import print_function
import vcf
import argparse
import csv
import collections

# Parse command line
parser = argparse.ArgumentParser(description='SV Heterozygosity.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-t', '--type', metavar='DEL', required=True, dest='svType', help='SV type [DEL, DUP, INV, INS, TRA] (required)')
parser.add_argument('-s', '--samples', metavar='sample.info', required=True, dest='sampleFile', help='input gender sample info (required)')
args = parser.parse_args()

# Parse gender
resolution = 4
allSuperPops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
gender = dict()
superpop = dict()
pop = dict()
if args.sampleFile:
    with open(args.sampleFile) as f:
        f_reader = csv.reader(f, delimiter="\t")
        for fields in f_reader:
            gender[fields[0]] = fields[3]
            superpop[fields[0]] = fields[1]
            pop[fields[0]] = fields[2]

# Parse vcf
missingCount = collections.Counter()
hetX = collections.Counter()
homX = collections.Counter()
refCount = collections.Counter()
hetCount = collections.Counter()
homCount = collections.Counter()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    for record in vcf_reader:
        if args.svType != record.INFO['SVTYPE']:
            continue
        for call in record.samples:
            if call.called:
                if (record.CHROM == 'X') or (record.CHROM == 'chrX'):
                    if call.gt_type == 1:
                        hetX[call.sample] += 1
                    else:
                        homX[call.sample] += 1
                if call.gt_type == 0:
                    refCount[call.sample] += 1
                elif call['GT'] == "0/1":
                    hetCount[call.sample] += 1
                elif call['GT'] == "1/1":
                    homCount[call.sample] += 1
            else:
                missingCount[call.sample] += 1

    # Summary table
    print("sample", "gender", "population", "superpopulation", "svtype", "missing", "homref", "het", "homalt", "hetX", "homX", sep="\t")
    for s in vcf_reader.samples:
        print(s, gender[s], pop[s], superpop[s], args.svType, missingCount[s], refCount[s], hetCount[s], homCount[s], hetX[s], homX[s], sep="\t")
