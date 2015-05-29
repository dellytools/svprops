#! /usr/bin/env python

from __future__ import print_function
import vcf
import argparse
import csv
import numpy

# Parse command line
parser = argparse.ArgumentParser(description='VCF SV site properties.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-s', '--samples', metavar='sample.info', required=True, dest='sampleFile', help='input gender sample info (required)')
args = parser.parse_args()

resolution = 8

# Parse gender
gender = dict()
if args.sampleFile:
    with open(args.sampleFile) as f:
        f_reader = csv.reader(f, delimiter="\t")
        for fields in f_reader:
            gender[fields[0]] = fields[3]

# Parse vcf
print("chr", "start", "end", "id", "svType", "vac", "af", "size", "sv2size", "ci", "hetperatio", "hetgq", "rarecarrier", "missingrate", sep="\t")
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    for record in vcf_reader:
        svLen = record.INFO['END'] - record.POS + 1
        mci = max([abs(ci) for ci in record.INFO['CIPOS'] + record.INFO['CIEND']])
        try:
            sv2info = record.INFO['SV2INFO']
            sv2size = int(sv2info[2]) - int(sv2info[1])
        except KeyError:
            sv2size = 0
        refAltAC = [0] * (len(record.ALT)+1)
        hetPERatio = list()
        gq = list()
        cs = list()
        uncalled = 0
        for call in record.samples:
            if call.called:
                try:
                    if (call.gt_type == 1) and (call['DV'] > 0) and (call['DR'] > 0):
                        hetPERatio.append(float(call['DV']) / float(call['DR']))
                        gq.append(call['GQ'])
                except AttributeError:
                    pass
                if call.gt_type != 0:
                    cs.append(call.sample)
                if ((record.CHROM == 'X') or (record.CHROM == 'chrX')) and (call.sample in gender.keys()) and (gender[call.sample] == 'male') and (len(call['GT'].split('/')) == 1):
                    gt = int(call['GT'])
                    refAltAC[gt] += 1
                else:
                    for i in [int(gVal) for gVal in call['GT'].split('/')]:
                        refAltAC[i] += 1
            else:
                uncalled += 1
        missingRate = float(uncalled) / float(len(record.samples))
        af = round(float(refAltAC[1])/float(sum(refAltAC)), resolution)
        if len(cs) <= 3:
            cs = ','.join(cs[0:min(3, len(cs))])
        else:
            cs = None
        if len(hetPERatio):
            hper = numpy.median(hetPERatio)
            mgq = numpy.median(gq)
        else:
            hper = None
            mgq = None
        svt = record.INFO['SVTYPE']
        if str(record.ALT[0]) == "<DDEL>":
            if sv2size == 0:
                svt = "DRDEL"
            else:
                svt = "DIDEL"
        elif str(record.ALT[0]) == "<INVDUP>":
            svt = "INVDUP"
        elif str(record.ALT[0]) == "<INVDEL>":
            svt = "INVDEL"
        elif str(record.ALT[0]) == "<PDUP>":
            svt = "PDUP"
        print(record.CHROM, record.POS, record.INFO['END'], record.ID, svt, refAltAC[1], af, svLen, sv2size, mci, hper, mgq, cs, missingRate, sep="\t")
