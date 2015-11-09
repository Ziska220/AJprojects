__author__ = 'kssindy'


import sys
import datetime
import time
import re
import difflib

BarcodeDict={'GTTATGAAGG':1, 'ATCACTTAAG':2, 'ATAGCTCAGA':3, 'CGCCCTCGCA':4, 'ACCAAAAAAC':5, 'GACGGGGGTG':6, 'TCGAAACATA':7, 'GAGTCGTCTG':8, 'ATCTACCTGA':9, 'TCTGCTAGTT':10}

Barcodes = BarcodeDict.keys()

Barcodes_File = ""

Barcode_list = []

AcceptRatio = 0.9

mspi = '[ACGTN]{11}(CGG)[ACGTN]{62}'
#MspI is CCGG
csp6i = '[ACGTN]{11}(CGGTAC)[ACGTN]{59}'
#Csp6I is GTAC
taqi = '[ACGTN]{11}(CGGTACTCGA)[ACGTN]{55}'
#TaqI = TCGA


ts1 = time.time()
st = datetime.datetime.fromtimestamp(ts1).strftime('%Y-%m-%d %H:%M:%S')

OutFileName1 = "taqi_forward.fastq"
OutFileName2 = "taqi_reverse.fastq"
OutFileName3 = "csp6i_forward.fastq"
OutFileName4 = "csp6i_reverse.fastq"
OutFileName5 = "mspi_forward.fastq"
OutFileName6 = "mspi_reverse.fastq"

OutFile1 = open(OutFileName1, 'w')
OutFile2 = open(OutFileName2, 'w')
OutFile3 = open(OutFileName3, 'w')
OutFile4 = open(OutFileName4, 'w')
OutFile5 = open(OutFileName5, 'w')
OutFile6 = open(OutFileName6, 'w')

sawtaqi = False
sawcsp6i = False
sawmspi = False

sawqual = False
sawqualcsp61 = False
sawqualmspi = False


KS01_R1 = open(sys.argv[1])
KS01_R2 = open(sys.argv[2])

for reverse, forward in zip(KS01_R2, KS01_R1):

    if reverse.startswith('@M00658'):
        reverse_name = reverse
        forward_name = forward


    for keys in Barcodes:

        ten_string = reverse[0:10]

        BarcodeRatio = difflib.SequenceMatcher(None, keys, ten_string).ratio()

        if BarcodeRatio >= AcceptRatio and re.search(taqi,reverse):
            barcode_seq_match_reverse = reverse
            barcode_seq_match_forward = forward

            sawtaqi = True

        elif sawtaqi and reverse.startswith('+\n'):
            reverse_plus = reverse
            forward_plus = forward

            sawtaqi = False
            sawqual = True

        elif sawqual and not reverse.startswith('@M00658'):
            if not reverse.startswith('+\n'):
                seq_qual_forward = forward
                seq_qual_reverse = reverse

                sawqual = False


                taqi_fastq_forward = forward_name + barcode_seq_match_forward + forward_plus + seq_qual_forward
                taqi_fastq_reverse = reverse_name + barcode_seq_match_reverse + reverse_plus + seq_qual_reverse
                OutFile1.write(taqi_fastq_forward)
                OutFile2.write(taqi_fastq_reverse)



        elif BarcodeRatio >= AcceptRatio and re.search(csp6i,reverse):
            barcode_seq_match_forward = forward
            barcode_seq_match_reverse = reverse

            sawcsp6i = True

        elif sawcsp6i and reverse.startswith('+\n'):
            forward_plus = forward
            reverse_plus = reverse

            sawcsp6i = False
            sawqualcsp61 = True

        elif sawqualcsp61 and not reverse.startswith('@M00658'):
            if not reverse.startswith('+\n'):
                seq_qual_forward = forward
                seq_qual_reverse = reverse
                sawqualcsp61 = False


                csp6i_fastq_forward = forward_name + barcode_seq_match_forward + forward_plus + seq_qual_forward
                csp6i_fastq_reverse = reverse_name + barcode_seq_match_reverse + reverse_plus + seq_qual_reverse
                OutFile3.write(csp6i_fastq_forward)
                OutFile4.write(csp6i_fastq_reverse)


        elif BarcodeRatio >= AcceptRatio and re.search(mspi,reverse):
            barcode_seq_match_forward = forward
            barcode_seq_match_reverse = reverse

            sawmspi = True

        elif sawmspi and reverse.startswith('+\n'):
            forward_plus = forward
            reverse_plus = reverse

            sawmspi = False
            sawqualmspi = True

        elif sawqualmspi and not reverse.startswith('@M00658'):
            if not reverse.startswith('+\n'):
                seq_qual_forward = forward
                seq_qual_reverse = reverse
                sawqualmspi = False


                mspi_fastq_forward = forward_name + barcode_seq_match_forward + forward_plus + seq_qual_forward
                mspi_fastq_reverse = reverse_name + barcode_seq_match_reverse + reverse_plus + seq_qual_reverse
                OutFile5.write(mspi_fastq_forward)
                OutFile6.write(mspi_fastq_reverse)
