__author__ = 'kssindy'

from Bio import SeqIO
import sys
import datetime
import time
import difflib
import os
import os.path
import shutil
import subprocess

import re

ts1 = time.time()
st = datetime.datetime.fromtimestamp(ts1).strftime('%Y-%m-%d %H:%M:%S')

mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-1'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-2'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-3'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-4'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-5'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-6'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-7'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-8'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-9'
if not os.path.exists(mypath): os.makedirs(mypath)
mypath = r'/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-10'
if not os.path.exists(mypath): os.makedirs(mypath)

OutFileName1 = "Barcode-1-R.fastq"
OutFileName2 = "Barcode-2-R.fastq"
OutFileName3 = "Barcode-3-R.fastq"
OutFileName4 = "Barcode-4-R.fastq"
OutFileName5 = "Barcode-5-R.fastq"
OutFileName6 = "Barcode-6-R.fastq"
OutFileName7 = "Barcode-7-R.fastq"
OutFileName8 = "Barcode-8-R.fastq"
OutFileName9 = "Barcode-9-R.fastq"
OutFileName10 = "Barcode-10-R.fastq"
OutFileName11 = "Barcode-1-F.fastq"
OutFileName12 = "Barcode-2-F.fastq"
OutFileName13 = "Barcode-3-F.fastq"
OutFileName14 = "Barcode-4-F.fastq"
OutFileName15 = "Barcode-5-F.fastq"
OutFileName16 = "Barcode-6-F.fastq"
OutFileName17 = "Barcode-7-F.fastq"
OutFileName18 = "Barcode-8-F.fastq"
OutFileName19 = "Barcode-9-F.fastq"
OutFileName20 = "Barcode-10-F.fastq"

OutFile1 = open(OutFileName1, 'w')
OutFile2 = open(OutFileName2, 'w')
OutFile3 = open(OutFileName3, 'w')
OutFile4 = open(OutFileName4, 'w')
OutFile5 = open(OutFileName5, 'w')
OutFile6 = open(OutFileName6, 'w')
OutFile7 = open(OutFileName7, 'w')
OutFile8 = open(OutFileName8, 'w')
OutFile9 = open(OutFileName9, 'w')
OutFile10 = open(OutFileName10, 'w')
OutFile11 = open(OutFileName11, 'w')
OutFile12 = open(OutFileName12, 'w')
OutFile13 = open(OutFileName13, 'w')
OutFile14 = open(OutFileName14, 'w')
OutFile15 = open(OutFileName15, 'w')
OutFile16 = open(OutFileName16, 'w')
OutFile17 = open(OutFileName17, 'w')
OutFile18 = open(OutFileName18, 'w')
OutFile19 = open(OutFileName19, 'w')
OutFile20 = open(OutFileName20, 'w')


KS01_R1 = open(sys.argv[1])
KS01_R2 = open(sys.argv[2])

BarcodeDict={'GTTATGAAGG':1, 'ATCACTTAAG':2, 'ATAGCTCAGA':3, 'CGCCCTCGCA':4, 'ACCAAAAAAC':5, 'GACGGGGGTG':6, 'TCGAAACATA':7, 'GAGTCGTCTG':8, 'ATCTACCTGA':9, 'TCTGCTAGTT':10}

Barcodes = BarcodeDict.keys()
OneBarcodeSeqReverse = ""
OneBarcodeSeqForward = ""

TwoBarcodeSeqReverse = ""
TwoBarcodeSeqForward = ""

ThreeBarcodeSeqReverse = ""
ThreeBarcodeSeqForward = ""

FourBarcodeSeqReverse = ""
FourBarcodeSeqForward = ""

FiveBarcodeSeqReverse = ""
FiveBarcodeSeqForward = ""

SixBarcodeSeqReverse = ""
SixBarcodeSeqForward = ""

SevenBarcodeSeqReverse = ""
SevenBarcodeSeqForward = ""

EightBarcodeSeqReverse = ""
EightBarcodeSeqForward = ""

NineBarcodeSeqReverse = ""
NineBarcodeSeqForward = ""

TenBarcodeSeqReverse = ""
TenBarcodeSeqForward = ""


#print Barcodes

Barcodes_File = ""

AcceptRatio = 0.9

Barcode_list = []

sawBarcodeOne = False
sawPlusOne = False

sawBarcodeTwo = False
sawPlusTwo = False

sawBarcodeThree = False
sawPlusThree = False

sawBarcodeFour = False
sawPlusFour = False

sawBarcodeFive = False
sawPlusFive = False

sawBarcodeSix = False
sawPlusSix = False

sawBarcodeSeven = False
sawPlusSeven = False

sawBarcodeEight = False
sawPlusEight = False

sawBarcodeNine = False
sawPlusNine = False

sawBarcodeTen = False
sawPlusTen = False




for reverse, forward in zip(KS01_R2, KS01_R1):
    ten_string = reverse[0:10]


    if reverse.startswith('@M00658'):
        reverse_name = reverse
        forward_name = forward

    for keys in Barcodes:

        BarcodeRatio = difflib.SequenceMatcher(None, keys, ten_string).ratio()
        if BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 1:
            OneBarcodeSeqReverse = reverse
            OneBarcodeSeqForward = forward
            sawBarcodeOne = True

        elif sawBarcodeOne and reverse.startswith('+\n'):
            forward_plus = forward
            reverse_plus = reverse
            sawBarcodeOne = False
            sawPlusOne = True

        elif sawPlusOne and not reverse.startswith('@M00658'):
            if not reverse.startswith('+\n'):
                seq_qual_forward = forward
                seq_qual_reverse = reverse
                sawPlusOne = False

                BarcodeOneFastqReverse = reverse_name + OneBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                OutFile1.write(BarcodeOneFastqReverse)
                BarcodeOneFastqForward = forward_name + OneBarcodeSeqForward + forward_plus + seq_qual_forward
                OutFile11.write(BarcodeOneFastqForward)


        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 2:
            TwoBarcodeSeqReverse = reverse
            TwoBarcodeSeqForward = forward
            sawBarcodeTwo = True

        elif sawBarcodeTwo and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeTwo = False
                sawPlusTwo = True

        elif sawPlusTwo and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusTwo = False

                    BarcodeTwoFastqReverse = reverse_name + TwoBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeTwoFastqForward = forward_name + TwoBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile2.write(BarcodeTwoFastqReverse)
                    OutFile12.write(BarcodeTwoFastqForward)
                    print BarcodeTwoFastqForward


        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 3:
            ThreeBarcodeSeqReverse = reverse
            ThreeBarcodeSeqForward = forward
            sawBarcodeThree = True

        elif sawBarcodeThree and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeThree = False
                sawPlusThree = True

        elif sawPlusThree and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusThree = False
                    BarcodeThreeFastqReverse = reverse_name + ThreeBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeThreeFastqForward = forward_name + ThreeBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile3.write(BarcodeThreeFastqReverse)
                    OutFile13.write(BarcodeThreeFastqForward)

        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 4:
            FourBarcodeSeqReverse = reverse
            FourBarcodeSeqForward = forward
            sawBarcodeFour = True

        elif sawBarcodeFour and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeFour = False
                sawPlusFour = True

        elif sawPlusFour and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusFour = False

                    BarcodeFourFastqReverse = reverse_name + FourBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeFourFastqForward = forward_name + FourBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile4.write(BarcodeFourFastqReverse)
                    OutFile14.write(BarcodeFourFastqForward)

        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 5:
            FiveBarcodeSeqReverse = reverse
            FiveBarcodeSeqForward = forward
            sawBarcodeFive = True

        elif sawBarcodeFive and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeFive = False
                sawPlusFive = True

        elif sawPlusFive and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusFive = False
                    BarcodeFiveFastqReverse = reverse_name + FiveBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeFiveFastqForward = forward_name + FiveBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile5.write(BarcodeFiveFastqReverse)
                    OutFile15.write(BarcodeFiveFastqForward)

        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 6:
            SixBarcodeSeqReverse = reverse
            SixBarcodeSeqForward = forward
            sawBarcodeSix = True

        elif sawBarcodeSix and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeSix = False
                sawPlusSix = True

        elif sawPlusSix and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusSix = False

                    BarcodeSixFastqReverse = reverse_name + SixBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeSixFastqForward = forward_name + SixBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile6.write(BarcodeSixFastqReverse)
                    OutFile16.write(BarcodeSixFastqForward)


        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 7:
            SevenBarcodeSeqReverse = reverse
            SevenBarcodeSeqForward = forward
            sawBarcodeSeven = True

        elif sawBarcodeSeven and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeSeven = False
                sawPlusSeven = True

        elif sawPlusSeven and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusSeven = False

                    BarcodeSevenFastqReverse = reverse_name + SevenBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeSevenFastqForward = forward_name + SevenBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile7.write(BarcodeSevenFastqReverse)
                    OutFile17.write(BarcodeSevenFastqForward)


        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 8:
            EightBarcodeSeqReverse = reverse
            EightBarcodeSeqForward = forward
            sawBarcodeEight = True

        elif sawBarcodeEight and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeEight = False
                sawPlusEight = True

        elif sawPlusEight and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusEight = False

                    BarcodeEightFastqReverse = reverse_name + EightBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeEightFastqForward = forward_name + EightBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile8.write(BarcodeEightFastqReverse)
                    OutFile18.write(BarcodeEightFastqForward)


        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 9:
            NineBarcodeSeqReverse = reverse
            NineBarcodeSeqForward = forward
            sawBarcodeNine = True

        elif sawBarcodeNine and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeNine = False
                sawPlusNine = True

        elif sawPlusNine and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusNine = False
                    BarcodeNineFastqReverse = reverse_name + NineBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeNineFastqForward = forward_name + NineBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile9.write(BarcodeNineFastqReverse)
                    OutFile19.write(BarcodeNineFastqForward)


        elif BarcodeRatio >= AcceptRatio and BarcodeDict.get(keys) == 10:
            TenBarcodeSeqReverse = reverse
            TenBarcodeSeqForward = forward
            sawBarcodeTen = True

        elif sawBarcodeTen and reverse.startswith('+\n'):
                forward_plus = forward
                reverse_plus = reverse
                sawBarcodeTen = False
                sawPlusTen = True

        elif sawPlusTen and not reverse.startswith('@M00658'):
                if not reverse.startswith('+\n'):
                    seq_qual_forward = forward
                    seq_qual_reverse = reverse
                    sawPlusTen = False

                    BarcodeTenFastqReverse = reverse_name + TenBarcodeSeqReverse + reverse_plus + seq_qual_reverse
                    BarcodeTenFastqForward = forward_name + TenBarcodeSeqForward + forward_plus + seq_qual_forward
                    OutFile10.write(BarcodeTenFastqReverse)
                    OutFile20.write(BarcodeTenFastqForward)


OutFile1.close()
OutFile11.close()


Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-1"
Source1 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-1-R.fastq"
Source2 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-1-F.fastq"
shutil.move(Source1, Destination)
shutil.move(Source2, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-1')
shutil.move('Barcode-1-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-1-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-2"
Source3 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-2-R.fastq"
Source4 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-2-F.fastq"
shutil.move(Source3, Destination)
shutil.move(Source4, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-2')
shutil.move('Barcode-2-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-2-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-3"
Source5 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-3-R.fastq"
Source6 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-3-F.fastq"
shutil.move(Source5, Destination)
shutil.move(Source6, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-3')
shutil.move('Barcode-3-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-3-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-4"
Source7 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-4-R.fastq"
Source8 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-4-F.fastq"
shutil.move(Source7, Destination)
shutil.move(Source8, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-4')
shutil.move('Barcode-4-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-4-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-5"
Source9 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-5-R.fastq"
Source10 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-5-F.fastq"
shutil.move(Source9, Destination)
shutil.move(Source10, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-5')
shutil.move('Barcode-5-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-5-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-6"
Source11 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-6-R.fastq"
Source12 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-6-F.fastq"
shutil.move(Source11, Destination)
shutil.move(Source12, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-6')
shutil.move('Barcode-6-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-6-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-7"
Source13 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-7-R.fastq"
Source14 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-7-F.fastq"
shutil.move(Source13, Destination)
shutil.move(Source14, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-7')
shutil.move('Barcode-7-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-7-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-8"
Source15 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-8-R.fastq"
Source16 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-8-F.fastq"
shutil.move(Source15, Destination)
shutil.move(Source16, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-8')
shutil.move('Barcode-8-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-8-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-9"
Source17 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-9-R.fastq"
Source18 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-9-F.fastq"
shutil.move(Source17, Destination)
shutil.move(Source18, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-9')
shutil.move('Barcode-9-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-9-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign')

Destination = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-10"
Source19 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-10-R.fastq"
Source20 = "/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-10-F.fastq"
shutil.move(Source19, Destination)
shutil.move(Source20, Destination)

os.chdir('/vol1/home/sindyk/2013-12-23CC_5201_3202-2014-03-07_MiSeq/2014-03-21_ACG-5/csp6i_ACG-5/10BPandAlign/Barcode-10')
shutil.move('Barcode-10-R.fastq', 'Barcode-R.fastq')
shutil.move('Barcode-10-F.fastq', 'Barcode-F.fastq')
subprocess.call("../runalign.sh", shell=True)

