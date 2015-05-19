__author__ = 'kssindy'


#from Bio import SeqIO
import sys




#Sets input file name
#argv[1] is original sequencing file with all barcodes mixed together

#KH122
OutFileName1 = “CAG_SU2.fastq”
#KH123
OutFileName2 = “GTC_SU2nodox.fastq"
#KH126
OutFileName3 = “CTA_SU2PB.fastq"

#KH122
OutFileName4 = “CAG_SU2-Pair.fastq"
#KH123
OutFileName5 = “GTC_SU2nodox-Pair.fastq"
#KH126
OutFileName6 = “CTA_SU2PB-Pair.fastq"


OutFileName7 = "NoMatch"


#Names output file for writing
#
OutFile1 = open(OutFileName1, 'w')
OutFile2 = open(OutFileName2, 'w')
OutFile3 = open(OutFileName3, 'w')

OutFile4 = open(OutFileName4, 'w')
OutFile5 = open(OutFileName5, 'w')
OutFile6 = open(OutFileName6, 'w')

OutFile7 = open(OutFileName7, 'w')



#Creates output file for writing


CAG_find = “CAG”
GTC_find = “GTC”
CTA_find = “CTA”


# for record in SeqIO.parse(KS01_R1, "fastq"):
#     #print record.seq
#     for recordR in SeqIO.parse(KS01_R2, "fastq"):
#         #print recordR.seq

forward_name = ""
reverse_name = ""

sawaCAG = False
sawaGTC = False
sawaCTA = False
#sawaTGAseq = False

KS01_R1 = open(sys.argv[1])
KS01_R2 = open(sys.argv[2])

for forward, reverse in zip(KS01_R1, KS01_R2):
    if (forward.startswith('@MSQ-M01247R')):
        forward_name = forward
        reverse_name = reverse


    elif forward.startswith(CAG_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaACG = True


    elif sawaCAG and (not forward.startswith('+\n')):
        if not forward.startswith("@MSQ-M01247R"):
            forward_fastq1 = forward_name
            forward_fastq2 = forward_seq
            forward_fastq3 = '+\n'
            forward_fastq4 = forward


            forward_fastq = forward_fastq1 + forward_fastq2 + forward_fastq3 + forward_fastq4
            OutFile1.write(forward_fastq)

            reverse_fastq1 = reverse_name
            reverse_fastq2 = reverse_seq
            reverse_fastq3 = '+\n'
            reverse_fastq4 = reverse

            reverse_fastq = reverse_fastq1 + reverse_fastq2 + reverse_fastq3 + reverse_fastq4
            OutFile4.write(reverse_fastq)
            sawaCAG = False

    elif forward.startswith(GTC_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaGTC = True

    elif sawaGTC and (not forward.startswith('+\n')):
        if not forward.startswith("@MSQ-M01247R"):
            forward_fastq1 = forward_name
            forward_fastq2 = forward_seq
            forward_fastq3 = '+\n'
            forward_fastq4 = forward


            forward_fastq = forward_fastq1 + forward_fastq2 + forward_fastq3 + forward_fastq4
            OutFile2.write(forward_fastq)

            reverse_fastq1 = reverse_name
            reverse_fastq2 = reverse_seq
            reverse_fastq3 = '+\n'
            reverse_fastq4 = reverse

            reverse_fastq = reverse_fastq1 + reverse_fastq2 + reverse_fastq3 + reverse_fastq4
            OutFile5.write(reverse_fastq)
            sawaGTC = False



    elif forward.startswith(CTA_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaCTA = True

    elif sawaCTA and (not forward.startswith('+\n')):
        if not forward.startswith("@M00658"):
            forward_fastq1 = forward_name
            forward_fastq2 = forward_seq
            forward_fastq3 = '+\n'
            forward_fastq4 = forward


            forward_fastq = forward_fastq1 + forward_fastq2 + forward_fastq3 + forward_fastq4
            OutFile3.write(forward_fastq)

            reverse_fastq1 = reverse_name
            reverse_fastq2 = reverse_seq
            reverse_fastq3 = '+\n'
            reverse_fastq4 = reverse

            reverse_fastq = reverse_fastq1 + reverse_fastq2 + reverse_fastq3 + reverse_fastq4
            OutFile6.write(reverse_fastq)
            sawaCTA = False



