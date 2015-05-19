__author__ = 'kssindy'


#from Bio import SeqIO
import sys




#Sets input file name
#argv[1] is original sequencing file with all barcodes mixed together

OutFileName1 = "ACG_SU.fastq"
OutFileName2 = "CTA_SUnodox_NoDox.fastq"
OutFileName3 = "GAT_SUPB.fastq"
OutFileName4 = "TGC_MDA.fastq"
OutFileName5 = "TGA_AL.fastq"

OutFileName6 = "ACG_SU-pair.fastq"
OutFileName7 = "CTA_SUnodox-pair.fastq"
OutFileName8 = "GAT_SUPB-pair.fastq"
OutFileName9 = "TGC_MDA-pair.fastq"
OutFileName10 = "TGA_AL-pair.fastq"

OutFileName11 = "NoMatch"


#Names output file for writing
#
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

OutFile11 = open(OutFileName10, 'w')

#Creates output file for writing


ACG_find = "ACG"
CTA_find = "CTA"
GAT_find = "GAT"
TGC_find = "TGC"
TGA_find = "TGA"


# for record in SeqIO.parse(KS01_R1, "fastq"):
#     #print record.seq
#     for recordR in SeqIO.parse(KS01_R2, "fastq"):
#         #print recordR.seq

forward_name = ""
reverse_name = ""

sawaACG = False
sawaCTAseq = False
sawaGATseq = False
sawaTGCseq = False
sawaTGAseq = False

KS01_R1 = open(sys.argv[1])
KS01_R2 = open(sys.argv[2])

for forward, reverse in zip(KS01_R1, KS01_R2):
    if (forward.startswith('@M00658')):
        forward_name = forward
        reverse_name = reverse


    elif forward.startswith(ACG_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaACG = True


    elif sawaACG and (not forward.startswith('+\n')):
        if not forward.startswith("@M00658"):
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
            OutFile6.write(reverse_fastq)
            sawaACG = False

    elif forward.startswith(CTA_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaCTAseq = True

    elif sawaCTAseq and (not forward.startswith('+\n')):
        if not forward.startswith("@M00658"):
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
            OutFile7.write(reverse_fastq)
            sawaCTAseq = False



    elif forward.startswith(GAT_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaGATseq = True

    elif sawaGATseq and (not forward.startswith('+\n')):
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
            OutFile8.write(reverse_fastq)
            sawaGATseq = False

    elif forward.startswith(TGC_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaTGCseq = True

    elif sawaTGCseq and (not forward.startswith('+\n')):
        if not forward.startswith('@M00658'):
            forward_fastq1 = forward_name
            forward_fastq2 = forward_seq
            forward_fastq3 = '+\n'
            forward_fastq4 = forward


            forward_fastq = forward_fastq1 + forward_fastq2 + forward_fastq3 + forward_fastq4
            OutFile4.write(forward_fastq)

            reverse_fastq1 = reverse_name
            reverse_fastq2 = reverse_seq
            reverse_fastq3 = '+\n'
            reverse_fastq4 = reverse

            reverse_fastq = reverse_fastq1 + reverse_fastq2 + reverse_fastq3 + reverse_fastq4
            OutFile9.write(reverse_fastq)
            sawaTGCseq = False

    elif forward.startswith(TGA_find) == True:
        reverse_seq = reverse
        forward_seq = forward
        sawaTGAseq = True

    elif sawaTGAseq and (not forward.startswith('+\n')):
        if not forward.startswith('@M00658'):
            forward_fastq1 = forward_name
            forward_fastq2 = forward_seq
            forward_fastq3 = '+\n'
            forward_fastq4 = forward


            forward_fastq = forward_fastq1 + forward_fastq2 + forward_fastq3 + forward_fastq4
            OutFile5.write(forward_fastq)

            reverse_fastq1 = reverse_name
            reverse_fastq2 = reverse_seq
            reverse_fastq3 = '+\n'
            reverse_fastq4 = reverse

            reverse_fastq = reverse_fastq1 + reverse_fastq2 + reverse_fastq3 + reverse_fastq4
            OutFile10.write(reverse_fastq)
            sawaTGCseq = False

