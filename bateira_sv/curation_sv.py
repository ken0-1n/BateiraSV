import os, sys
import subprocess
import pysam
import edlib
from . import utils

def curation1(input_file, in_bam1, in_bam2, output_file, margin, reference_genome, max_depth, validate_sequence_length, ed_threas = 0.1):

    bamfile1 = pysam.Samfile(in_bam1, 'rb')
    bamfile2 = pysam.Samfile(in_bam2, 'rb')
    
    hout = open(output_file, 'w')
    with open(input_file, "r") as hin:
        for line in hin: 
            line = line.rstrip('\n')
            if line.startswith("Chr_1"):
                print(line+"\ttarget_soft_clipping1\ttarget_softclipping2\tdominant1\tdominant2\tsoftclip_in_normal2\tsoftclip_in_normal2", file=hout)
                continue  # for geomonSV results
            F = line.split('\t')

            juncChr1 = F[0]
            juncPos1 = F[1]
            juncDir1 = F[2]
            juncChr2 = F[3]
            juncPos2 = F[4]
            juncDir2 = F[5]
            juncSeq = F[6]
            target_seq1, count_junc1, count_other1, target_seq2, count_junc2, count_other2 = \
            geted(juncChr1,juncPos1,juncDir1,juncChr2,juncPos2,juncDir2,juncSeq, bamfile1, max_depth, reference_genome, validate_sequence_length, ed_threas)
            
            dominant1 = "---"
            if count_junc1+count_other1 > 0:
                diminant1 = count_junc1 / (count_junc1+count_other1)
            dominant2 = "---"
            if count_junc2+count_other2 > 0:
                dominant2 = count_junc2 / (count_junc2+count_other2)

            target_seq1_2, count_junc1_2, count_other1_2, target_seq2_2, count_junc2_2, count_other2_2 = \
            geted(juncChr1,juncPos1,juncDir1,juncChr2,juncPos2,juncDir2,juncSeq, bamfile2, max_depth, reference_genome, validate_sequence_length, ed_threas)

            print(line +"\t"+target_seq1+"\t"+target_seq2+"\t"+str(round(diminant1,4))+"\t"+str(round(dominant2,4))+"\t"+str(count_junc1_2) +"\t"+ str(count_junc2_2),file=hout)
          
    hout.close  
            

def geted(juncChr1,juncPos1,juncDir1,juncChr2,juncPos2,juncDir2,juncSeq, bamfile, max_depth, reference_genome, validate_sequence_length, ed_threas):
    # if the #sequence read is over the `maxDepth`, then that key is ignored
    depthFlag = 0
    if bamfile.count(juncChr1, int(juncPos1) - 1, int(juncPos1) + 1) >= max_depth: depthFlag = 1
    if bamfile.count(juncChr2, int(juncPos2) - 1, int(juncPos2) + 1) >= max_depth: depthFlag = 1
    if depthFlag == 1:
        print("sequence depth exceeds the threshould for: " + ','.join([juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2]), file = sys.stderr)
        return 1 

    f_inversion = False
    if (juncDir1 == "+" and juncDir2 == "+") or (juncDir1 == "-" and juncDir2 == "-"):
        f_inversion = True

    target_seq2 = getTargetSeq(juncChr2, juncPos2, juncDir2, juncSeq, reference_genome, validate_sequence_length, f_inversion)
    target_seq1 = getTargetSeq(juncChr1, juncPos1, juncDir1, juncSeq, reference_genome, validate_sequence_length, f_inversion)
                
    tmp_seq1 = ""
    tmp_ed1 = 100
    target_softclip1 = ""
    count_other1 = 0     
    count_junc1 = 0   
    for read in bamfile.fetch(juncChr1, max(0, int(juncPos1) - 1), int(juncPos1) + 1):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue
    
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")

        # M	BAM_CMATCH	0
        # I BAM_CINS	1
        # D	BAM_CDEL	2
        # S	BAM_CSOFT_CLIP	4
        # [(0, 52), (4, 24)]
        # result = 'even' if a % 2 == 0 else 'odd'
        cigars= read.cigartuples
        cigar_target = cigars[-1] if juncDir1 == "+" else cigars[0]
        clipped_size = cigar_target[1] if int(cigar_target[0]) == 4 else 0 
        
        if clipped_size < 10: continue
        
        # TODO
        # is the softclipping breakpoint within 10 bp of the junction?
        tmp_seq = ""
        if juncDir1 == "+":
            tmp_seq = read.seq[-(clipped_size):]
        else:
            tmp_seq = read.seq[:clipped_size]
        seq = ""
        # if flags[4] == "1":
        #     seq = utils.reverseComplement(str(tmp_seq))
        # else:
        seq = tmp_seq

        ret = edlib.align(seq, target_seq2, mode="HW", task="path")
        ed = ret["editDistance"]

        print("----------")
        print(ed)
        print(tmp_ed1)
        print(seq)
        print(tmp_seq1)
        if len(tmp_seq1) == 0 or ((int(ed) / len(tmp_seq1)) < tmp_ed1) or ((int(ed) / len(tmp_seq1)) == tmp_ed1 and len(seq) >= len(tmp_seq1)):
            if len(seq) >= 10:
                tmp_seq1 = seq
                tmp_ed1 = (int(ed) / len(tmp_seq1))
                
        # print (cigars)
        # print (clipped_size)
        # print (seq)
        # print (target_seq2)
        # print (ed)
        if ed != 0 and (ed/clipped_size) > ed_threas:
            count_other1 += 1
        else:
            count_junc1 += 1

    tmp_seq2 = ""
    tmp_ed2 = 100        
    count_other2 = 0
    count_junc2 = 0     
    for read in bamfile.fetch(juncChr2, max(0, int(juncPos2) - 1), int(juncPos2) + 1):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue
    
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")

        # M	BAM_CMATCH	0
        # I BAM_CINS	1
        # D	BAM_CDEL	2
        # S	BAM_CSOFT_CLIP	4
        # [(0, 52), (4, 24)]
        # result = 'even' if a % 2 == 0 else 'odd'
        cigars= read.cigartuples
        cigar_target = cigars[-1] if juncDir2 == "+" else cigars[0]
        clipped_size = cigar_target[1] if int(cigar_target[0]) == 4 else 0 

        if clipped_size < 10: continue

        # TODO
        # consider insert seq.
        # is the softclipping breakpoint within 10 bp of the junction?
        tmp_seq = ""
        if juncDir2 == "+":
            tmp_seq = read.seq[-(clipped_size):]
        else:
            tmp_seq = read.seq[:clipped_size]
        seq = ""
        # if flags[4] == "1":
        #     seq = utils.reverseComplement(str(tmp_seq))
        # else:
        seq = tmp_seq

        ret = edlib.align(seq, target_seq1, mode="HW", task="path")
        ed = ret["editDistance"]

        # print (cigars)
        # print (clipped_size)
        # print (seq)
        # print (target_seq1)
        # print (ed)
        if ed != 0 and (ed/clipped_size) > ed_threas:
            count_other2 += 1
        else:
            count_junc2 += 1

        if len(tmp_seq2) == 0 or ((int(ed) / len(tmp_seq2)) < tmp_ed2) or ((int(ed) / len(tmp_seq2)) == tmp_ed2 and len(seq) >= len(tmp_seq2)):
                tmp_seq2 = seq
                tmp_ed2 = (int(ed) / len(tmp_seq2))
                
    return([tmp_seq1, count_junc1, count_other1, tmp_seq2, count_junc2, count_other2]) 



def getTargetSeq(juncChr, juncPos, targetDir, juncSeq, reference_genome, validate_sequence_length, f_inversion):

    if juncSeq == "---": juncSeq = ""
    
    if targetDir == "+":
        seq = utils.get_seq(reference_genome, juncChr, int(juncPos) - validate_sequence_length, int(juncPos) ) 
        seq = seq + juncSeq
    else:
        seq = utils.get_seq(reference_genome, juncChr, int(juncPos), int(juncPos) + validate_sequence_length) 
        seq = juncSeq + seq

    if f_inversion == True:
        seq = utils.reverseComplement(seq)


    return seq

def curation_main(args):

    curation1(args.in_sv, args.in_bam1, args.in_bam2, args.output, args.margin, args.ref_genome, args.max_depth, args.validate_sequence_length)

