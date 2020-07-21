import os, sys
import subprocess
import pysam
import edlib
from . import utils

def curation1(input_file, in_bam1, in_bam2, output_file, margin, \
            reference_genome, max_depth, \
            validate_sequence_length, validate_sequence_minus_length, \
            ed_threas, min_mapping_quality, f_bedpe):

    bamfile1 = pysam.Samfile(in_bam1, 'rb')
    if in_bam2 != None:
        bamfile2 = pysam.Samfile(in_bam2, 'rb')
    
    hout = open(output_file, 'w')
    with open(input_file, "r") as hin:
        for line in hin: 
            line = line.rstrip('\n')
            if line.startswith("Chr_1"):
                print(line+"\ttarget_softclipping1\ttarget_softclipping2\tdominant1\tdominant2\tsoftclip_in_tumor1\tsoftclip_in_tumor2\tsoftclip_in_normal1\tsoftclip_in_normal2", file=hout)
                continue  # for geomonSV results
            F = line.split('\t')

            juncChr1,juncPos1,juncDir1  = F[0], F[1], F[2]
            juncChr2,juncPos2,juncDir2  = F[3], F[4], F[5]
            if f_bedpe == True:
                juncChr1,juncPos1,juncDir1  = F[0], F[2], F[9]
                juncChr2,juncPos2,juncDir2  = F[3], F[5], F[10]
            juncSeq = F[6]
            
            f_inversion = False
            if (juncDir1 == "+" and juncDir2 == "+") or (juncDir1 == "-" and juncDir2 == "-"):
                f_inversion = True
            
            target_seq2 = getTargetSeq(juncChr2, juncPos2, juncDir2, juncSeq, reference_genome, \
                        validate_sequence_length, validate_sequence_minus_length, f_inversion)
            target_seq1 = getTargetSeq(juncChr1, juncPos1, juncDir1, juncSeq, reference_genome, \
                        validate_sequence_length, validate_sequence_minus_length, f_inversion)
    
            ret_code = checkBamDepth(bamfile1, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, max_depth)
            
            dataframe1 = getTargetDataFrame(bamfile1, juncChr1, juncPos1, juncDir1, \
                target_seq2, margin, validate_sequence_length, validate_sequence_minus_length, juncSeq, min_mapping_quality)
    
            dominant_group1, dominant1, target_softclip1= getDominantGroup(dataframe1, ed_threas)

            target_count1 = getDominantReadCount(dataframe1, dominant_group1, ed_threas)
            
            dataframe2 = getTargetDataFrame(bamfile1, juncChr2, juncPos2, juncDir2, \
                target_seq1, margin, validate_sequence_length, validate_sequence_minus_length, juncSeq, min_mapping_quality)

            dominant_group2, dominant2, target_softclip2 = getDominantGroup(dataframe2, ed_threas)

            target_count2 = getDominantReadCount(dataframe2, dominant_group2, ed_threas)
            
            
            # print("----------------------------------------------")
            # print(juncChr1)
            # print(juncPos1)            
            # print(juncDir1)            
            # print("dominant_group1: "+str(dominant_group1))
            # print("dominant1: "+str(dominant1))
            # print("target_softclip1: "+str(target_softclip1))
            
            # print("----------------------------------------------")
            # print(juncChr2)
            # print(juncPos2)
            # print(juncDir2)
            # print("dominant_group2: "+str(dominant_group2))
            # print("dominant2: "+str(dominant2))
            # print("target_softclip2: "+str(target_softclip2))

            if in_bam2 != None:
                ret_code = checkBamDepth(bamfile2, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, max_depth)

                dataframe1_2 = getTargetDataFrame(bamfile2, juncChr1, juncPos1, juncDir1, \
                    target_seq2, margin, validate_sequence_length, validate_sequence_minus_length, juncSeq, min_mapping_quality)
            
                target_count1_2 = getDominantReadCount(dataframe1_2, dominant_group1, ed_threas)
            
                dataframe2_2 = getTargetDataFrame(bamfile2, juncChr2, juncPos2, juncDir2, \
                    target_seq1, margin, validate_sequence_length, validate_sequence_minus_length, juncSeq, min_mapping_quality)
                
                target_count2_2 = getDominantReadCount(dataframe2_2, dominant_group2, ed_threas)
                
            # print(line +"\t"+target_seq1+"\t"+target_seq2+"\t"+str(round(dominant1,4))+"\t"+str(round(dominant2,4)),file=hout)
            print_line = line +"\t"+target_softclip1+"\t"+target_softclip2+"\t"+str(round(dominant1,4))+"\t"+str(round(dominant2,4))+"\t"+str(target_count1) +"\t"+str(target_count2)
           
            if in_bam2 != None:
                print_line = print_line +"\t"+str(target_count1_2) +"\t"+ str(target_count2_2)

            print(print_line, file=hout)
          
    hout.close  
            
def getTargetSeq(juncChr, juncPos, targetDir, juncSeq, reference_genome, validate_sequence_length, validate_sequence_minus_length, f_inversion):

    if juncSeq == "---": juncSeq = ""
    
    if targetDir == "+":
        seq = utils.get_seq(reference_genome, juncChr, int(juncPos) - validate_sequence_length + 1 , int(juncPos) ) 
        seq = seq + juncSeq
        seq = seq + utils.get_seq(reference_genome, juncChr, int(juncPos) + 1, int(juncPos) + validate_sequence_minus_length) 
    else:
        seq = utils.get_seq(reference_genome, juncChr, int(juncPos) - validate_sequence_minus_length , int(juncPos) -1 ) 
        seq = seq + juncSeq
        seq = seq + utils.get_seq(reference_genome, juncChr, int(juncPos), int(juncPos) + validate_sequence_length - 1) 
    if f_inversion == True:
        seq = utils.reverseComplement(seq)

    return seq


def checkBamDepth(bamfile, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, max_depth):

    # if the #sequence read is over the `maxDepth`, then that key is ignored
    depthFlag = 0
    if bamfile.count(juncChr1, int(juncPos1) - 1, int(juncPos1) + 1) >= max_depth: depthFlag = 1
    if bamfile.count(juncChr2, int(juncPos2) - 1, int(juncPos2) + 1) >= max_depth: depthFlag = 1
    if depthFlag == 1:
        print("sequence depth exceeds the threshould for: " + ','.join([juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2]), file = sys.stderr)
        return 1 


def getDominantGroup(dataframe, ed_threas):
    
    h_group = {}
    for df in dataframe:
        exact_pos = df["exactposition"]
        softclip_len = len(df["softclipseq"])
        basematch = df["basematch"]
        if (basematch / softclip_len) < ed_threas:
            exact_pos = exact_pos-1000 
        if exact_pos in h_group:
            h_group[exact_pos] += 1
        else:
            h_group[exact_pos] = 1
    
    target_group = "---"
    target_value = 0
    other_value = 0
    f_target = True
    for key, value in sorted(h_group.items(), key=lambda x: -x[1]):
        if  f_target == True:
            target_group = key
            target_value = int(value)
            f_target = False
        else:
            other_value += int(value)

    dominant = 0
    if target_value > 0 or other_value > 0:
        dominant = target_value / (target_value + other_value)  

    ret_seq = "---"    
    for df in dataframe:
        if df["exactposition"] == target_group:
            softclip_len = len(df["softclipseq"])
            basematch = df["basematch"]
            if (basematch / softclip_len) < ed_threas: continue
        
            if len(df["softclipseq"]) > len(ret_seq):
                ret_seq = df["softclipseq"]

    return([target_group, dominant, ret_seq])


def getDominantReadCount(dataframe, target_group, ed_threas):
    
    h_group = {}
    for df in dataframe:
        exact_pos = df["exactposition"]
        softclip_len = len(df["softclipseq"])
        basematch = df["basematch"]
        if (basematch / softclip_len) < ed_threas:
            exact_pos = exact_pos-1000 
        if exact_pos in h_group:
            h_group[exact_pos] += 1
        else:
            h_group[exact_pos] = 1
    
    read_count = 0
    if target_group in h_group:
        read_count = h_group[target_group]

    return read_count


def getTargetDataFrame(bamfile, juncChr, juncPos, juncDir, target_seq, mergin, validate_sequence_length, validate_sequence_minus_length, juncSeq, min_mapping_quality):

    if juncSeq == "---": juncSeq = ""
    
    target_dataframe = []
    for read in bamfile.fetch(juncChr, max(0, int(juncPos) - mergin), int(juncPos) + mergin):

        # get mapq
        if read.mapq <= min_mapping_quality: continue

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue
    
        # softclip sequence 
        # M	BAM_CMATCH	0
        # I BAM_CINS	1
        # D	BAM_CDEL	2
        # S	BAM_CSOFT_CLIP	4
        # [(0, 52), (4, 24)]
        # result = 'even' if a % 2 == 0 else 'odd'
        cigars= read.cigartuples
        cigar_target = cigars[-1] if juncDir == "+" else cigars[0]
        clipped_size = cigar_target[1] if int(cigar_target[0]) == 4 else 0 
        
        if clipped_size < 10: continue
    
        # softclip seq, softclip pos
        softclipseq = ""
        softclip_bp_pos = 0
        if juncDir == "+":
            softclipseq = read.seq[-(clipped_size):]
            softclip_bp_pos = read.aend
        else:
            softclipseq = read.seq[:clipped_size]
            softclip_bp_pos = read.pos + 1

        if abs(softclip_bp_pos - int(juncPos)) > mergin: continue 

        # edlib
        ret = edlib.align(softclipseq, target_seq, mode="HW", task="path")
        ed = ret["editDistance"]
        ed_locat = ret["locations"]
        # print(ret)
        # nice = edlib.getNiceAlignment(ret, softclipseq, target_seq)
        # print(nice)
        # target alignment position
        ed_idx = -1
        for tmp_index, locate in enumerate(ed_locat):
            if len(softclipseq) == (ed_locat[tmp_index][1] - ed_locat[tmp_index][0] + 1):
                ed_idx = tmp_index
                break
                
        # exact position        
        exact_position = 0
        if juncDir == "+":
            alignment_start_pos = ed_locat[ed_idx][0] - validate_sequence_minus_length
            alignment_end_pos = ed_locat[ed_idx][1] - validate_sequence_minus_length
            exact_position = int(alignment_start_pos) + (int(juncPos) - softclip_bp_pos)
        else:
            alignment_start_pos = (ed_locat[ed_idx][1] - (validate_sequence_length - 1)) * -1
            alignment_end_pos = (ed_locat[ed_idx][0] - (validate_sequence_length - 1)) * -1 
            exact_position = int(alignment_start_pos) + len(juncSeq) - (int(juncPos) - softclip_bp_pos)
            
        # data frame
        h_ret = {"readid":read.qname, "softclipseq":softclipseq, "direction":juncDir, "juncseq":juncSeq, \
        "junctionpos":juncPos, "softclipstart": softclip_bp_pos,  \
        "basematch": len(softclipseq) - ed, "editdistance": ed, "exactposition": exact_position}
        
        target_dataframe.append(h_ret)
            
    return target_dataframe


def getColorSupportRead(bamfile, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, mergin):

    colors_read = 0
    for read in bamfile.fetch(juncChr1, max(0, int(juncPos1) - mergin), int(juncPos1) + mergin):

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

        if juncChr2 == chr_pair and juncDir2 == dir_pair:
            if dir_pair == "+":
                if juncPos2 - 1000 < pos_pair and pos_pair < juncPos2 + 20:
                    colors_read += 1
            else:
                if juncPos2 - 20 < pos_pair and pos_pair < juncPos2 + 1000:
                    colors_read += 1

def curation_main(args):

    curation1(args.in_sv, args.in_bam1, args.in_bam2, args.output, args.margin, \
            args.ref_genome, args.max_depth, \
            args.validate_sequence_length, args.validate_sequence_minus_length,
            args.ed_threashold, args.min_mapping_quality, args.bedpe)

