import os, sys
import subprocess
import pysam
import edlib
from . import utils

def curation1(input_file, in_bam1, in_bam2, output_file, margin, \
            reference_genome, max_depth, \
            validate_sequence_length, validate_sequence_minus_length, \
            ed_threas):

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

            juncChr1,juncPos1,juncDir1  = F[0], F[1], F[2]
            juncChr2,juncPos2,juncDir2  = F[3], F[4], F[5]
            juncSeq = F[6]
            
            f_inversion = False
            if (juncDir1 == "+" and juncDir2 == "+") or (juncDir1 == "-" and juncDir2 == "-"):
                f_inversion = True
            
            target_seq2 = getTargetSeq(juncChr2, juncPos2, juncDir2, juncSeq, reference_genome, \
                        validate_sequence_length, validate_sequence_minus_length, f_inversion)
            target_seq1 = getTargetSeq(juncChr1, juncPos1, juncDir1, juncSeq, reference_genome, \
                        validate_sequence_length, validate_sequence_minus_length, f_inversion)
    
            ## bam1 start ############

            ret_code = checkBamDepth(bamfile1, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, max_depth)
            
            dataframe1 = getTargetDataFrame(bamfile1, juncChr1, juncPos1, juncDir1, \
                target_seq2, margin, validate_sequence_length, validate_sequence_minus_length, juncSeq)
    
            print("----------")
            
            h_group = {}
            for df in dataframe1:
                # print(df)
                if len(df["soctclipseq"]) < 10: continue
                if df["exactposition"] in h_group:
                    h_group[df["exactposition"]] += 1
                else:
                    h_group[df["exactposition"]] = 1
            
            print(h_group)        

    
            dataframe2 = getTargetDataFrame(bamfile1, juncChr2, juncPos2, juncDir2, \
                target_seq1, margin, validate_sequence_length, validate_sequence_minus_length, juncSeq)
            
            ## bam1 end ############
            
            # target_seq1_2, count_junc1_2, count_other1_2, target_seq2_2, count_junc2_2, count_other2_2 = \
            # geted(juncChr1,juncPos1,juncDir1,juncChr2,juncPos2,juncDir2,juncSeq, bamfile2, max_depth, reference_genome, validate_sequence_length, ed_threas, margin)
            # print(line +"\t"+target_seq1+"\t"+target_seq2+"\t"+str(round(dominant1,4))+"\t"+str(round(dominant2,4)),file=hout)
            # print(line +"\t"+target_seq1+"\t"+target_seq2+"\t"+str(round(dominant1,4))+"\t"+str(round(dominant2,4))+"\t"+str(count_junc1_2) +"\t"+ str(count_junc2_2),file=hout)
          
    hout.close  
            
def getTargetSeq(juncChr, juncPos, targetDir, juncSeq, reference_genome, validate_sequence_length, validate_sequence_minus_length, f_inversion):

    if juncSeq == "---": juncSeq = ""
    
    if targetDir == "+":
        seq = utils.get_seq(reference_genome, juncChr, int(juncPos) - validate_sequence_length, int(juncPos) ) 
        seq = seq + juncSeq
        seq = seq + utils.get_seq(reference_genome, juncChr, int(juncPos) + 1, int(juncPos) + validate_sequence_minus_length) 
    else:
        seq = utils.get_seq(reference_genome, juncChr, int(juncPos) - validate_sequence_minus_length, int(juncPos) -1 ) 
        seq = seq + juncSeq
        seq = seq + utils.get_seq(reference_genome, juncChr, int(juncPos), int(juncPos) + validate_sequence_length) 

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


def getTargetDataFrame(bamfile, juncChr, juncPos, juncDir, target_seq, mergin, validate_sequence_length, validate_sequence_minus_length, juncSeq):

    if juncSeq == "---": juncSeq = ""
    
    target_dataframe = []
    for read in bamfile.fetch(juncChr, max(0, int(juncPos) - mergin), int(juncPos) + mergin):

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
        
        if clipped_size == 0: continue

        softclipseq = ""
        if juncDir == "+":
            softclipseq = read.seq[-(clipped_size):]
        else:
            softclipseq = read.seq[:clipped_size]

        # softcclip pogision
        softclip_bp_pos = 0
        if juncDir == "+":
            softclip_bp_pos = read.aend
        else:
            softclip_bp_pos = read.pos + 1

        # edlib
        ret = edlib.align(softclipseq, target_seq, mode="HW", task="path")
        ed = ret["editDistance"]
        ed_locations = ret["locations"]

        # target alignment position
        alignment_start_pos = ""
        alignment_end_pos = ""
        number_of_bases_matched = -1000
        multi_blocks = False 
        # if len(ed_locations) == 1:
        ed_idx = -1
        ed_match = -1
        for index, locate in enumerate(ed_locations):
            alignment_match = ed_locations[index][1] - ed_locations[index][0]
            if alignment_match > ed_match:
                ed_idx = index
                ed_match = alignment_match
            
        if ((ed_locations[ed_idx][1] - ed_locations[ed_idx][0]) / len(softclipseq)) > 0.8:
            if juncDir == "+":
                alignment_start_pos = ed_locations[ed_idx][0] - validate_sequence_minus_length
                alignment_end_pos = ed_locations[ed_idx][1] - validate_sequence_minus_length
            else:
                alignment_start_pos = (ed_locations[ed_idx][1] - validate_sequence_length)*-1
                alignment_end_pos = (ed_locations[ed_idx][0] - validate_sequence_length)*-1
            number_of_bases_matched = alignment_end_pos - alignment_start_pos
        else:
            multi_blocks = True
    
        # exact position        
        exact_position = -1000
        if multi_blocks == False:
            if juncDir == "+":
                exact_position = int(alignment_start_pos) + (int(juncPos) - softclip_bp_pos)
            else:
                exact_position = int(alignment_start_pos) + len(juncSeq) - (int(juncPos) - softclip_bp_pos)
            
        h_ret = {"readid":read.qname, "soctclipseq":softclipseq, "direction":juncDir, "juncseq":juncSeq, \
        "junctionpos":juncPos, "softclipstart": softclip_bp_pos, "multiblocks": multi_blocks, \
        "location": (alignment_start_pos, alignment_end_pos), "cigars": cigars, "locations":ed_locations, \
        "basematch": number_of_bases_matched, "editdistance": ed, "exactposition": exact_position}
        
        target_dataframe.append(h_ret)
            
    return target_dataframe


def curation_main(args):

    curation1(args.in_sv, args.in_bam1, args.in_bam2, args.output, args.margin, \
            args.ref_genome, args.max_depth, \
            args.validate_sequence_length, args.validate_sequence_minus_length,
            args.ed_threashold)
