#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .curation_sv import curation_main 

def create_parser():
    prog = "bateira_sv"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
    
    def _create_curate_parser(subparsers):
        
        curate_parser = subparsers.add_parser("curation", help = "Curate the SV results")
        curate_parser.add_argument("--in_sv", help = "The result of SV", type = str, required=True)
        curate_parser.add_argument("--in_bam1", help = "The first of the bum files (usually a tumor)", type = str, required=True)
        curate_parser.add_argument("--in_bam2", help = "The second of the bum files (usually a normal)", type = str, required=True)
        curate_parser.add_argument("--output", help = "The output file", type = str, required=True)
        curate_parser.add_argument("--margin", help = "The margin for comparing SVs and SVs", type = int, default = 3)
        curate_parser.add_argument("--ref_genome", help = "The path to the reference genomoe sequence", type=str)
        curate_parser.add_argument("--max_depth", help = "Candidates having coverages more than specified value are ignored (default: %(default)s)", type = int, default = 1000)
        curate_parser.add_argument("--validate_sequence_length", help = "Length of sequences for validation (each from the breakpoint) (default: %(default)s)", type = int, default = 200)
        curate_parser.add_argument("--validate_sequence_minus_length", help = "Length of sequences for validation (each from the breakpoint) (default: %(default)s)", type = int, default = 20)
        curate_parser.add_argument("--ed_threashold", help = "The threashols of the edit distance", type = float, default = 0.05)
   
        return curate_parser

    curate_parser = _create_curate_parser(subparsers)
    curate_parser.set_defaults(func = curation_main)
    return parser
