#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from ont_fast5_api.fast5_interface import get_fast5_file

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    exit(error_type)

def log(string):
    sys.stderr.write(f'LOG: {string}\n')

#####

modbases = None
get_all = False

# add input parsing here

mod_min_thresh = 0.5

if modbases is None:
    get_all = True

f5files = sys.argv[1:]


#####
# extract

read_ids = set()

for f5file in f5files:

    with get_fast5_file(f5file, mode="r") as f5:

        for read_id in f5.get_read_ids():

            if read_id in read_ids:
                log(f'Skipping duplicate read: {read_id}')
                continue
            else:
                read_ids.add(read_id)

            read = f5.get_read(read_id)
            latest_basecall = read.get_latest_analysis('Basecall_1D')
            mod_base_table = read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs')
            fastq_seq = read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/Fastq').split('\n')[1]

            assert mod_base_table.shape[0] == len(fastq_seq)
            n_bases = len(fastq_seq)

            table_path = '{}/BaseCalled_template/ModBaseProbs'.format(latest_basecall)
            metadata = read.get_analysis_attributes(table_path)

            assert metadata['modified_base_long_names'] == '6mA 5mC'
            assert metadata['output_alphabet'] == 'AYCZGT'

            modA = []
            modC = []

            for pos0 in range(n_bases):
                if mod_base_table[pos0, 1]/255 >= mod_min_thresh and fastq_seq[pos0] == 'A':
                    # add 1-based position in read
                    modA.append(pos0+1)
                    # log(f'{fastq_seq[pos0]}  -  {mod_base_table[pos0, :]}')
                if mod_base_table[pos0, 3]/255 >= mod_min_thresh and fastq_seq[pos0] == 'C':
                    modC.append(pos0+1)
                    # log(f'{fastq_seq[pos0]}  -  {mod_base_table[pos0, :]}')


            # output
            print(f'>{read_id}', metadata['modified_base_long_names'], 'mod_min_threshold:', mod_min_thresh)
            print(f'6mA\t{",".join([str(p) for p in modA])}')
            print(f'5mC\t{",".join([str(p) for p in modC])}')
            
            

            