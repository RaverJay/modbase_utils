#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pathlib
import argparse
from ont_fast5_api.fast5_interface import get_fast5_file

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    exit(error_type)

def log(string):
    sys.stderr.write(f'LOG: {string}\n')

#####

parser = argparse.ArgumentParser(description='Extract modified bases from guppy fast5 basecalls.\nThis only applies if --modbase was enabled for basecalling.',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("fast5dir", metavar='DIR', help="folder containing fast5 files with modbase basecalls")
parser.add_argument("-t", "--threshold_min_mod", metavar='FLOAT', type=float, help="minimum caller confidence of modified bases that will be extracted", default=0.5)
parser.add_argument("-m", "--modbases_string", metavar='STRING', help="modified base(s) to extract\nneed to be comma-separated strings with structure:  name,symbol,canonicalbase",
                    action='append', default=["6mA,Y,A", "5mC,Z,C"])
args = parser.parse_args()

modbases = {}
for mb in args.modbases_string:
    mbname, mbsymbol, mbcanon = mb.split(',')
    assert mbname not in modbases
    modbases[mbname] = {'symbol': mbsymbol, 'canonical': mbcanon}
    log(f'Modified base: {mbname} {mbsymbol} {mbcanon}')   

#####
# extract

read_ids = set()

f5path = pathlib.Path(args.fast5dir)

for f5file in f5path.rglob('*.fast5'):

    with get_fast5_file(str(f5file), mode="r") as f5:

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

            # gather data
            modpos = {}
            modcol = {}
            for mb in modbases:
                assert mb in metadata['modified_base_long_names']
                assert modbases[mb]['symbol'] in metadata['output_alphabet']
                modcol[mb] = metadata['output_alphabet'].index(modbases[mb]['symbol'])
                modpos[mb] = []

            for pos0 in range(n_bases):
                for mb in modbases:
                    if fastq_seq[pos0] == modbases[mb]['canonical'] and mod_base_table[pos0, modcol[mb]]/255 >= args.threshold_min_mod:
                        # add 1-based position in read
                        modpos[mb].append(pos0+1)


            # output
            print(f'>{read_id}', ','.join(modbases.keys()), f'mod_min_threshold:{args.threshold_min_mod}')
            for mb in modbases:
                print(f'{mb}\t{",".join([str(p) for p in modpos[mb]])}')
