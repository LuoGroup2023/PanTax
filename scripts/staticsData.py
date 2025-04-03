#!/usr/bin/env python3
import sys, os, argparse, fcntl, gzip
from tqdm import tqdm
# from prettytable import PrettyTable
import concurrent.futures
from toolkits import Logger

usage = "Stat basic information of nucleic acid sequence file in fasta format."

def main():
    parser = argparse.ArgumentParser(prog="python staticsData.py", description=usage)

    parser.add_argument("--filename", default=None, help="input a path of a fasta file ")  
    parser.add_argument("--filelist", default=None, help="input a path of a file with paths of fasta file")
    parser.add_argument("-t", "--threads", dest="threads", default=1, type=int, help="threads")
    parser.add_argument("-o", "--out", dest="out", default="genome_statics.txt", type=str, help="Output file path.")
    args = parser.parse_args()
    log = Logger()
    filename = args.filename
    filelist = args.filelist
    if filename is not None:
        single_fasta_statics(filename)
    elif filelist is not None:
        parallel_write(filelist, args.threads, args.out)
        log.logger.info("Stat basic information of nucleic acid sequence file in fasta format successfully.")
    return

def open_file(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")

def read_data(filename):
    all_sequence = {}
    sequence_name = None
    sequence = []
    with open_file(filename) as file:
        for line in file:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if sequence_name:
                    all_sequence[sequence_name] = ''.join(sequence)     
                sequence_name = line[1:]
                sequence = [] 
            else:
                sequence.append(line)
        if sequence_name:
            all_sequence[sequence_name] = ''.join(sequence)
    return all_sequence

def contig_delete_gap(sequences):
    contig_list = []
    for sequence in sequences:
        segments = sequence.split('N')
        segments_without_N = [segment for segment in segments if segment]
        contig_list.extend(segments_without_N)
    return contig_list

def scaffold_statics(sequences):
    # total number
    total_number = len(sequences)
    # max length
    max_length = max(sequences,key=len,default='')
    # min length
    min_length = min(sequences,key=len,default='')
    sequences.sort(key = lambda i:len(i),reverse=True)
        # total length
    total_length = sum(len(seq) for seq in sequences)
    # N50
    N50_length = 0
    for sequence in sequences:
        N50_length += len(sequence)
        if N50_length >= (total_length * 0.5):
            N50 = len(sequence)
            break
    # N90
    N90_length = 0
    for sequence in sequences:
        N90_length += len(sequence)
        if N90_length >= (total_length * 0.9):  
            N90 = len(sequence)  
            break
    # avg length
    avg_length = total_length / total_number
    sequence_without_N_length = sum(len(sequence.replace('N', '')) for sequence in sequences)
    # gap
    gap_length = total_length - sequence_without_N_length
    GC_sum = sum(sequence.count('G') + sequence.count('C') for sequence in sequences)
    # GC 
    GC = (GC_sum / sequence_without_N_length) * 100
    return total_number, total_length, gap_length, avg_length, N50, N90, max_length, min_length, GC

def contig_statics(sequences):
    contig_list = contig_delete_gap(sequences)
    # total number
    total_number = len(contig_list)
    # total length
    contig_without_N_length = sum(len(contig) for contig in contig_list)
    # avg length
    avg_length = contig_without_N_length / total_number
    # min length
    min_length = min(contig_list,key=len,default='')
    # max length
    max_length = max(contig_list,key=len,default='')
    contig_list.sort(key = lambda i:len(i),reverse=True)
    # N50
    N50_length = 0
    for contig in contig_list:
        N50_length += len(contig)
        if N50_length >= (contig_without_N_length * 0.5):
            N50 = len(contig)   
            break
    # N90
    N90_length = 0
    for contig in contig_list:
        N90_length += len(contig)
        if N90_length >= (contig_without_N_length * 0.9):
            N90 = len(contig) 
            break
    # GC
    GC_sum = sum(sequence.count('G') + sequence.count('C') for sequence in contig_list) 
    GC = (GC_sum/contig_without_N_length) * 100
    # gap length
    gap_length = 0
    return total_number, contig_without_N_length, gap_length, avg_length, N50, N90, max_length, min_length, GC

def contig_limi_statics(sequences):
    contig_list = contig_delete_gap(sequences)
    contig_limi = [contig for contig in contig_list if len(contig) >= 500]
    # total number
    total_number = len(contig_limi)
    # total length
    total_length = sum(len(contig) for contig in contig_limi)
    # avg length
    avg_length = total_length / total_number
    # max length
    max_length = max(contig_limi,key=len,default='')
    # min length
    min_length = min(contig_limi,key=len,default='')
    contig_limi.sort(key = lambda i:len(i), reverse=True)
    N50_length = 0
    # N50
    for contig in contig_limi:
        N50_length += len(contig)
        if N50_length >= (total_length* 0.5):
            N50= len(contig)  
            break
    # N90
    N90_length = 0
    for contig in contig_limi:
        N90_length += len(contig)
        if N90_length >= (total_length * 0.9):
            N90 = len(contig) 
            break
    # GC
    GC_sum= sum(sequence.count('G') + sequence.count('C') for sequence in contig_limi) 
    GC = (GC_sum / total_length) * 100
    # gap length
    gap_length = 0
    return total_number, total_length, gap_length, avg_length, N50, N90, max_length, min_length, GC

def original_statics(sequences, limi = False):
    sca_total_number, genome_total_length, sca_gap_length, sca_avg_length, sca_N50, sca_N90, sca_max_length, sca_min_length, sca_GC = scaffold_statics(sequences)
    contig_total_number, contig_without_N_length, contig_gap_length, contig_avg_length, contig_N50, contig_N90, contig_max_length, contig_min_length, contig_GC = contig_statics(sequences)
    contig_limi_total_number, contig_limi_total_length, contig_limi_gap_length, contig_limi_avg_length, contig_limi_N50, contig_limi_N90, contig_limi_max_length, contig_limi_min_length, contig_limi_GC = contig_limi_statics(sequences)
    x = PrettyTable()
    if limi:
        x.field_names = ['statistical level:','cutoff scaffold > (500)','contig','contig > (500)']
    else:
        x.field_names = ['statistical level:','original scaffold','contig','contig > (500)']
    x.add_row(['Total number (>)',sca_total_number,contig_total_number,contig_limi_total_number])
    x.add_row(['Total length of (bp)',genome_total_length,contig_without_N_length,contig_limi_total_length])
    x.add_row(['Gap number (bp)',sca_gap_length,contig_gap_length,contig_limi_gap_length])
    x.add_row(['Average length (bp)',('%.2f' % sca_avg_length),('%.2f' % contig_avg_length),('%.2f' % contig_limi_avg_length)])
    x.add_row(['N50 length (bp)',sca_N50,contig_N50,contig_limi_N50])
    x.add_row(['N90 length (bp)',sca_N90,contig_N90,contig_limi_N90])
    x.add_row(['Maximum length of (bp)',len(sca_max_length),len(contig_max_length),len(contig_limi_max_length)])
    x.add_row(['Minimum length of (bp)',len(sca_min_length),len(contig_min_length),len(contig_limi_min_length)])
    x.add_row(['GC content is (%)',('%.2f' % sca_GC),('%.2f' % contig_GC),('%.2f' % contig_limi_GC)])
    print(x)

def limi_statics(sequences):
    sequences = [sequence for sequence in sequences if len(sequence) >= 500]
    original_statics(sequences, limi=True)

def single_fasta_statics(filename):
    all_sequence = read_data(filename)
    sequences = list(all_sequence.values())
    original_statics(sequences)
    limi_statics(sequences)

def statics_and_write(file_path):
    all_sequence = read_data(file_path)
    sequences = list(all_sequence.values())   
    sca_total_number, genome_total_length, sca_gap_length, sca_avg_length, sca_N50, sca_N90, sca_max_length, sca_min_length, sca_GC = scaffold_statics(sequences)
    filename = os.path.basename(file_path)
    lock_file = f'{file_path}.lock'
    with open(lock_file, 'w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX)
                with open("genome_statics.txt", 'a') as file:
                    file.write(f"{filename}\t{sca_total_number}\t{genome_total_length}\t{sca_gap_length}\t{'%.2f' %sca_avg_length}\t{sca_N50}\t{sca_N90}\t{len(sca_max_length)}\t{len(sca_min_length)}\t{'%.2f' %sca_GC}\n")
            finally:
                fcntl.flock(lock, fcntl.LOCK_UN)
                os.remove(lock_file)

# def parallel_write(filelist):
#     with open(filelist, "r") as f:
#         all_paths = [line.strip() for line in f]
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         executor.map(statics_and_write, all_paths)

def statics_and_collect(file_path):
    all_sequence = read_data(file_path)
    sequences = list(all_sequence.values())   
    sca_total_number, genome_total_length, sca_gap_length, sca_avg_length, sca_N50, sca_N90, sca_max_length, sca_min_length, sca_GC = scaffold_statics(sequences)
    filename = os.path.basename(file_path)
    return f"{filename}\t{sca_total_number}\t{genome_total_length}\t{sca_gap_length}\t{'%.2f' % sca_avg_length}\t{sca_N50}\t{sca_N90}\t{len(sca_max_length)}\t{len(sca_min_length)}\t{'%.2f' % sca_GC}\n"

def parallel_write(filelist, threads, out_path):
    with open(filelist, "r") as f:
        all_paths = [line.strip() for line in f]
    
    futures = []
    chunk_size = 10000
    if len(all_paths) <= chunk_size: chunk_size = len(all_paths)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for i in tqdm(range(0, len(all_paths), chunk_size), desc="Processing chunks"):
            chunk = all_paths[i:i + chunk_size]
            futures = [executor.submit(statics_and_collect, path) for path in chunk]
            
            with open(f"{out_path}", 'a') as file:
                for future in concurrent.futures.as_completed(futures):
                    file.write(future.result())
    print("All tasks completed and results written to genome_statics.txt.")

if __name__ == "__main__":
    sys.exit(main())


    
