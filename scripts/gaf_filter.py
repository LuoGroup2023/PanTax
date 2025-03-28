#!/usr/bin/env python3
from collections import defaultdict
from read_classification import ReadClassification
import argparse, sys, json

def filter_max_alignment(gaf_file):
    readid_to_max_alignment_10 = defaultdict(float)
    readid_to_max_alignment_16 = defaultdict(float)

    with open(gaf_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) >= 16:
                read_id = fields[0]
                alignment_value_10 = int(fields[9])             

                alignment_field_16 = fields[15]
                alignment_value_16 = float(alignment_field_16.split(':')[-1])

                if alignment_value_10 > readid_to_max_alignment_10[read_id]:
                    readid_to_max_alignment_10[read_id] = alignment_value_10
                    readid_to_max_alignment_16[read_id] = alignment_value_16
                elif alignment_value_10 == readid_to_max_alignment_10[read_id] and alignment_value_16 > readid_to_max_alignment_16[read_id]:
                    readid_to_max_alignment_16[read_id] = alignment_value_16

    output_file = gaf_file.replace('.gaf', '_filtered.gaf')
    unique_read_ids = set()
    with open(gaf_file, 'r') as file, open(output_file, 'w') as output:
        for line_number,line in enumerate(file, start=1):
            fields = line.strip().split('\t')
            read_id = fields[0]
            alignment_value_10 = int(fields[9])

            alignment_field_16 = fields[15]
            alignment_value_16 = float(alignment_field_16.split(':')[-1])

            #awk -F '\t' '$12>20 && $4-$3>1000' $gaf_file > $gaf_baseName"_f.gaf"
            if int(fields[11]) > 20 and (int(fields[3]) - int(fields[2])) > 1000:
                if alignment_value_10 == readid_to_max_alignment_10[read_id] and alignment_value_16 == readid_to_max_alignment_16[read_id]:
                    if read_id not in unique_read_ids:
                        output.write(line)
                        unique_read_ids.add(read_id)

    print(f"Filtered GAF file written to: {output_file}")

def json_filter(gaf_file, json_file, species_range_file):
    read_cls = ReadClassification(species_range_file, gaf_file)
    species = read_cls.read_species_range_file()
    reads = read_cls.read_mapped_gaf_file()
    results = read_cls.thread_parallel_batch_read_classification(reads)
    reads_info = {}
    for result in results:
        reads_info[result[0]] = result[2]
    filter_json_file = json_file.replace(".json", "_filtered.json")
    
    with open(json_file, "r") as fr, open(filter_json_file, "w") as fw:
        for line in fr:
            aln = json.loads(line)
            try:
                read_id = aln["name"]
                path = aln["path"]
                mapping = path["mapping"]
                read_info = reads_info[read_id]
                if len(aln["sequence"]) < 200:
                    continue
            except KeyError:
                # read unmapped
                continue
            for node_info in mapping:
                position = node_info["position"]
                node_id = int(position["name"])
                if node_id in read_info:
                    fw.write(line)
                    del reads_info[read_id]
                    break

def main():
    parser = argparse.ArgumentParser(prog='gaf_filter.py', description='Calculate and output read lengths.')
    parser.add_argument('-i', '--gaf_file', dest="gaf_file", help="gaf file from graphAligner", required=True)
    parser.add_argument('-j', '--json_file', dest="json_file", help="json file from graphAligner", default=None)
    parser.add_argument("-s", "--species_range_file", dest="species_range_file", help="Species range file", default=None)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    if not args.json_file:
        filter_max_alignment(args.gaf_file)
    else:
        json_filter(args.gaf_file, args.json_file, args.species_range_file)


if __name__ == "__main__":
    sys.exit(main())
