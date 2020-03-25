#!/usr/bin/python

import sys
import os
import re

with open(sys.argv[1], 'r') as f:
    tenx_contigs = f.readlines()

with open(sys.argv[2], 'r') as f:
    pipeline_contigs = f.readlines()

new_file = "compare.fasta"

tenx_names = tenx_contigs.copy()
del tenx_names[1::2]

pipeline_names = pipeline_contigs.copy()
del pipeline_names[1::2]

del tenx_contigs[::2]
num_tenx_contigs = len(tenx_contigs)
del pipeline_contigs[::2]
num_pipeline_contigs = len(pipeline_contigs)

count = 0

with open(new_file, "w") as f:
    f.write("10x: " + str(num_tenx_contigs) + " contigs found\n" + "consensus: " + str(num_pipeline_contigs) + " contigs found\n\n")

for p in range(len(pipeline_contigs)):
    pipeline_contigs[p] = re.sub(r'([^ACTG])', "", pipeline_contigs[p])

for i in range(len(tenx_contigs)):
    tenx_contig = tenx_contigs[i]
    name = tenx_names[i]
    tenx_contig = re.sub(r'([^ACTG])', "", tenx_contig)
    mid = int(len(tenx_contig) / 2)
    start_indices = [20, (mid//4), (3*mid//8), mid-5, (5*mid//8), (3*mid)//4, len(tenx_contig)-40]
    chunks = [tenx_contig[20:30], tenx_contig[mid//4:mid//4+10], tenx_contig[(3*mid//8):(3*mid//8)+10], tenx_contig[mid-5:mid+5], tenx_contig[(5*mid//8):(5*mid//8)+10], tenx_contig[(3*mid)//4:(3*mid)//4 + 10], tenx_contig[len(tenx_contig)-40:len(tenx_contig)-30]]

    j = 0
    matches = []
    for pp_contig_num in range(len(pipeline_contigs)):
        for chunk_num in range(len(chunks)):
            if chunks[chunk_num] in pipeline_contigs[pp_contig_num]:
                matches.append([pp_contig_num, chunk_num])
    if len(matches) == 0:
        with open(new_file, "a") as f:
            f.write(name + "No matching contig" + "\n\n")
        with open("unmatched.fasta", "a") as f:
            f.write(name + tenx_contig + "\n")
        continue

    max_length = []
    for j in matches:
        pipeline_contig = pipeline_contigs[j[0]]
        offset = pipeline_contig.index(chunks[j[1]]) - start_indices[j[1]]
        start = start_indices[j[1]] + 5
        end = start_indices[j[1]] + 5

        if offset >= 0:
            while (start >= 4) and (tenx_contig[start] == pipeline_contig[start + offset] or tenx_contig[start - 4:start - 1] == pipeline_contig[start + offset - 4:start + offset - 1]):
                start -= 1
            while (end+offset < len(pipeline_contig)-4 and end+4 < len(tenx_contig)) and (tenx_contig[end] == pipeline_contig[end + offset] or tenx_contig[end + 1:end + 4] == pipeline_contig[end + offset + 1:end + offset + 4]):
                end += 1
        elif offset < 0:
            while (start+offset) >= 4 and (tenx_contig[start] == pipeline_contig[start + offset] or tenx_contig[start - 4:start - 1] == pipeline_contig[start + offset - 4:start + offset - 1]):
                start -= 1
            while (end+offset) < len(pipeline_contig)-4 and end+4 < len(tenx_contig) and (tenx_contig[end] == pipeline_contig[end + offset] or tenx_contig[end + 1:end + 4] == pipeline_contig[end + offset + 1:end + offset + 4]):
                end += 1

        if len(max_length) == 0 or (end-start) > max_length[0]:
            max_length = [end-start, start, end, j[0]]

    start = max_length[1]
    end = max_length[2]
    if (end-2-start) / len(tenx_contig) * 100 < 5:
        with open(new_file, "a") as f:
            f.write(name + "No matching contig" + "\n\n")
        with open("unmatched.fasta", "a") as f:
            f.write(name + tenx_contig + "\n")
        continue

    pipeline_contig = pipeline_contigs[max_length[3]]
    # del pipeline_contigs[j]
    count += 1
    with open(new_file, "a") as f:
        f.write(name + pipeline_names[max_length[3]])
        f.write(tenx_contig[start + 1:end - 1] + "\n")
        f.write("% match - 10x contig-->consensus: " + str((end-2-start) / len(tenx_contig) * 100) + "\n"
                + "% match - consensus-->10x contig: " + str((end-2-start) / len(pipeline_contig) * 100) + "\n\n")

with open(new_file, "a") as f:
    f.write("Overlapping contigs: " + str(count) + "\n\n")