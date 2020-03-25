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
    mid = (int)(len(tenx_contig) / 2)
    chunk = tenx_contig[mid - 10:mid + 10]

    j = 0
    matches = []
    for k in range(len(pipeline_contigs)):
        if chunk in pipeline_contigs[k]:
            matches.append(k)
    if len(matches) == 0:
        with open(new_file, "a") as f:
            f.write(name + "No matching contig" + "\n\n")
        with open("unmatched.fasta", "a") as f:
            f.write(name + tenx_contig + "\n")
        continue

    max_length = []
    for j in matches:
        pipeline_contig = pipeline_contigs[j]
        offset = pipeline_contig.index(chunk) - mid + 10
        start = mid
        end = mid

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
            max_length = [end-start, start, end, j]

    start = max_length[1]
    end = max_length[2]
    if (end-2-start) / len(tenx_contig) * 100 < 20:
        with open(new_file, "a") as f:
            f.write(name + "No matching contig" + "\n\n")
        continue

    pipeline_contig = pipeline_contigs[max_length[3]]
    # del pipeline_contigs[j]
    count += 1
    with open(new_file, "a") as f:
        f.write(name + pipeline_names[max_length[3]])
        f.write(tenx_contig[start + 1:end - 2] + "\n")
        f.write("% match - 10x contig-->consensus: " + str((end-2-start) / len(tenx_contig) * 100) + "\n"
                + "% match - consensus-->10x contig: " + str((end-2-start) / len(pipeline_contig) * 100) + "\n\n")

with open(new_file, "a") as f:
    f.write("Overlapping contigs: " + str(count) + "\n\n")