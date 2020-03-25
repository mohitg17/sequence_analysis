#!/usr/bin/python

import sys
import os

ids = {}
match = False

def main(referencefile, targetfile, targetdir):
    with open(referencefile) as f:
        for line in f:
            id = line.split("_")
            hashvalue = hash(id[0])
            ids[hashvalue] = id[0]

    with open(targetfile) as f, open(targetdir+"/"+referencefile.split(".")[0]+"_R1.fastq.gz", "w") as f2:
        count = 0
        for line in f:
            if(count == 0):
                match = False
            if(match):
                f2.write(line)
                count -= 1
            else:
                if("@") in line:
                    id = line.split()
                    hashvalue = hash(id[0])
                    if hashvalue in ids:
                        f2.write(line)
                        match = True
                        count = 3


if __name__ == "__main__":
    targetdir = sys.argv[1].split(".")[0]+"_fastqs"
    os.mkdir(targetdir)
    main(sys.argv[1], sys.argv[2], targetdir)
    

