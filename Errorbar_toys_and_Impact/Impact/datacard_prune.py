#!/usr/bin/env python3
import re

dc = open('../hzg_datacard_v1p4p0_cleaned.txt','r')
outfile = open('sys_list.txt','w')
for line in dc:
    nentry = 0
    end = 0
    start = 0
    if len(line) < 1000:
        continue
    for i in range(len(line)):
        if not line[i] == ' ' and line[i+1] == ' ':
            outfile.write(line[:i+1] + '\n')
            break
dc.close()
outfile.close()
