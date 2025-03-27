#!/bin/env python3

import sys
import re
import subprocess
import glob

pref = sys.argv[1]

count=0
curchart = 1
mxcount = 27


argmat0 = ["/programs/Krona-2.8/bin/ktImportText","-o"]

argmat =argmat0.copy()
curname = pref + "_krona"+str(curchart)+".html"
argmat.append(curname)

argline = " -o " + curname

for fle in sorted(glob.glob("*.report.krona")):
  count += 1
  if count > mxcount:
    print(argline)
    subprocess.run(argmat,stdout=subprocess.PIPE)
    count = 0
    curchart += 1
    argmat =argmat0.copy()
    curname = pref + "_krona"+str(curchart)+".html"
    argmat.append(curname)
    argline = " -o " + curname

  print(fle)
  sample = re.split(r'_',fle)[0]
  print(sample)
  argline += ' ' + fle + ',' + sample
  argmat.append(fle + ',' + sample)

print(argline)
subprocess.run(argmat,stdout=subprocess.PIPE)