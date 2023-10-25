import ROOT
import re
import glob
import numpy as np
from os.path import exists, expandvars, join, basename


ROOT.gSystem.Load("/pbs/home/g/gdiazlop/Software/HK_Software/WCSim/install-Linux_x86_64-gcc_9-python_3.10.13/lib/libWCSimRoot.so")

if exists("makeChargePDFplot_C.so"): ROOT.gSystem.Load("makeChargePDFplot_C.so")
else: ROOT.gROOT.LoadMacro("makeChargePDFplot.C++")

get_mu_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))
infiles = glob.glob(join(expandvars("$LUSTRE/Charge/out"), "*.root"))
infiles = [filename for filename in infiles if re.match("out_\d+(?:\.\d+)?_1.root", basename(filename))]
infiles = sorted(infiles, key=get_mu_and_index)

# mus = np.array([get_mu_and_index(f)[0] for f in infiles])
# for mu in mus: 
#     if mu.is_integer(): mu = int(mu)
#     print(mu)

for file in infiles:
    mu = get_mu_and_index(file)[0]
    if mu.is_integer(): mu = int(mu)

    print(f"Processing mu={mu}")

    ROOT.makeChargePDFplot(file, f"./chargePDF/{mu}_pdf.root")