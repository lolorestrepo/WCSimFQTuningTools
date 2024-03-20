#!/usr/bin/python

import os
if os.environ['TERM'] == 'xterm':
    os.environ['TERM'] = 'linux'
# Now it's OK to import readline :)

import sys
from ROOT import TFile, TTree, gROOT
def main() :

    gROOT.ProcessLine("gErrorIgnoreLevel = 3000;")
    inFile = TFile(sys.argv[1])
    inTree = inFile.Get("wcsimT")
    entries = inTree.GetEntries()
    inFile.Close()
    print str(entries).strip()
    return entries

if __name__ == '__main__' :
    main()
