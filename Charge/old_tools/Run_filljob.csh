#!/bin/csh
# $1: SK version [1-4]

# Run in output directory!

if ($#argv == 0) then
echo "Specify SK version!"
exit
endif

echo "SK " $1

ln -s ../detsimjob/qbins_sk$1.txt qbins.txt
ln -s ../mutbl.txt

ln -s $FITQUN_ROOT/fQChrgPDF* .

# Skip this part if you already have the 1d charge pdf files [mu]_pdf.root
foreach i (`cat mutbl.txt`)
#  ../genchrgpdf/genchrgpdf $i.dat # produce 1D normalized f(q) at each mu
end

root -b -q ../gen2d.cc # Combine 1D pdfs into 2d histogram f(q|mu)

root -b -q ../plotChrgPDF.cc # Define mu threshold as function of q

# Fit f(q|mu) wrt. mu at different q range
root -b -q ../fitpdf.cc'('0,$1')'
root -b -q ../fitpdf.cc'('1,$1')'

root -b -q ../MakecPDFparFile.cc

mv cPDFpar.root cPDFpar_sk$1.root
