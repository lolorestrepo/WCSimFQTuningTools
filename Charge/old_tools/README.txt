
makeChargePDFplot.C:

  gSystem->Load("/pbs/home/g/gdiazlop/Software/WCSim/install/lib/libWCSimRoot.so");
  .L makeChargePDFplot.C++
  makeChargePDFplot("../../Simulation/local/out/out.root", "ChargePDF.root")


gen2d.cc:
  .L gen2d.cc
  gen2d()

fitpdf.cc:
  .L /pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/fQChrgPDF.cc++
  .L /pbs/home/g/gdiazlop/Software/HK_Software/fiTQun/fQChrgPDF_cc.so
  .L fitpdf.cc
  fitpdf(0, 1)
