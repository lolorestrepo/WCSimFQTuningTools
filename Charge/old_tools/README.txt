
makeChargePDFplot.C:

  gSystem->Load("/pbs/home/g/gdiazlop/Software/WCSim/install/lib/libWCSimRoot.so");
  .L makeChargePDFplot.C++
  makeChargePDFplot("../../Simulation/local/out/out.root", "ChargePDF.root")


gen2d.cc:
  .L gen2d.cc++
  gen2d()

fitpdf.cc:
  .L /pbs/home/g/gdiazlop/Software/fiTQun/fiTQun/fQChrgPDF.cc++
  .L /pbs/home/g/gdiazlop/Software/fiTQun/fiTQun/fQChrgPDF_cc.so
  .L fitpdf.cc++
