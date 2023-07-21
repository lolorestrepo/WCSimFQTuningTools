
// Check the tails of fitted charge PDF

int CheckcPDF(TString ifName, int iPMTType=0){
  
  gROOT->ProcessLine(".L fQChrgPDF.cc+");
  
  fQChrgPDF::Get()->LoadParams(ifName.Data());
  
  fQChrgPDF::Get()->flgStrict = true;
  
  int nPoint = 50000;
  
  double qRange[2]={0.05,1500.};
  
  cout << Form("Checking charge PDF for Type %d PMT: ",iPMTType) << endl;
  for (int i=0; i<nPoint; i++) {
    double qval = exp((log(qRange[1])-log(qRange[0]))*(double)i/(double)nPoint+log(qRange[0]));
    
//    cout << " q=" << qval << endl;
    
    int nparam;
    double par[10];
    double muthr[2];
    double coeff[6];
    
    fQChrgPDF::Get()->SetCPDFParams(qval,iPMTType,nparam,par,muthr,coeff);
  }
  
  cout << "Done!" << endl;
}
