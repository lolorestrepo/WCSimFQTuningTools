
// Plot the raw and the fitted charge pdf

int plotpdf(TString ifName, Double_t prjval, int fPrjax=1, int iPMTType=0){//fPrjax: Plot as a function of 0:mu 1:q
  int i,j,k;
  
  gROOT->ProcessLine(".L fQChrgPDF.cc+");
  
  fQChrgPDF::Get()->LoadParams(ifName.Data());
  
//  fQChrgPDF::Get()->flgStrict = true;
  
  TFile *parf= new TFile(ifName.Data());
  
  TH2D *hpdf2d =(TH2D*)(parf->Get(Form("hst2d_type%d",iPMTType)));
  
  TString strAxis[]={"q","mu"};
  
  int prjbin;// Bin at wich histogram is evaluated at
  if (fPrjax==1) {
    prjbin=hpdf2d->GetXaxis()->FindBin(prjval);
    
    cout << "mu=" << prjval << endl;
    cout << "Hist. at mu=" << hpdf2d->GetXaxis()->GetBinLowEdge(prjbin) << endl;
    
    TH1D *hprjct=hpdf2d->ProjectionY("proj",prjbin,prjbin,"");//function of q
  }
  else {
    prjbin=hpdf2d->GetYaxis()->FindBin(prjval);
    
    cout << "q=" << prjval << endl;
    cout << "Hist. at q=" << hpdf2d->GetYaxis()->GetBinCenter(prjbin) << endl;
    
    TH1D *hprjct=hpdf2d->ProjectionX("proj",prjbin,prjbin,"");//function of mu
  }
  
  Double_t bEdgs,bCont,bErr,*gEntX,*gEntY,*gEntYlg,*gErr,*gErru,*gErrd,Errtmp;
  
  int nbns=hprjct->GetNbinsX();
  gEntX = new Double_t[nbns];
  gEntY = new Double_t[nbns];
  gEntYlg = new Double_t[nbns];
  gErr = new Double_t[nbns];
  gErru = new Double_t[nbns];
  gErrd = new Double_t[nbns];
  
  double gEntYmax=-1.;
  
  int ngent=0;
  for (i=0; i<nbns; i++) {
    
    if (fPrjax==1) {//function of q
      bEdgs=hprjct->GetBinCenter(i+1);
    }
    else {//function of mu
      bEdgs=hprjct->GetBinLowEdge(i+1);
    }
    bCont=hprjct->GetBinContent(i+1);
    bErr=hprjct->GetBinError(i+1);
    if (bCont>0.) {
      gEntX[ngent] = bEdgs;
      gEntY[ngent] = bCont;
      if (gEntYmax<gEntY[ngent]) gEntYmax=gEntY[ngent];
      gEntYlg[ngent] = log(bCont);
      gErr[ngent] = bErr;
      gErru[ngent] = log(1+bErr/bCont);
      Errtmp = 1-bErr/bCont;
      if (!(Errtmp>0.)) Errtmp=1e-100; 
      gErrd[ngent] = -log(Errtmp);
      ngent++;
    }
  }
  TGraphAsymmErrors *glogPDF = new TGraphAsymmErrors(ngent,gEntX,gEntYlg,NULL,NULL,gErrd,gErru);
  TGraphErrors *gPDF = new TGraphErrors(ngent,gEntX,gEntY,NULL,gErr);
  
  
  if (fPrjax==1) {//function of q
    TF1 *fplot = new TF1("fPlot",Form("f2dpdf(%7.3f,x,%d)",prjval,iPMTType),0,1200);
    gPDF->GetXaxis()->SetTitle("q (p.e.)");
    gPDF->SetTitle(Form("#mu=%5.2f p.e.",prjval));
  }
  else {
    TF1 *fplot = new TF1("fPlot",Form("f2dpdf(x,%7.3f,%d)",prjval,iPMTType),0,1200);
    gPDF->GetXaxis()->SetTitle("#mu (p.e.)");
    gPDF->SetTitle(Form("q=%5.2f p.e.",prjval));
  }
  
//  cout << "Integral: " << fplot->Integral(0.01,1000) << endl;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TCanvas *c = new TCanvas("c1","",0,0,400,400);
  
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.15);
  
  cout << gEntYmax << endl;
  
  gPDF->GetXaxis()->SetTitleSize(0.06);
  gPDF->GetYaxis()->SetTitleSize(0.06);
  gPDF->GetXaxis()->SetTitleOffset(0.8);
  gPDF->GetYaxis()->SetTitleOffset(1.0);
  
  gPDF->Draw("AP");
  gPDF->GetXaxis()->SetRangeUser(0,25.);
  gPDF->GetYaxis()->SetRangeUser(0,0.7);
  
  fplot->SetNpx(1000);
  fplot->Draw("SAME");
  
  TArrow *arw = new TArrow();
  arw->SetLineWidth(1);
  
  if (fPrjax==0) {
    int iRang,nparam;
    int nRang = fQChrgPDF::Get()->GetCPDFRange(prjval,iPMTType,iRang,nparam);
    for (int k=0; k<2; k++) {
      double thrval = fQChrgPDF::Get()->gmuthr[iPMTType][iRang][k]->Eval(prjval);
      arw->DrawLine(thrval,0.,thrval,1.);
    }
  }
  
  TLatex *ltx = new TLatex;
  ltx->SetTextAlign(11);
  ltx->SetNDC();
  ltx->SetTextSize(0.08);
  
  ltx->DrawLatex(0.55,0.75,gPDF->GetTitle());
  
  gPad->SetLogy(1);
  
  c->Print(Form("cpdf_%s_%.2f_type%d.pdf",strAxis[fPrjax].Data(),prjval,iPMTType));
  
}

Double_t f2dpdf(Double_t mu, Double_t q, int iPMTType){//charge pdf
  return exp(fQChrgPDF::Get()->flogcPDF(mu,q,iPMTType));
}
