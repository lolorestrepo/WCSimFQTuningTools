#include "TScatTable.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

ClassImp(TScatTable);

//double interpvals[3][3][3][3][3][3];

TScatTable::TScatTable(){

}

TScatTable::TScatTable(TScatTable* tst){

  //  std::cout << "Inside copy constructor, (nbs,abs) have size (" << nbs.size() << "," << abs.size() << ")" << std::endl;
  TString name = tst->GetName();
  TString title = tst->GetTitle();
  name += "_2";
  title += "_2";
  Initialize(name,title,
	     tst->GetNBins(0),tst->GetAxisBound(0,0),tst->GetAxisBound(0,1),
	     tst->GetNBins(1),tst->GetAxisBound(1,0),tst->GetAxisBound(1,1),
	     tst->GetNBins(2),tst->GetAxisBound(2,0),tst->GetAxisBound(2,1),
	     tst->GetNBins(3),tst->GetAxisBound(3,0),tst->GetAxisBound(3,1),
	     tst->GetNBins(4),tst->GetAxisBound(4,0),tst->GetAxisBound(4,1),
	     tst->GetNBins(5),tst->GetAxisBound(5,0),tst->GetAxisBound(5,1));
  for (int iele=0; iele<nbinstot; iele++) {
    table[iele] = tst->GetElement(iele);
  }

}

TScatTable::TScatTable(const char *name, const char *title, int nzsbins, double minzs, double maxzs, int nrsbins, double minrs, double maxrs, int ntbins, double mint, double maxt, int nastbins, double minast, double maxast, int nctbins, double minct, double maxct, int nphibins, double minphi, double maxphi){

  Initialize(name, title, nzsbins, minzs, maxzs, nrsbins, minrs, maxrs, ntbins, mint, maxt, nastbins, minast, maxast, nctbins, minct, maxct, nphibins, minphi, maxphi);

}

TScatTable::~TScatTable(){
  
}

bool TScatTable::Initialize(const char *name, const char *title, int nzsbins, double minzs, double maxzs, int nrsbins, double minrs, double maxrs, int ntbins, double mint, double maxt, int nastbins, double minast, double maxast, int nctbins, double minct, double maxct, int nphibins, double minphi, double maxphi){

  SetName(name);
  SetTitle(title);

  bool badentries = false;
  if (nzsbins<=0) {
    std::cerr << "ERROR in TScatTable::Initialize: nzsbins = " << nzsbins << ", which is <= 0" << std::endl;
    badentries = true;
  }
  if (nrsbins<=0) {
    std::cerr << "ERROR in TScatTable::Initialize: nrsbins = " << nrsbins << ", which is <= 0" << std::endl;
    badentries = true;
  }
  if (ntbins<=0) {
    std::cerr << "ERROR in TScatTable::Initialize: ntbins = " << ntbins << ", which is <= 0" << std::endl;
    badentries = true;
  }
  if (nastbins<=0) {
    std::cerr << "ERROR in TScatTable::Initialize: nastbins = " << nastbins << ", which is <= 0" << std::endl;
    badentries = true;
  }
  if (nctbins<=0) {
    std::cerr << "ERROR in TScatTable::Initialize: nctbins = " << nctbins << ", which is <= 0" << std::endl;
    badentries = true;
  }
  if (nphibins<=0) {
    std::cerr << "ERROR in TScatTable::Initialize: nphibins = " << nphibins << ", which is <= 0" << std::endl;
    badentries = true;
  }
  if (badentries) {
    return badentries;
  }

  int idim = 0;
  nbins[idim] = nzsbins;
  axisbounds[idim][0] = minzs;
  axisbounds[idim][1] = maxzs;

  idim = 1;
  nbins[idim] = nrsbins;
  axisbounds[idim][0] = minrs;
  axisbounds[idim][1] = maxrs;

  idim = 2;
  nbins[idim] = ntbins;
  axisbounds[idim][0] = mint;
  axisbounds[idim][1] = maxt;

  idim = 3;
  nbins[idim] = nastbins;
  axisbounds[idim][0] = minast;
  axisbounds[idim][1] = maxast;

  idim = 4;
  nbins[idim] = nctbins;
  axisbounds[idim][0] = minct;
  axisbounds[idim][1] = maxct;

  idim = 5;
  nbins[idim] = nphibins;
  axisbounds[idim][0] = minphi;
  axisbounds[idim][1] = maxphi;  

  nbinstot = 1;
  for (idim=0; idim<maxdimensions; idim++) {
    if (nbins[idim] > 0) {
      nbinstot *= nbins[idim];
    }
  }
  std::cout << "Scattering table has " << nbinstot << " elements." << std::endl;
  
  fillbininfo();
  
  table.clear();
  table.resize(nbinstot,0.);
  for (int iele=0; iele<=nbinstot; iele++) {
    table[iele] = 0.;
  }

  return badentries;

}

void TScatTable::fillbininfo(){
  
  gbinincr[0]=1;
  gbinincr[1]=nbins[0];
  gbinincr[2]=nbins[0]*nbins[1];
  gbinincr[3]=nbins[0]*nbins[1]*nbins[2];
  gbinincr[4]=nbins[0]*nbins[1]*nbins[2]*nbins[3];
  gbinincr[5]=nbins[0]*nbins[1]*nbins[2]*nbins[3]*nbins[4];
  
  for (int idim=0; idim<maxdimensions; idim++) {    
    if (nbins[idim]>maxnbins) {
      std::cout << "nbins exceeds maxnbins!" << std::endl;
      exit(-1);
    }
    bnwdth[idim]=(axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
    for (int j=0; j<nbins[idim]; j++) {
      binCtr[idim][j]=axisbounds[idim][0]+bnwdth[idim]*(j+0.5);
    }
  }
  
}

void TScatTable::Reset(){
  for (int iele=0; iele<=nbinstot; iele++) {
    table[iele] = 0.;
  }  
}

int TScatTable::GetIndex(int* subindices){

  int index=0;
  for (int idim=maxdimensions-1; idim>=0; idim--) {
    if (subindices[idim] >= nbins[idim]) {
      std::cerr << "ERROR in TScatTable::GetIndex: subindex " << idim << " (" << subindices[idim] << ") is out of bounds (range = [0," << nbins[idim]-1 << "])" << std::endl;
      char a[10];
      std::cin >> a;
      return -1;
    }
    if (idim < maxdimensions-1) {
      index *= nbins[idim];
    }
    index += subindices[idim];
  }
  return index;

}

int TScatTable::GetIndexFast(int* subindices){  
  int index=0;
  for (int idim=maxdimensions-1; idim>=0; idim--) {
    index *= nbins[idim];
    index += subindices[idim];
  }
  return index;
}

int TScatTable::GetIndex(int izs, int irs, int it, int iast, int ict, int iphi){

  int subindices[maxdimensions] = {izs, irs, it, iast, ict, iphi};

  return GetIndex(subindices);

}

void TScatTable::GetSubIndices(int* subindices, int iglobalindex){

  for (int idim=0; idim<maxdimensions; idim++) {
    subindices[idim] = 0;
  }

  int prevnbins = 1;
  for (int idim=0; idim<maxdimensions; idim++) {
    if (idim > 0) {
      prevnbins *= nbins[idim-1];
    }
    subindices[idim] = ( (iglobalindex - GetIndex(subindices)) / prevnbins ) % nbins[idim];
  }

}

bool TScatTable::IsInABin(double zs, double rs, double t, double ast, double ct, double phi){

  double inputs[maxdimensions] = {zs, rs, t, ast, ct, phi};
  for (int idim=0; idim<maxdimensions; idim++) {
    if ((inputs[idim] > axisbounds[idim][1])||(inputs[idim] < axisbounds[idim][0])) {
      return false;
    }
  }

    return true;

}

double TScatTable::GetBinCenter(int idimension, int ibin) {

  if (nbins[idimension] <= 1) {
    return 0.;
  }
  return axisbounds[idimension][0] + (ibin + 0.5)*(axisbounds[idimension][1] - axisbounds[idimension][0]);

}

void TScatTable::FindClosestSubBins(int* ibinarray, double zs, double rs, double t, double ast, double ct, double phi){

  double val[maxdimensions] = {zs, rs, t, ast, ct, phi};

  for (int idim=0; idim<maxdimensions; idim++) {
    if (nbins[idim] <= 1) {
      ibinarray[idim] = 0;
      continue;
    }
    if (val[idim] >= axisbounds[idim][1]){
      ibinarray[idim] = nbins[idim] - 1;
      continue;
    } else if (val[idim] < axisbounds[idim][0]) {
      ibinarray[idim] = 0;
      continue;
    }
    double binwidth = (axisbounds[idim][1] - axisbounds[idim][0])/nbins[idim];
    for (ibinarray[idim]=0; ibinarray[idim]<nbins[idim]; ibinarray[idim]++) {
      //      std::cout << "idim, val, bin edge = " << idim << ", " << val[idim] << ", " << axisbounds[idim][0] + binwidth*(ibinarray[idim]+1) << std::endl;
      if (val[idim] < axisbounds[idim][0] + binwidth*(ibinarray[idim]+1)) {
	break;
      }
    }
  }

}

void TScatTable::FindSurroundingBins(int& min, int& max, int idim, double val){

  if (idim >= maxdimensions) {
    std::cerr << "ERROR in TScatTable::FindSurroundingBins: idim (=" << idim << ") > maxdimensions (=" << maxdimensions << ")" << std::endl;
    min = 0;
    max = 0;
    return;
  }

  if (nbins[idim] <= 1) {
    min = 0;
    min = 0;
    //    std::cout << "!!!!! FSB:  <= 1 bin, returning zeros" << std::endl;
    return;
  } else if (val < axisbounds[idim][0]) {
    min = 0;
    max = 1;
    //    std::cout << "!!!!! FSB:  val (" << val << ") < low bound (" << axisbounds[idim][0] << "), returning (0,1)" << std::endl;
    return;
  } else if (val > axisbounds[idim][1]) {
    min = nbins[idim]-2;
    max = nbins[idim]-1;
    //    std::cout << "!!!!! FSB:  val (" << val << ") > high bound (" << axisbounds[idim][1] << "), returning (" << min << "," << max << ")" << std::endl;
    return;
  }

  double binwidth = (axisbounds[idim][1] - axisbounds[idim][0])/nbins[idim];
  for (int ibin=1; ibin<nbins[idim]; ibin++) {
    double lowedge = ibin*binwidth + axisbounds[idim][0];
    min = ibin-1;
    max = ibin;
    if (val < lowedge) {
      break;
    }
  }

  //  std::cout << "!!!!! FSB:  val (" << val << "), returning (" << min << "," << max << "), which correspond to (" << min*binwidth + axisbounds[idim][0] << "," << max*binwidth + axisbounds[idim][0] << ")" << std::endl;

}

void TScatTable::SetElement(int iglobalindex, double val){

  table[iglobalindex] = val;

}

void TScatTable::SetElement(int izs, int irs, int it, int iast, int ict, int iphi, double val){

  if (table.size() > 0) {
    int index = GetIndex(izs, irs, it, iast, ict, iphi);
    if (index >= 0) {
      table[index] = val;
    }
  } else {
    std::cerr << "ERROR in TScatTable::SetElement: table has not yet been initialized" << std::endl;
  }

}

void TScatTable::SetElement(int izs, int irs, int it, int iast, double val){

  SetElement(izs, irs, it, iast, 0, 0, val);

}

double TScatTable::GetElement(int iglobalindex){

  return table[iglobalindex];

}

double TScatTable::GetElement(int izs, int irs, int it, int iast, int ict, int iphi){

  int subindices[maxdimensions] = {izs, irs, it, iast, ict, iphi};

  return GetElement(subindices);

}

double TScatTable::GetElement(int* subindices){

  if (table.size() > 0) {
    int index = GetIndex(subindices);
    if (index >= 0) {
      return table[index];
    }
  } else {
    std::cerr << "ERROR in TScatTable::GetElement: table has not yet been initialized" << std::endl;
  }

  return -1.;

}

double TScatTable::GetElement6DFast(int izs, int irs, int it, int iast, int ict, int iphi){

  int index = izs + nbins[0]*(irs + nbins[1]*(it + nbins[2]*(iast + nbins[3]*(ict + nbins[4]*(iphi)))));
  return table[index];

}

double TScatTable::GetElement6DFast(int* subindices){

  int index = subindices[0] + nbins[0]*(subindices[1] + nbins[1]*(subindices[2] + nbins[2]*(subindices[3] + nbins[3]*(subindices[4] + nbins[4]*(subindices[5])))));
  return table[index];

}

double TScatTable::GetInterpVal(double zs, double rs, double t, double ast, double ct, double phi){

  double vals[maxdimensions] = {zs, rs, t, ast, ct, phi};

//  for (int idim=0; idim<maxdimensions; idim++) {
//    if ((vals[idim] < axisbounds[idim][0])||
//	(vals[idim] > axisbounds[idim][1])) {
//      std::cout << "For idim " << idim << ", val," << vals[idim] << ", is out of bounds [" << axisbounds[idim][0]<< "," << axisbounds[idim][1] << "]" << std::endl;
//    }
//  }

  int ibins[maxdimensions];
  double binwidths[maxdimensions];
  double highedges[maxdimensions];
  for (int idim=0; idim<maxdimensions; idim++) {

    double lowbound = axisbounds[idim][0];
    double highbound = axisbounds[idim][1];
    binwidths[idim] = (highbound - lowbound) / nbins[idim];
    double lowcenter = lowbound + 0.5*binwidths[idim];
    double highcenter = highbound - 0.5*binwidths[idim];

    if (vals[idim] < lowcenter) {
      vals[idim] = lowcenter;
      ibins[idim] = 0;
      highedges[idim] = lowbound + binwidths[idim];
      continue;
    }
    if (vals[idim] > highcenter) {
      vals[idim] = highcenter;
      ibins[idim] = nbins[idim] - 1;
      highedges[idim] = highbound;
      continue;
    }

    for (ibins[idim]=0; ibins[idim]<nbins[idim]; ibins[idim]++) {
      highedges[idim] = (ibins[idim]+1)*binwidths[idim]+lowbound;
      if (vals[idim] < highedges[idim]) {
	break;
      }
    }

  }

  double closestfunc = GetElement6DFast(ibins);

  double interpsum = 0.;
  for (int idim=0; idim<maxdimensions; idim++) {

    int closestbin = ibins[idim];
    double closestcenter = highedges[idim] - 0.5*binwidths[idim];

    int nextbin;
    double nextcenter;
    if (closestbin == 0) {
      nextbin = 1;
      nextcenter = closestcenter + binwidths[idim];
    } else if (closestbin == nbins[idim] - 1) {
      nextbin = closestbin - 1;
      nextcenter = closestcenter - binwidths[idim];
    } else if (vals[idim] < closestcenter) {
      nextbin = closestbin - 1;
      nextcenter = closestcenter - binwidths[idim];
    } else {
      nextbin = closestbin + 1;
      nextcenter = closestcenter + binwidths[idim];
    }

    ibins[idim] = nextbin;
    double nextfunc = GetElement6DFast(ibins);
    ibins[idim] = closestbin;

    double interpfunc = (nextfunc - closestfunc) / (nextcenter - closestcenter) * (vals[idim] - closestcenter) + closestfunc;
    if (interpfunc < 0.) {
      std::cout << "ERROR in TScatTable::GetInterpVal(). interpfunc = " << interpfunc << std::endl;
      std::cout << "   nextfunc, closestfunc, nextcenter, closestcenter, val = " << nextfunc << ", " << closestfunc << ", " << nextcenter << ", " << closestcenter << ", " << vals[idim] << std::endl;
    }
    interpsum += interpfunc;

  }

  return interpsum / maxdimensions;

}


double TScatTable::GetInterpVal2(double zs, double rs, double t, double ast, double ct, double phi){

  double interpvals[3][3][3][3][3][3];
  
//  std::cout << std::endl << "============================================================================" << std::endl << std::endl;
//  std::cout << "zs, rs, t, ast, ct, phi = " << zs << ", " << rs << ", " << t << ", " << ast << ", " << ct << ", " << phi << std::endl;
  
  int idim = 0;
  int zsbins[2];
  FindSurroundingBins(zsbins[0],zsbins[1],idim,zs);
  double zswidth = (axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
  idim = 1;
  int rsbins[2];
  FindSurroundingBins(rsbins[0],rsbins[1],idim,rs);
  double rswidth = (axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
  idim = 2;
  int tbins[2];
  FindSurroundingBins(tbins[0],tbins[1],idim,t);
  double twidth = (axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
  idim = 3;
  int astbins[2];
  FindSurroundingBins(astbins[0],astbins[1],idim,ast);
  double astwidth = (axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
  idim = 4;
  int ctbins[2];
  FindSurroundingBins(ctbins[0],ctbins[1],idim,ct);
  double ctwidth = (axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
  idim = 5;
  int phibins[2];
  FindSurroundingBins(phibins[0],phibins[1],idim,phi);
  double phiwidth = (axisbounds[idim][1]-axisbounds[idim][0])/nbins[idim];
  
  double zsvals[2];
  double rsvals[2];
  double tvals[2];
  double astvals[2];
  double ctvals[2];
  double phivals[2];
  
  for (int izsbin=zsbins[0]; izsbin<= zsbins[1]; izsbin++) {
    idim = 0;
    int izs = izsbin - zsbins[0];
    zsvals[izs] = zswidth*izsbin+axisbounds[idim][0];
    
    for (int irsbin=rsbins[0]; irsbin<= rsbins[1]; irsbin++) {
      idim = 1;
      int irs = irsbin - rsbins[0];
      rsvals[irs] = rswidth*irsbin+axisbounds[idim][0];
      
      for (int itbin=tbins[0]; itbin<= tbins[1]; itbin++) {
        idim = 2;
        int it = itbin - tbins[0];
        tvals[it] = twidth*itbin+axisbounds[idim][0];
        
        for (int iastbin=astbins[0]; iastbin<= astbins[1]; iastbin++) {
          idim = 3;
          int iast = iastbin - astbins[0];
          astvals[iast] = astwidth*iastbin+axisbounds[idim][0];
          
          for (int ictbin=ctbins[0]; ictbin<= ctbins[1]; ictbin++) {
            idim = 4;
            int ict = ictbin - ctbins[0];
            ctvals[ict] = ctwidth*ictbin+axisbounds[idim][0];
            
            for (int iphibin=phibins[0]; iphibin<= phibins[1]; iphibin++) {
              idim = 5;
              int iphi = iphibin - phibins[0];
              phivals[iphi] = phiwidth*iphibin+axisbounds[idim][0];
              //	      interpvals[izs][irs][it][iast][ict][iphi] = GetElement6DFast(izsbin, irsbin, itbin, iastbin, ictbin, iphibin);
              interpvals[izs][irs][it][iast][ict][iphi] = table[izs + nbins[0]*(irs + nbins[1]*(it + nbins[2]*(iast + nbins[3]*(ict + nbins[4]*(iphi)))))];
            }
            interpvals[izs][irs][it][iast][ict][2] = Interp(phi, phivals[0], phivals[1], interpvals[izs][irs][it][iast][ict][0], interpvals[izs][irs][it][iast][ict][1]);
          }
          interpvals[izs][irs][it][iast][2][2] = Interp(ct, ctvals[0], ctvals[1], interpvals[izs][irs][it][iast][0][2], interpvals[izs][irs][it][iast][1][2]);
        }
        interpvals[izs][irs][it][2][2][2] = Interp(ast, astvals[0], astvals[1], interpvals[izs][irs][it][0][2][2], interpvals[izs][irs][it][1][2][2]);
      }
      interpvals[izs][irs][2][2][2][2] = Interp(t, tvals[0], tvals[1], interpvals[izs][irs][0][2][2][2], interpvals[izs][irs][1][2][2][2]);
    }
    interpvals[izs][2][2][2][2][2] = Interp(rs, rsvals[0], rsvals[1], interpvals[izs][0][2][2][2][2], interpvals[izs][1][2][2][2][2]);
  }
  interpvals[2][2][2][2][2][2] = Interp(zs, zsvals[0], zsvals[1], interpvals[0][2][2][2][2][2], interpvals[1][2][2][2][2][2]);
  
//  if (interpvals[2][2][2][2][2][2] < 0.) {
//    // Test output
//    int ibins[maxdimensions];
//    FindClosestSubBins(ibins, zs, rs, t, ast, ct, phi);
//    double inival = GetElement(ibins);
//    std::cout << "closest bin value = " << inival << ", and interp value = " << interpvals[2][2][2][2][2][2] << std::endl;
//  }
  
  return interpvals[2][2][2][2][2][2];

}

double TScatTable::GetInterpVal6DFast(double zs, double rs, double t, double ast, double ct, double phi){
  
  int binarr[6],binincrar[6];
  int iglbidx[2][2][2][2][2];
  double interpvals[3][3][3][3][3][3];
  double vals[maxdimensions] = {zs, rs, t, ast, ct, phi};
  double wt[6];
  
  for (int i=0; i<6; i++) {
    binarr[i]=int((nbins[i]-1)*(vals[i]-binCtr[i][0])/(binCtr[i][nbins[i]-1]-binCtr[i][0]));//index of closest bin center smaller than val[i]
    if (binarr[i]<0) {
      binarr[i]=0;
      binincrar[i]=0;
      wt[i]=0.;
    }
    else if (binarr[i]>nbins[i]-2) {
      binarr[i]=nbins[i]-1;
      binincrar[i]=0;
      wt[i]=0.;
    }
    else {
      binincrar[i]=gbinincr[i];
      wt[i]=(vals[i]-binCtr[i][binarr[i]])/bnwdth[i];
    }
  }
  
  iglbidx[0][0][0][0][0]=binarr[0] + nbins[0]*(binarr[1] + nbins[1]*(binarr[2] + nbins[2]*(binarr[3] + nbins[3]*(binarr[4] + nbins[4]*(binarr[5])))));
  
  
//loops here are intentionally unrolled for speed optimization!
  
  int izs,irs,it;
  
  izs=0;
  iglbidx[1][0][0][0][0]=iglbidx[0][0][0][0][0]+binincrar[0];
  do{
    irs=0;
    iglbidx[izs][1][0][0][0]=iglbidx[izs][0][0][0][0]+binincrar[1];
    do{
      it=0;
      iglbidx[izs][irs][1][0][0]=iglbidx[izs][irs][0][0][0]+binincrar[2];
      do{
        
        //---------------- ast=0 -----------------------
        interpvals[izs][irs][it][0][0][0] = table[iglbidx[izs][irs][it][0][0]];
        interpvals[izs][irs][it][0][0][1] = table[iglbidx[izs][irs][it][0][0]+binincrar[5]];
        interpvals[izs][irs][it][0][0][2] = interpvals[izs][irs][it][0][0][0]*(1.-wt[5])+interpvals[izs][irs][it][0][0][1]*wt[5];
        
        iglbidx[izs][irs][it][0][1]=iglbidx[izs][irs][it][0][0]+binincrar[4];
        interpvals[izs][irs][it][0][1][0] = table[iglbidx[izs][irs][it][0][1]];
        interpvals[izs][irs][it][0][1][1] = table[iglbidx[izs][irs][it][0][1]+binincrar[5]];
        interpvals[izs][irs][it][0][1][2] = interpvals[izs][irs][it][0][1][0]*(1.-wt[5])+interpvals[izs][irs][it][0][1][1]*wt[5];
        
        interpvals[izs][irs][it][0][2][2] = interpvals[izs][irs][it][0][0][2]*(1.-wt[4])+interpvals[izs][irs][it][0][1][2]*wt[4];
        
        
        iglbidx[izs][irs][it][1][0]=iglbidx[izs][irs][it][0][0]+binincrar[3];
        
        //---------------- ast=1 -----------------------
        interpvals[izs][irs][it][1][0][0] = table[iglbidx[izs][irs][it][1][0]];
        interpvals[izs][irs][it][1][0][1] = table[iglbidx[izs][irs][it][1][0]+binincrar[5]];
        interpvals[izs][irs][it][1][0][2] = interpvals[izs][irs][it][1][0][0]*(1.-wt[5])+interpvals[izs][irs][it][1][0][1]*wt[5];
        
        iglbidx[izs][irs][it][1][1]=iglbidx[izs][irs][it][1][0]+binincrar[4];
        interpvals[izs][irs][it][1][1][0] = table[iglbidx[izs][irs][it][1][1]];
        interpvals[izs][irs][it][1][1][1] = table[iglbidx[izs][irs][it][1][1]+binincrar[5]];
        interpvals[izs][irs][it][1][1][2] = interpvals[izs][irs][it][1][1][0]*(1.-wt[5])+interpvals[izs][irs][it][1][1][1]*wt[5];
        
        interpvals[izs][irs][it][1][2][2] = interpvals[izs][irs][it][1][0][2]*(1.-wt[4])+interpvals[izs][irs][it][1][1][2]*wt[4];
        
        
        interpvals[izs][irs][it][2][2][2] = interpvals[izs][irs][it][0][2][2]*(1.-wt[3])+interpvals[izs][irs][it][1][2][2]*wt[3];
        it++;
      }while(it<2); 
      interpvals[izs][irs][2][2][2][2] = interpvals[izs][irs][0][2][2][2]*(1.-wt[2])+interpvals[izs][irs][1][2][2][2]*wt[2];
      irs++;
    }while(irs<2); 
    interpvals[izs][2][2][2][2][2] = interpvals[izs][0][2][2][2][2]*(1.-wt[1])+interpvals[izs][1][2][2][2][2]*wt[1];
    izs++;
  }while(izs<2); 
  interpvals[2][2][2][2][2][2] = interpvals[0][2][2][2][2][2]*(1.-wt[0])+interpvals[1][2][2][2][2][2]*wt[0];
  
  return interpvals[2][2][2][2][2][2];
  
}

double TScatTable::Interp(double x, double x1, double x2, double y1, double y2){
//  double y = (y2 - y1) / (x2 - x1) * (x - x1) + y1;
//  std::cout << "(" << x1 << "," << x2 << "," << x << "," << y1 << "," << y2 << "," << y << ")" << std::endl;
  if (x < x1) {
    return y1;
  } else if (x > x2) {
    return y2;
  } else {
    return (y2 - y1) / (x2 - x1) * (x - x1) + y1;
  }
}

void TScatTable::Fill(double zs, double rs, double t, double ast, double ct, double phi, double wgt){

  double vals[maxdimensions] = {zs, rs, t, ast, ct, phi};
  for (int idim=0; idim<maxdimensions; idim++) {
    if (nbins[idim] <= 1) {
      continue;
    }
    if ((vals[idim] > axisbounds[idim][1])||(vals[idim] < axisbounds[idim][0])) {
      if (idim >= 4) {
        std::cout << "value (" << vals[idim] << ") outside boundaries [" << axisbounds[idim][0] << "," << axisbounds[idim][1] << "] in dimension " << idim << std::endl;
      }
      return;
    }
  }
  int ibins[maxdimensions];
  FindClosestSubBins(ibins, zs, rs, t, ast, ct, phi);
  int index = GetIndex(ibins);
  //  std::cout << "index = " << index << std::endl;
  if (index >= 0) {
    table[index] += wgt;
  }

}

void TScatTable::Fill(double zs, double rs, double t, double ast, double wgt){

  Fill(zs, rs, t, ast, 0, 0, wgt);

}

TH1F* TScatTable::MakeTH1F(const char *name, const char *title, int subindex){

  if ((subindex >= maxdimensions)||(subindex < 0)) {
    std::cerr << "ERROR in TScatTable::MakeTH1F: subindex (" << subindex << ") is out of bounds" << std::endl;
    return NULL;
  }
  TH1F* histo = new TH1F(name, title, nbins[subindex], axisbounds[subindex][0], axisbounds[subindex][1]);
  for (int igbin=0; igbin<nbinstot; igbin++) {
    int subindices[maxdimensions];
    GetSubIndices(subindices,igbin);
    histo->SetBinContent(subindices[subindex]+1,histo->GetBinContent(subindices[subindex]+1)+table[igbin]);
  }

  return histo;

}

TH2F* TScatTable::MakeTH2F(const char *name, const char *title, int subindex1, int subindex2){

  if ((subindex1 >= maxdimensions)||(subindex1 < 0)) {
    std::cerr << "ERROR in TScatTable::MakeTH2F: subindex1 (" << subindex1 << ") is out of bounds" << std::endl;
    return NULL;
  }
  if ((subindex2 >= maxdimensions)||(subindex2 < 0)) {
    std::cerr << "ERROR in TScatTable::MakeTH2F: subindex2 (" << subindex2 << ") is out of bounds" << std::endl;
    return NULL;
  }

  TH2F* histo = new TH2F(name, title, nbins[subindex1], axisbounds[subindex1][0], axisbounds[subindex1][1], nbins[subindex2], axisbounds[subindex2][0], axisbounds[subindex2][1]);
  for (int igbin=0; igbin<nbinstot; igbin++) {
    int subindices[maxdimensions];
    GetSubIndices(subindices,igbin);
    histo->SetBinContent(subindices[subindex1]+1,subindices[subindex2]+1,histo->GetBinContent(subindices[subindex1]+1,subindices[subindex2]+1)+table[igbin]);
  }

  return histo;

}

TScatTable* TScatTable::Make4DScatTable(){

  TString tname = GetName();
  tname += "4d";
  TString ttitle = GetTitle();
  ttitle += " 4D Projection";

  TScatTable* tst = new TScatTable(tname,ttitle,
				   nbins[0],axisbounds[0][0],axisbounds[0][1],
				   nbins[1],axisbounds[1][0],axisbounds[1][1],
				   nbins[2],axisbounds[2][0],axisbounds[2][1],
				   nbins[3],axisbounds[3][0],axisbounds[3][1]);
  for (int igbin=0; igbin<nbinstot; igbin++) {
    int sis[maxdimensions];
    GetSubIndices(sis,igbin);
    double fillval = GetElement(sis);
    double currentval = tst->GetElement(sis[0],sis[1],sis[2],sis[3]);
    tst->SetElement(sis[0],sis[1],sis[2],sis[3],currentval+fillval);
  }

  return tst;

}

TScatTable* TScatTable::MakeNormalized4DScatTable(){
  
//  int nconsolidatedbins = nbins[4] + nbins[5]; should be multiplication
  int nconsolidatedbins = nbins[4] * nbins[5];
  TScatTable* t4d = Make4DScatTable();
  for (int ibin=0; ibin<t4d->GetNBinsTot(); ibin++) {
    t4d->SetElement(ibin,t4d->GetElement(ibin)/nconsolidatedbins);
  }
  return t4d;

}

int TScatTable::GetNBinsTot(){

  return nbinstot;

}

int TScatTable::GetNBins(int idim){

  return nbins[idim];

}

double TScatTable::GetAxisBound(int idim, int ibound){

  return axisbounds[idim][ibound];

}

double TScatTable::GetIntegral(){
  double sum = 0.;
  for (int ibin=0; ibin<nbinstot; ibin++) {
    sum += table[ibin];
  }
  return sum;
}

void TScatTable::Add(TScatTable* tst, double c1){
  
  if (nbinstot != tst->GetNBinsTot()) {
    std::cerr << "ERROR in TScatTable::Add. Both scattering tables must be the same length" << std::endl;
    return;
  }
  
  for (int ibin=0; ibin<nbinstot; ibin++) {
    table[ibin] += c1*tst->GetElement(ibin);
  }
}

void TScatTable::Subtract(TScatTable* tst, int cutValue){

  if (nbinstot != tst->GetNBinsTot()) {
    std::cerr << "ERROR in TScatTable::Subtract. Both scattering tables must be the same length" << std::endl;
    return;
  }

  int nmoredirect = 0;
  for (int ibin=0; ibin<nbinstot; ibin++) {
    if ((table[ibin] < cutValue)||
	(tst->GetElement(ibin) < cutValue)) {
      table[ibin] = 0.;
      continue;
    }
    if (table[ibin] < tst->GetElement(ibin)) {
      //      std::cout << "Strange: In TScatTable::Subtract (a - b), a = " << table[ibin] << " and b = " << tst->GetElement(ibin) << std::endl;
      nmoredirect++;
      table[ibin] = 0.;
      continue;
    }

    table[ibin] -= tst->GetElement(ibin);

  }
  std::cout << "In TScatTable::Subtract, found more direct light than scattered light " << nmoredirect << " times." << std::endl;

}

void TScatTable::Divide4D(TScatTable* tst){

  for (int ibin=0; ibin<nbinstot; ibin++) {
    int sis[maxdimensions];
    GetSubIndices(sis,ibin);
    double denom = tst->GetElement(sis[0],sis[1],sis[2],sis[3]);
    if (denom <= 0.) {
      table[ibin] = 0.;
      continue;
    }
    table[ibin] /= denom;
  }

}

void TScatTable::DivideUnnormalized4D(TScatTable* tst){
  
  int nconsolidatedbins = nbins[4] * nbins[5];
  
  for (int ibin=0; ibin<nbinstot; ibin++) {
    int sis[maxdimensions];
    GetSubIndices(sis,ibin);
    double denom = tst->GetElement(sis[0],sis[1],sis[2],sis[3])/nconsolidatedbins;
    if (denom <= 0.) {
      table[ibin] = 0.;
      continue;
    }
    table[ibin] /= denom;
  }
  
}


template < typename T > void TScatTable::PrintVector(std::vector<T> vec){

  std::cout << "(";
  int vsize = vec.size();
  for (int iele=0; iele<vsize; iele++) {
    std::cout << vec[iele];
    if (iele < vsize-1) {
      std::cout << ",";
    }
  }
  std::cout << ")";

}

void TScatTable::PrintMinElement(){
  double minelement = 1.e10;
  int counter = 0;
  for (int iele=0; iele<nbinstot; iele++) {
    //    std::cout << "table[" << iele << "] = " << table[iele] << std::endl;
    if (table[iele] <= minelement) {
      if (fabs(table[iele] - minelement)<1.e-5) {
	counter++;
      } else {
	minelement = table[iele];
	counter = 1;
      }
    }
  }
  std::cout << "Found value of " << minelement << " a total of " << counter << " times." << std::endl;
}

TH1F* TScatTable::MakeNPhotHisto(const char* name, const char* title, int nxbins, double xlow, double xhigh) {

  TH1F* hnphot = new TH1F(name,title,nxbins,xlow,xhigh);
  for (int ibin=0; ibin<nbinstot; ibin++) {
    hnphot->Fill(table[ibin]);
  }
  return hnphot;

}

void TScatTable::CheckDiff(TScatTable* tst){

  if (nbinstot != tst->GetNBinsTot()){
    std::cout << "These 2 tables have a different number of bins" << std::endl;
    return;
  }
  for (int ibin=0; ibin<nbinstot; ibin++) {
    if (fabs(table[ibin]-tst->GetElement(ibin)) > 1.e-5) {
      int sis[maxdimensions];
      GetSubIndices(sis,ibin);
      std::cout << "Index " << ibin << " (" << sis[0] << "," << sis[1] << "," << sis[2] << "," << sis[3] << "," << sis[4] << "," << sis[5] << ") has values of " << table[ibin] << " and " << tst->GetElement(ibin) << std::endl;
    }
  }

}
