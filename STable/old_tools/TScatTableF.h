#ifndef TScatTableF_h_seen
#define TScatTableF_h_seen

#include "TNamed.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TScatTable.h"
#include <vector>

class TScatTableF : public TNamed {

public:

  TScatTableF();
  TScatTableF(TScatTableF* tst);
  TScatTableF(TScatTable* tst);
  TScatTableF(const char *name, const char *title, int nzsbins, double minzs, double maxzs, int nrsbins, double minrs, double maxrs, int ntbins, double mint, double maxt, int nastbins, double minast, double maxast, int nctbins=1, double minct=0., double maxct=0., int nphibins=1, double minphi=0., double maxphi=0.);
  ~TScatTableF();

  bool Initialize(const char *name, const char *title, int nzsbins, double minzs, double maxzs, int nrsbins, double minrs, double maxrs, int ntbins, double mint, double maxt, int nastbins, double minast, double maxast, int nctbins=1, double minct=0., double maxct=0., int nphibins=1, double minphi=0., double maxphi=0.);
  void Reset();

  int GetIndex(int* subindices);
  int GetIndexFast(int* subindices);
  int GetIndex(int izs, int irs, int it, int iast, int ict=0, int iphi=0);
  double GetBinCenter(int idimension, int ibin);
  void GetSubIndices(int* subindices, int iglobalindex);
  bool IsInABin(double zs, double rs, double t, double ast, double ct=0., double phi=0.);
  void FindClosestSubBins(int* ibinarray, double zs, double rs, double t, double ast, double ct=0., double phi=0.);
  void FindSurroundingBins(int& min, int& max, int idim, double val);
  void SetElement(int iglobalindex, double val);
  void SetElement(int izs, int irs, int it, int iast, int ict, int iphi, double val);
  void SetElement(int izs, int irs, int it, int iast, double val);
  double GetElement(int iglobalindex);
  double GetElement(int* subindices);
  double GetElement(int izs, int irs, int it, int iast, int ict=0, int iphi=0);
  double GetElement6DFast(int* subindices);
  double GetElement6DFast(int izs, int irs, int it, int iast, int ict, int iphi);
  double GetInterpVal(double zs, double rs, double t, double ast, double ct, double phi);
  double GetInterpVal2(double zs, double rs, double t, double ast, double ct=0., double phi=0.);
  double GetInterpVal6DFast(double zs, double rs, double t, double ast, double ct=0., double phi=0.);
  double Interp(double x, double x1, double x2, double y1, double y2);
  void Fill(double zs, double rs, double t, double ast, double ct, double phi, double wgt=1.);
  void Fill(double zs, double rs, double t, double ast, double wgt=1.);
  TH1F* MakeTH1F(const char *name, const char *title, int subindex);
  TH2F* MakeTH2F(const char *name, const char *title, int subindex1, int subindex2);
  TScatTableF* Make4DScatTable();
  TScatTableF* MakeNormalized4DScatTable();

  int GetNBinsTot();
  int GetNBins(int idim);
  double GetAxisBound(int idim, int ibound);
  double GetIntegral();
  void Add(TScatTableF* tst, double c1=1.);
  void Subtract(TScatTableF* tst, int cutValue=0);
  void Divide4D(TScatTableF* tst);
  void DivideUnnormalized4D(TScatTableF* tst);

  template < typename T > void PrintVector(std::vector<T> vec);
  void PrintMinElement();
  TH1F* MakeNPhotHisto(const char* name, const char* title, int nbins, double xlow, double xhigh);
  void CheckDiff(TScatTableF* tst);
  
  void fillbininfo();

protected:

  static const int maxdimensions = 6;
  static const int maxnbins = 50;

  int nbins[maxdimensions];
  int nbinstot;
  double axisbounds[maxdimensions][2];
  //  std::vector< std::vector<int> > axisbounds(maxdimensions, std::vector<int>(2,0));
  
  std::vector<float> table;

  //  double interpvals[3][3][3][3][3][3];
  
  int gbinincr[maxdimensions];
  double bnwdth[maxdimensions];//bin width
  double binCtr[maxdimensions][maxnbins];//bin center
  double bnwdthrecp[maxdimensions];//reciprocal of bin width
  
  ClassDef(TScatTableF,2);

};

#endif
