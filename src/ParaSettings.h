// File with small helper classes
#ifndef ParaSettingsH
#define ParaSettingsH

#include <map>
#include <vector>
#include <string>
#include <cmath>

//---------------------------------------------------------------------------
class CSpecSquare
{
public:
	int iX;
	int iY;
	std::map<int,int> Spec;

	CSpecSquare(){};

	~CSpecSquare() {Spec.clear();};

	void InitspecSquare(int x, int y)
	{
		iX=x;
		iY=y;
		Spec.clear();
	};
};

//--------------------------------------------------------------------------------
//model settings that apply to several simulation runs and are usually not varied
//in parameter optimization
class CModelSettings
{
public:
   //technical settings
	int m_nGen;         //number generations (# complete turnover of community)
	bool m_stepsOut;
	double m_cellSize;  // size of neighborhood grid

   int m_Jm;           //metacommunity size in number of individuals
   int m_metaSAD;      //model for metacommunity
                       //0 ... uniform (metaSR)
                       //1 ... lognormal (metaSR, metaCV)

   //local community size - total extent of the simulated forest
   int m_nTrees;           //number of trees in local community
   double m_Xext;
   double m_Yext;

   //point pattern output
   double m_rmax;      //maximum neighborhood radius for point pattern
   double m_bw1;     //bandwidth for community level patterns

   CModelSettings();
   CModelSettings(int nGen,
                  bool stepsOut,
                  double cellSize,
                  int Jm,
                  int metaSAD,
                  int nTrees,
                  double Xext,
                  double Yext,
                  double rmax,
                  double bw1
                  );

   ~CModelSettings();
};

//---------------------------------------------------------------------------
//parameters that are varied in optimization approaches
class CPara
{
public:
	int    metaSR;  //species richness in case of uniform or log-normal metacommunity
	double metaCV;  //cv of abundances for log-normal metacommunity
	double m;
	double r_max;
	double aRec;
	double disp_mean;
	double disp_sd;
	double CNDD_mean;     // factor to calculate heterospecific competition relative to conspecific competition
	double CNDD_sd;
	double pRec_mean;      //mean of species-specific recruitment rates
	double pRec_sd;
	int    trade1_CNDD_pRec;   //trade-off between CNDD and pRec;
	double a_CNDD_pRec;
	double b_CNDD_pRec;
	double c_CNDD_pRec;
	int    trade2_CNDD_abund;  //trade-off between CNDD and abundance in the metacommunity
	double a_CNDD_abund;       //parameters for linear relationship between log(abundance) and CNDD
	double b_CNDD_abund;
	double c_CNDD_abund;
	int    trade3_disp_pRec;   //trade-off between dispersal distance and pRec;
	double a_disp_pRec;
	double b_disp_pRec;
	double c_disp_pRec;

	CPara(){};
	CPara(int    metaSR1,
         double metaCV1,
         double m1,
         double r_max1,
         double aRec1,
         double disp_m1,
         double disp_cv1,
         double CNDD_m1,
         double CNDD_cv1,
         double pRec_m1,
         double pRec_cv1,
         int    trade_CNDD_pRec1,
         double a_CNDD_pRec1,
         double b_CNDD_pRec1,
         double c_CNDD_pRec1,
         int    trade_CNDD_abund1,
         double a_CNDD_abund1,
         double b_CNDD_abund1,
         double c_CNDD_abund1,
         int trade_disp_pRec1,
         double a_disp_pRec1,
         double b_disp_pRec1,
         double c_disp_pRec1
        );
	~CPara(){};
};

//---------------------------------------------------------------------------
class CSpecPara
{
public:
   double meanDisp; //mean dispersal distance
   double probRec;  //recruitment probability without competition

   CSpecPara(){};

   CSpecPara(double mDisp, double pRec) :
      meanDisp{mDisp},
      probRec{pRec}
   {};

   ~CSpecPara(){};
};

//---------------------------------------------------------------------------
#endif

