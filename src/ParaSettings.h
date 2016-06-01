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
	double m_CNDDspec;     // factor to calculate heterospecific competition relative to conspecific competition
	double cv_CNDDspec;
	double sd_CNDDspec;

	CPara(){};
	CPara(int    metaSR1,
         double metaCV1,
         double m1,
         double r_max1,
         double aRec1,
         double disp1,
         double m_CNDD,
         double cv_CNDD1
        );
	~CPara(){};
};

//---------------------------------------------------------------------------
#endif

