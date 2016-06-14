//---------------------------------------------------------------------------

#ifndef ForestH
#define ForestH

#include "Tree.h"
#include "ParaSettings.h"
#include "randomc.h"
#include "stocc.h"
#include <fstream>
#include <sstream>
#include <Rcpp.h>

//---------------------------------------------------------------------------
class CForest
{
public:
	CPara *pPars;
	CModelSettings *pSettings;

	//Immigration rate
	double m;

	//Landscape
	double Xmax;
	double Ymax;

	//Trees
	int NTrees;
	std::vector<CTree*> TreeList;
	int64_t TreeID;

	//Species
	int SpecMax;
	std::vector<double> CumRelAbundMeta;

	std::map<int,int> SpecAbund;       // map with first --> key, second --> abund
	std::map<int,CSpecPara> SpecPars;  // species specific parameters

	double **InteractMat;      //matrix with species-specific interaction coefficients

	//Runtime
	//int64_t BD_max;

	//Grid
	int XCells;
	int YCells;
	int grid_steps;

	CCell **Grid;

	//random number generators
	CRandomMersenne *RandGen1;
	StochasticLib1 *RandGen2;

	//Point pattern variables
	int nBins1;
	double *rvec1;

	int *CountAll;
	int *CountCon;
	double *PropCon;
	double *PCF_all;
	double *DenomPCF;
	double *SAR;

	// Point pattern variables single species
	std::map<int,int>::iterator spec_it;

	// Variables for quadrat-based SAR
	std::vector<double> SARq_scales; // square side length
	int SARq_n;      //number of scales;
	double *SARq_m;  // mean species richness

	//index variables
	int isim;
	int irep;

	//Output Variables
	int64_t BD_total;

	static const int MaxSAD = 12;
	int SAD[MaxSAD];      //Species abundance distributions as octave curve 2^0 - 2^11

	//Private functions
	inline int GetRandSpec(); //draw a species from the species pool

	inline double Distance(double x1, double y1, double x2, double y2);

	inline void PeriodBound(double& xx, double& yy);
	inline void BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax);
	inline void BoundGrid(int& xx, int& yy, double& xb, double& yb);

	inline double FracRingInWin(double x, double y, double R);
	inline void fKK( double D1, double D2, double r, double* ra);

	double FracRingInWin2( double x, double y, double R);
	double FracCircleInWin2( double x, double y, double R);

	void GetSRLocal(double sq_size, double& m_SR);

	double GetProbRecruit(double x0, double y0, unsigned int spec_id);

	void AddTree(CTree* pTree1);
	void RemoveTree(CTree* pTree1);

	//int GetSpecMaxAbund();
	void GetNewXY(double &x1, double &y1, int idspec);
	std::vector<double> UniformSAD(unsigned int nSpecies);
	std::vector<double> LognormSAD(unsigned int nSpecies, unsigned int nIndividuals, double cv_abund);

//public:
	CForest(int seed, CModelSettings* pset);
	~CForest();

	void initTrees();
	void initSpecies();
	bool BirthDeathAsync();
	void UpdateTrees();
	void clearTrees();
	void clearSpecies();
	void GetPPA();
	double GetShannon();
	void GetDiversity(int& nspec, double& shannon, double& simpson);
	int GetSAD();
	void GetSARq();
	//void OneRun(int isim, int irep);

	std::string IntToString(int i)
	{
	  std::ostringstream os;
	  os<<i;
	  return os.str();
	}

};
//---------------------------------------------------------------------------
#endif
