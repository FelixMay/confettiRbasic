// ---------------------------------------------------------------------------
#include "Forest.h"
#include <time.h>
#include <fstream>
using std::ifstream;
using std::ofstream;

#include <algorithm>
using std::min;

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <map>
using std::map;

#include <string>
using std::string;

//#include <iomanip>

#include <iostream>
using std::cout;
using std::endl;

#include <deque>
using std::deque;

#include <Rcpp.h>
using namespace Rcpp;

// ---------------------------------------------------------------------------
CForest::CForest(int seed, CModelSettings* pset){

	pSettings = pset;

	// Landscape - forest plot
	Xmax = pSettings->m_Xext;
	Ymax = pSettings->m_Yext;

	// Grid
	XCells = Xmax / pSettings->m_cellSize;
	YCells = Ymax / pSettings->m_cellSize;

	Grid = new CCell*[XCells];
	for (int x = 0; x < XCells; ++x)
	   Grid[x] = new CCell[YCells];

 	RandGen1 = new CRandomMersenne(seed);
 	RandGen2 = new StochasticLib1(seed);

	// Point pattern variables all trees
 	nBins1 = ceil(pSettings->m_rmax / pSettings->m_bw1);

	rvec1 = new double[nBins1];
	for (int i = 0; i < nBins1; ++i)
		rvec1[i] = pSettings->m_bw1 * i + 0.5 * pSettings->m_bw1;

	CountAll = new int[nBins1];
	CountCon = new int[nBins1];
	PropCon = new double[nBins1];

	PCF_all = new double[nBins1];
	DenomPCF = new double[nBins1];
	SAR = new double[nBins1];

	int dmin;
	if (Xmax <= Ymax) dmin = floor(Xmax);
	else              dmin = floor(Ymax);

	SARq_scales.clear();
	for (int i = 1; i <= dmin; ++i)
		if (dmin % i == 0)
			SARq_scales.push_back(i*i);

	if (Xmax > Ymax) SARq_scales.push_back(Xmax*Ymax);

	SARq_n = SARq_scales.size();

	SARq_m = new double[SARq_n];

	// Input-Output
	isim = 1;
	irep = 1;
}

// ---------------------------------------------------------------------------
CForest::~CForest() {
	for (TreeIterV itree = TreeList.begin(); itree != TreeList.end(); ++itree)
		delete(*itree);
	TreeList.clear();

	CumRelAbundMeta.clear();
	SpecAbund.clear();
	SARq_scales.clear();

	for (int x = 0; x < XCells; ++x)
		delete[] Grid[x];
	delete[] Grid;

	delete RandGen1;
	delete RandGen2;

   delete[] rvec1;

	delete[] CountAll;
	delete[] CountCon;
	delete[] PropCon;

	delete[] PCF_all;
	delete[] DenomPCF;

	delete[] SARq_m;
}

// ---------------------------------------------------------------------------
void CForest::clearTrees()
{
	for (TreeIterV itree = TreeList.begin(); itree != TreeList.end(); ++itree)
		delete (*itree);
	TreeList.clear();

	SpecAbund.clear();
}

// ---------------------------------------------------------------------------
void CForest::clearSpecies()
{
	for (int ispec = 0; ispec < SpecMax; ++ispec)
		delete[] InteractMat[ispec];
   delete[] InteractMat;

	CumRelAbundMeta.clear();
}

// ---------------------------------------------------------------------------
inline double CForest::Distance(double x1, double y1, double x2, double y2) {
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}


// ---------------------------------------------------------------------------
inline void CForest::PeriodBound(double& xx, double& yy) {
	xx = Xmax * (xx / Xmax - floor(xx / Xmax));
	yy = Ymax * (yy / Ymax - floor(yy / Ymax));
}

// ---------------------------------------------------------------------------
inline void CForest::BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax) {
	xx = xx % Xmax;
	if (xx < 0)
		xx = Xmax + xx;

	yy = yy % Ymax;
	if (yy < 0)
		yy = Ymax + yy;
}

// ---------------------------------------------------------------------------
inline void CForest::BoundGrid(int& xx, int& yy, double& xb, double& yb) {
	if (xx >= XCells) {
		xx = xx - XCells;
		xb = +Xmax;
	}
	else if (xx < 0) {
		xx = xx + XCells;
		xb = -Xmax;
	}

	if (yy >= YCells) {
		yy = yy - YCells;
		yb = +Ymax;
	}
	else if (yy < 0) {
		yy = yy + YCells;
		yb = -Ymax;
	}
}

// ---------------------------------------------------------------------------
inline int CForest::GetRandSpec() {

	int ispec = 0;
	double r1 = RandGen1->Random();
	bool choose = false;

	while (choose == false) {
		if (r1 <= CumRelAbundMeta[ispec])
			choose = true;
		// if (r1 <= CumProbImmi[ispec]) choose = true;
		else
			++ispec;
	}

	return(ispec);
}


// ---------------------------------------------------------------------------
vector<double> CForest::UniformSAD(unsigned int nSpecies)
{
   vector<double>RelAbund(nSpecies, 0);
   vector<double>CumRelAbund(nSpecies, 0);

   for (unsigned int is = 0; is < nSpecies; ++is) {
		RelAbund[is] = 1.0/nSpecies;
		if (is == 0)
			CumRelAbund[is] = RelAbund[is];
		else
			CumRelAbund[is] = CumRelAbund[is - 1] + RelAbund[is];
      //cout<<CumRelAbund[is]<<endl;
	}

	return(CumRelAbund);
}

// ---------------------------------------------------------------------------
vector<double> CForest::LognormSAD(unsigned int nSpecies, unsigned int nIndividuals, double cv_abund)
{
   vector<double> Abund(nSpecies, 0);
   vector<double> RelAbund(nSpecies, 0);
   vector<double> CumRelAbund(nSpecies, 0);

   double meanAbund = static_cast<double>(nIndividuals)/nSpecies;
   double sdAbund = meanAbund*cv_abund;

   double sigma = sqrt(log( (sdAbund*sdAbund)/(meanAbund*meanAbund) + 1.0));
   double mu = log(meanAbund) - 0.5*sigma*sigma;

   int sum = 0;

   for (unsigned int is = 0; is < nSpecies; ++is) {
		Abund[is] = exp(RandGen2->Normal(mu,sigma));
		sum += Abund[is];
   }

   for (unsigned int is = 0; is < nSpecies; ++is) {
		RelAbund[is] = static_cast<double>(Abund[is])/sum;
		if (is == 0)
			CumRelAbund[is] = RelAbund[is];
		else
			CumRelAbund[is] = CumRelAbund[is - 1] + RelAbund[is];
	}

	return(CumRelAbund);
}

// ---------------------------------------------------------------------------
void CForest::initSpecies()
{
	// Init metacommunity if not read from file
	CumRelAbundMeta.clear();

   switch (pSettings->m_metaSAD) {
   case 0:
      CumRelAbundMeta = UniformSAD(pPars->metaSR);
      break;
   case 1:
      CumRelAbundMeta = LognormSAD(pPars->metaSR, pSettings->m_Jm, pPars->metaCV);
      break;
   default:
      CumRelAbundMeta = UniformSAD(pPars->metaSR);
   }

	SpecMax = CumRelAbundMeta.size();

	double disp_spec{0.0}, CNDD_spec{0.0}, pRec_spec{0.0}, meta_relabund{0};

	//Log-normal distribution for mean dispersal distance
	double disp_sigma = sqrt(log(1.0 + (pPars->disp_sd*pPars->disp_sd)/
                                      (pPars->disp_mean*pPars->disp_mean)));
	double disp_mu = log(pPars->disp_mean) - 0.5 * disp_sigma*disp_sigma;


	//Construct species interaction matrix
	InteractMat = new double*[SpecMax];
   for (int ispec = 0; ispec < SpecMax; ++ispec)
      InteractMat[ispec] = new double[SpecMax];

   //Initialize interaction matrix with reference value for all values
   for (int ispec = 0; ispec < SpecMax; ++ispec)
      for (int jspec = 0; jspec < SpecMax; ++jspec)
         InteractMat[ispec][jspec] = 1.0;

	//init species traits
	for (int ispec = 0; ispec < SpecMax; ++ispec){

	   SpecAbund[ispec] = 0;

		//1. Recruitment probability
		pRec_spec = RandGen2->Normal(pPars->pRec_mean, pPars->pRec_sd);
		if (pRec_spec < 0.001)
		   pRec_spec = 0.001;
		if (pRec_spec > 1.0)
		   pRec_spec = 1.0;

		//2. CNDD - conspecific negative density dependence

		//trade-off with recruitment rate
		switch (pPars->trade1_CNDD_pRec) {

   		case 0: //CNDD - normal distribution
   		   CNDD_spec = RandGen2->Normal(pPars->CNDD_mean,pPars->CNDD_sd);
   		   break;
   		case 1: //linear relationship with pRec
   		   CNDD_spec = pPars->a_CNDD_pRec + pPars->b_CNDD_pRec * pRec_spec;
   		   break;
   		case 2:
   		   CNDD_spec = pPars->a_CNDD_pRec + pPars->b_CNDD_pRec * log(pRec_spec);
   		   break;
		   case 3:
		      CNDD_spec = pPars->a_CNDD_pRec + pPars->b_CNDD_pRec * exp(pPars->c_CNDD_pRec * pRec_spec);
		      break;
   		default:
   		   CNDD_spec = RandGen2->Normal(pPars->CNDD_mean,pPars->CNDD_sd);
		}

		if (ispec == 0) meta_relabund = CumRelAbundMeta[0];
		else            meta_relabund = CumRelAbundMeta[ispec] - CumRelAbundMeta[ispec-1];

		//trade-off with metacommunity abundance
		switch (pPars->trade2_CNDD_abund) {
   		case 1:
            CNDD_spec = pPars->a_CNDD_abund + pPars->b_CNDD_abund * meta_relabund;
   		   break;
		   case 2:
		      CNDD_spec = pPars->a_CNDD_abund + pPars->b_CNDD_abund * log(meta_relabund);
		      break;
		   case 3:
		      CNDD_spec = pPars->a_CNDD_abund + pPars->b_CNDD_abund * exp(pPars->c_CNDD_abund * meta_relabund);
		      break;
		   default:
		      ;
		}

		if (CNDD_spec < 1.0) //exclude positive density dependence
		   CNDD_spec = 1.0;
		InteractMat[ispec][ispec] = CNDD_spec;

		//3. Dispersal

		//trade-off with recruitment rate
		switch (pPars->trade3_disp_pRec) {
		case 0: //CNDD - lognormal distribution
		   disp_spec = exp(RandGen2->Normal(disp_mu, disp_sigma));
		   break;
		case 1: //linear relationship with pRec
		   disp_spec = pPars->a_disp_pRec + pPars->b_disp_pRec * pRec_spec;
		   break;
		case 2: //log-linear relationship with pRec
		   disp_spec = pPars->a_disp_pRec + pPars->b_disp_pRec * log(pRec_spec);
		   break;
		case 3:
		   disp_spec = pPars->a_disp_pRec + pPars->b_disp_pRec * exp(pPars->c_disp_pRec * pRec_spec);
		   break;
		default:
		   disp_spec = exp(RandGen2->Normal(disp_mu, disp_sigma));
		}

		if (disp_spec < 0.1)
		   disp_spec = 0.1;

		SpecPars[ispec] = CSpecPara(disp_spec, pRec_spec);
	}

	// Init immigration rate
	if (pPars->m < 0.0)
		m = (2 * Xmax + 2 * Ymax) * pPars->disp_mean / (Pi * Xmax * Ymax);
	else
		m = pPars->m;
}

// ---------------------------------------------------------------------------
void CForest::initTrees()
{
	CTree *pTree1;

	int iX1, iY1;

	BD_total = 0;

	// Init Grid Cells
	for (iX1 = 0; iX1 < XCells; ++iX1)
		for (iY1 = 0; iY1 < YCells; ++iY1)
			Grid[iX1][iY1].InitCell(iX1, iY1,1);

	//grid_steps = (int) ceil(2.0 * Pars->r_max / CellSize);
	grid_steps = (int) ceil(pPars->r_max / pSettings->m_cellSize);

	// random init
	TreeID = 0;
	NTrees = pSettings->m_nTrees;

	double x,y;
	int SpecID;

	for (int i = 0; i < NTrees; ++i) {
		x = RandGen1->Random() * Xmax;
		y = RandGen1->Random() * Ymax;
		SpecID = GetRandSpec();

		pTree1 = new CTree(TreeID, x, y, SpecID, pPars->r_max);
		TreeList.push_back(pTree1);
		++TreeID;

		iX1 = (int) floor(x / pSettings->m_cellSize);
		iY1 = (int) floor(y / pSettings->m_cellSize);
		Grid[iX1][iY1].TreeList.push_back(pTree1);

		++SpecAbund[pTree1->SpecID];
	}
}

// ---------------------------------------------------------------------------
double CForest::GetProbRecruit(double x1, double y1, unsigned int spec_id)
{
	int iX1 = floor(x1/pSettings->m_cellSize);
	int iY1 = floor(y1/pSettings->m_cellSize);

	int iX2, iY2;
	double xbound, ybound, d12;
	CTree* pTree2;

	// calculate NCI
	double NCI{0.0};
	double densNCI{0.0};

	//global neighborhood
	if (pPars->r_max > 99.0){
      for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it)
         NCI = NCI + spec_it->second * InteractMat[spec_id][spec_it->first];

      densNCI = NCI/(Xmax*Ymax);
	}

	//local neighborhood
	else {

      for (int dxi = -grid_steps; dxi <= grid_steps; ++dxi) {
         for (int dyi = -grid_steps; dyi <= grid_steps; ++dyi) {

            iX2 = iX1 + dxi;
            iY2 = iY1 + dyi;

            xbound = 0.0;
            ybound = 0.0;

            // torus boundary
            BoundGrid(iX2, iY2, xbound, ybound);

            // loop over trees in cell 2
            if (Grid[iX2][iY2].TreeList.size() > 0) {
               for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
                  itree2 != Grid[iX2][iY2].TreeList.end(); ++itree2) {

                  pTree2 = (*itree2);
                  d12 = Distance(x1, y1, pTree2->X + xbound,pTree2->Y + ybound);

                  //distance based
                  if (d12 < pPars->r_max) { // do trees overlap?

                        //hyperbolic interaction kernel following Uriarte et al. 2004 JEcol
                        if (d12<0.0001) d12 = 0.0001;
                        NCI = NCI + InteractMat[spec_id][pTree2->SpecID]/d12;

                     } // if overlap
               } // end tree 2
            } // if TreeList > 0

         } // end dx
      } // end dy

      //Standardize NCI by neighborhood area
      if (pPars->r_max > 0.0)
         densNCI = NCI/(pPars->r_max * pPars->r_max * Pi);
	} //if rmax < 99

	// double prob_rec1{0.0};
	// if (NCI == 0)
	//    prob_rec1 = 1.0;
	// else prob_rec1 = (1.0 - densNCI/(pPars->aRec + densNCI));

	//Rcout<<SpecPars[spec_id].probRec<<std::endl;

	double prob_rec1 = (1.0 - densNCI/(pPars->aRec + densNCI)) * SpecPars[spec_id].probRec;

	return(prob_rec1);
}

// ---------------------------------------------------------------------------
// add tree to competition neighborhood and to recruitment grid
void CForest::AddTree(CTree* pTree1)
{
	// add tree to grid
	int iX1 = floor(pTree1->X / pSettings->m_cellSize);
	int iY1 = floor(pTree1->Y / pSettings->m_cellSize);

	Grid[iX1][iY1].TreeList.push_back(pTree1);
}

// ---------------------------------------------------------------------------
// remove tree from competition neighborhood and from recruitment grid
void CForest::RemoveTree(CTree* pTree1) {

	// remove tree from grid
	int iX1 = floor(pTree1->X / pSettings->m_cellSize);
	int iY1 = floor(pTree1->Y / pSettings->m_cellSize);

	pTree1->NCI = 0;

	Grid[iX1][iY1].TreeList.remove(pTree1);
}

// ---------------------------------------------------------------------------
void CForest::GetNewXY(double &x1, double &y1, int idspec) {

	double x0 = x1;
	double y0 = y1;

	if (pPars->disp_mean < 999.0) {

		//log-normal dispersal kernel
		//double r_dist = exp(RandGen2->Normal(SpecPars[idspec].muDisp, SpecPars[idspec].sigmaDisp));

		//Gaussian dispersal kernel following Clark et al. 1999
      double r_dist = sqrt(Pi/2.0) * std::abs(RandGen2->Normal(0.0,pPars->disp_mean)); //Gaussian kernel with the same mean,
      double r_angle = RandGen1->Random() * 2.0 * Pi;

		x1 = x0 + cos(r_angle) * r_dist;
		y1 = y0 + sin(r_angle) * r_dist;

		PeriodBound(x1, y1);
	}

	else {
		// global dispersal
		x1 = RandGen1->Random() * Xmax; // random coordinates
		y1 = RandGen1->Random() * Ymax;
	}
}

// ---------------------------------------------------------------------------
bool CForest::BirthDeathAsync() {

	CTree *pTree = nullptr;
	CTree *pMotherTree = nullptr;

	int ID_Spec;
	double xnew, ynew;

	bool recruit = false;
	bool stop = false;

	double ProbRecruit;

	// choose random tree and kill this tree
	pTree = TreeList[RandGen1->IRandom(0, TreeList.size() - 1)];

	double pkill = 1.0;

	if (RandGen1->Random() < pkill) {

		++BD_total;

		RemoveTree(pTree);

		if (SpecAbund[pTree->SpecID] > 0) --SpecAbund[pTree->SpecID];
		//if (SpecAbund[pTree->SpecID] == 0)  SpecAbund.erase(pTree->SpecID);

		// recruitment and speciation
      recruit = false;
		int ntrials = 0;
		const int max_trials = 999;

		// immigration from metacommunity pool
      if (RandGen1->Random() < m) {

         do {
            ID_Spec = GetRandSpec();

            xnew = RandGen1->Random() * Xmax; // random coordinates
            ynew = RandGen1->Random() * Ymax;

            ProbRecruit = GetProbRecruit(xnew,ynew,ID_Spec);

            ++ntrials;
            if (RandGen1->Random() < ProbRecruit)
               recruit = true;
         } while ((recruit == false) && (ntrials < max_trials));
      }

      // mother tree within the plot - local recruitment
      else {
         do {
            pMotherTree = TreeList[RandGen1->IRandom(0, TreeList.size() - 1)];
            ID_Spec = pMotherTree->SpecID;

            xnew = pMotherTree->X;
            ynew = pMotherTree->Y;
            GetNewXY(xnew, ynew, ID_Spec);

            ProbRecruit = GetProbRecruit(xnew,ynew,ID_Spec);

            ++ntrials;
            if (RandGen1->Random() < ProbRecruit)
               recruit = true;

         } while ((recruit == false) && (ntrials < max_trials));
      }

		if (recruit == true) {
			++TreeID;
			pTree->TreeID = TreeID;
			pTree->X = xnew;
			pTree->Y = ynew;
			pTree->SpecID = ID_Spec;
			++SpecAbund[pTree->SpecID];
			AddTree(pTree);
		}

		else {
			stop = true;
		}
	}

	return(stop);
}

// // ---------------------------------------------------------------------------
// void CForest::OneRun(int isim, int irep)
// {
// 	BD_max = pSettings->m_nGen * TreeList.size();
//
// 	bool stoprun = false;
//
//    int nspec = SpecAbund.size();
//
//    	while (BD_total < BD_max && nspec > 1) {
//
// 		for (int i = 0; i < NTrees; ++i) {
// 			stoprun = BirthDeathAsync();
// 			if (stoprun == true) break;
// 		}
//
//       nspec = SpecAbund.size();
// 		if (stoprun == true) break;
//
// 	}  // while BD_total < BD_max
// }

// ---------------------------------------------------------------------------
void CForest::GetDiversity(int& nspec, double& shannon, double& simpson) {

	double relabund{0.0};

   nspec = 0;
   shannon = 0.0;
   simpson = 0.0;

	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
		if (spec_it->second > 0) {
			++nspec;
			relabund = static_cast<double>(spec_it->second) / NTrees;
			shannon += log(relabund) * relabund;
			simpson += relabund * relabund;
		}
	}

	shannon = -shannon;
}

// ---------------------------------------------------------------------------
int CForest::GetSAD() {

	double Abund;
	int log2Abund;
	int nspec = 0;

	for (int i = 0; i < MaxSAD; ++i)
		SAD[i] = 0;

	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
		Abund = spec_it->second;
		if (Abund > 0) {
			++nspec;
			log2Abund = (int)floor(log(Abund) / log(2.0));
			// log2Abund = (int) ceil(log(Abund)/log(2.0));
			if (log2Abund >= MaxSAD)
				++SAD[MaxSAD - 1];
			else
				++SAD[log2Abund];
		}
	}

	return(nspec);
}

// Point pattern functions ---------------------------------------------------

// ---------------------------------------------------------------------------
// wird benutzt um Anteil Kreisbogen in W und Anteil Kreis in W auszurechnen
inline void CForest::fKK(double D1, double D2, double r, double* ra) {
	double D3, W;
	double S1, S2;
	if (D1 < r) {
		D3 = sqrt(r * r - D1 * D1);
		if (D3 < D2) {
			if (D1 == 0)
				S1 = Pi;
			else
				S1 = atan(D2 / D1);
			if (D1 == 0)
				S2 = Pi;
			else
				S2 = atan(D3 / D1);
			W = S1 - S2;
			ra[0] = ra[0] + 0.5 * (D1 * D3 + r * r * W);
			ra[1] = ra[1] + W * r;
		}
		else {
			ra[0] = ra[0] + 0.5 * D1 * D2;
		}
	}
	else {
		if (D1 == 0)
			W = Pi;
		else
			W = atan(D2 / D1);
		ra[0] = ra[0] + 0.5 * r * r * W;
		ra[1] = ra[1] + W * r;
	}
}

// einmal global definieren: ra:array[0..1] of real

// ------------------------------------------------------------------------------
double CForest::FracRingInWin2(double x, double y, double R) {
	double result;
	double a = Xmax, b = Ymax;
	double ra[2];

	if (R == 0)
		result = 1.0;
	else {
		/*
		if ( Xmax < Ymax )
		{
		a = Xmax;
		b = Ymax;
		}
		else
		{
		a = Ymax;
		b = Xmax;
		}
		 */
		ra[0] = 0;
		ra[1] = 0;
		fKK(a - x, b - y, R, ra);
		fKK(a - x, y, R, ra);
		fKK(y, a - x, R, ra);
		fKK(y, x, R, ra);
		fKK(x, y, R, ra);
		fKK(x, b - y, R, ra);
		fKK(b - y, x, R, ra);
		fKK(b - y, a - x, R, ra);
		result = ra[1] / (2 * Pi * R);
	}
	return result;
}

// ------------------------------------------------------------------------------
double CForest::FracCircleInWin2(double x, double y, double R) {
	double result;
	double a = Xmax, b = Ymax;
	double ra[2];

	if (R == 0)
		result = 1.0;
	else {
		/*
		if ( Xmax < Ymax )
		{
		a = Xmax;
		b = Ymax;
		}
		else
		{
		a = Ymax;
		b = Xmax;
		}
		 */
		ra[0] = 0;
		ra[1] = 0;
		fKK(a - x, b - y, R, ra);
		fKK(a - x, y, R, ra);
		fKK(y, a - x, R, ra);
		fKK(y, x, R, ra);
		fKK(x, y, R, ra);
		fKK(x, b - y, R, ra);
		fKK(b - y, x, R, ra);
		fKK(b - y, a - x, R, ra);
		result = ra[0] / (Pi * R * R);
	}
	return result;
}

// ---------------------------------------------------------------------------
// calculates the fraction of the circle with radius r and center x,y in Window
inline double CForest::FracRingInWin(double x, double y, double R) {
	// see Diggle 2003, pg. 51 and Goreaud & Pelissier 1999. J. Veg. Sci.
	double d1 = min(x, Xmax - x);
	double d2 = min(y, Ymax - y);
	double d3 = sqrt(d1 * d1 + d2 * d2);

	double alpha_out;

	if (R < d3) // corner outside the circle
		alpha_out = 2 * acos(min(d1, R) / R) + 2 * acos(min(d2, R) / R);
	else // corner inside
		alpha_out = Pi * 0.5 + acos(d1 / R) + acos(d2 / R);

	return 1.0 - (alpha_out / (2 * Pi));
}

// ---------------------------------------------------------------------------
void CForest::GetPPA() {

	CTree* pTree1;
	CTree* pTree2;


	double rmax = pSettings->m_rmax;
	double bw1  = pSettings->m_bw1;


	// ----------------------
	int Rgrid = (int) ceil(rmax / pSettings->m_cellSize);

	map<int, int>::iterator spec_it1;
	map<int, int>::iterator spec_it2;


	for (int ibin1 = 0; ibin1 < nBins1; ++ibin1) {
		CountAll[ibin1] = 0;
		CountCon[ibin1] = 0;
		PropCon[ibin1] = 0.0;
		DenomPCF[ibin1] = 0.0;
		PCF_all[ibin1] = 0.0;
	}

	int ibin1;
	int iX2, iY2;
	int Xcells2 = Xmax / pSettings->m_cellSize;
	int Ycells2 = Ymax / pSettings->m_cellSize;

	double xydist, xbound, ybound;


	//main nested loop from tree i to tree 2
	for (int iX1 = 0; iX1 < Xcells2; ++iX1) {
		for (int iY1 = 0; iY1 < Ycells2; ++iY1) {

			// loop over trees in cell 1
			for (TreeIterL itree1 = Grid[iX1][iY1].TreeList.begin();
				itree1 != Grid[iX1][iY1].TreeList.end(); ++itree1) {

				pTree1 = *itree1;

				// loop over neighboring cells up to maximum distance
				for (int diX = -Rgrid; diX <= Rgrid; ++diX) {
					for (int diY = -Rgrid; diY <= Rgrid; ++diY) {

						iX2 = iX1 + diX;
						iY2 = iY1 + diY;

						xbound = 0;
						ybound = 0;

						// no edge correcion
						if ((iX2 < 0 || iY2 < 0) || (iX2 >= XCells) || (iY2 >= YCells))
							continue;

						// loop over trees in cell 2
						for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
							itree2 != Grid[iX2][iY2].TreeList.end(); ++itree2) {

							pTree2 = *itree2;

							if (pTree1 != pTree2) {
								xydist = Distance(pTree1->X, pTree1->Y,
									pTree2->X + xbound, pTree2->Y + ybound);
								if (xydist <= rmax) {
									ibin1 = floor(xydist / bw1);
									++CountAll[ibin1];

									if (pTree1->SpecID == pTree2->SpecID)
										++CountCon[ibin1];

								} // if (d < rmax)

							} // if tree1 != tree2

						} // end tree 2

					} // end dy
				} // end dx

			} // end for tree 1

		} // end for y1
	} // end for x1
	// end main loop

	// calculate Wiegand-Moloney edge correction
	for (int ibin1 = 0; ibin1 < nBins1; ++ibin1) {
		for (unsigned int itree = 0; itree < TreeList.size(); ++itree) {
			pTree1 = TreeList[itree];
				DenomPCF[ibin1] = DenomPCF[ibin1] +
                              (2.0 * Pi * rvec1[ibin1] * bw1) * // area of annuli
                              FracRingInWin(pTree1->X, pTree1->Y, rvec1[ibin1]);
		}

		PCF_all[ibin1] = (Xmax * Ymax / (NTrees - 1))* CountAll[ibin1] / DenomPCF[ibin1];

		if (CountAll[ibin1] > 0)
			PropCon[ibin1] = static_cast<double>(CountCon[ibin1]) / CountAll[ibin1];
	}
}

// ------------------------------------------------------------------------------
void CForest::GetSRLocal(double sq_size, double& m_SR) {

	int Xsq = floor(Xmax / sq_size);
	int Ysq = floor(Ymax / sq_size);

	// 3D array
	CSpecSquare** SRgrid;

	SRgrid = new CSpecSquare*[Xsq];
	for (int iX = 0; iX < Xsq; ++iX)
		SRgrid[iX] = new CSpecSquare[Ysq];

	for (int iX = 0; iX < Xsq; ++iX)
		for (int iY = 0; iY < Ysq; ++iY)
			SRgrid[iX][iY].InitspecSquare(iX, iY);

	CTree* pTree1;

	int iX, iY;

	for (unsigned int itree = 0; itree < TreeList.size(); ++itree) {
		pTree1 = TreeList[itree];
		iX = (int)floor(pTree1->X / sq_size);
		iY = (int)floor(pTree1->Y / sq_size);
		++SRgrid[iX][iY].Spec[pTree1->SpecID];
	}

	// calculate mean and sd
	int sum_SR = 0;
	for (int iX = 0; iX < Xsq; ++iX)
		for (int iY = 0; iY < Ysq; ++iY)
			sum_SR = sum_SR + SRgrid[iX][iY].Spec.size();

	m_SR = static_cast<double>(sum_SR) / (Xsq * Ysq);

	// free storage
	for (int iX = 0; iX < Xsq; ++iX)
		delete[] SRgrid[iX];
	delete[] SRgrid;

}

// ------------------------------------------------------------------------------
void CForest::GetSARq() {

	for (int i = 0; i < (SARq_n - 1); ++i) {
		GetSRLocal(sqrt(SARq_scales[i]), SARq_m[i]);
	}

	int nspec = 0;

	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it)
		if (spec_it->second > 0)
			++nspec;

	SARq_m[SARq_n - 1] = nspec;
}

// ------------------------------------------------------------------------------
