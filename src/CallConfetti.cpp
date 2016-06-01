//---------------------------------------------------------------------------

#include "Forest.h"
#include <Rcpp.h>
using namespace Rcpp;


//' @useDynLib confettiRbasic
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List EvalConfetti(NumericVector pars,
                  int ngen = 100,
                  int ntrees = 10000,
                  double xext = 500,
                  double yext = 500,
                  double rmax = 100,
                  double bw1 = 1.0
                 )


{
	CModelSettings* pSettings = new CModelSettings(ngen,
                                                  false,
                                                  5,
                                                  ntrees*100,
                                                  ntrees,
                                                  xext,
                                                  yext,
                                                  rmax,
                                                  bw1
                                                  );

	CPara* pPara = new CPara(pars[0], pars[1], pars[2],
                            pars[3], pars[4], pars[5],
                            pars[6], pars[7]
                            );

	int seed = as<int>(runif(1,0,999999));

   CForest* pForest = new CForest(seed, pSettings);
   pForest->pPars = pPara;

   pForest->initSpecies();
   pForest->initTrees();

   pForest->OneRun(1,1);

   int nspec{0};
   double shannon{0.0}, simpson{0.0};

   pForest->GetDiversity(nspec, shannon, simpson);

 	int64_t BDtotal = pForest->BD_total;

   IntegerVector SAD(12);

   pForest->GetSAD();
   for (int iclass = 0; iclass < SAD.size(); ++iclass)
      SAD[iclass] = pForest->SAD[iclass];

   pForest->GetPPA();
   pForest->GetSARq();

   NumericVector x(ntrees);
   NumericVector y(ntrees);
   IntegerVector specID(ntrees);
   IntegerVector abund(nspec);

   for (int itree = 0; itree < ntrees; ++itree){
      x[itree]      = pForest->TreeList[itree]->X;
      y[itree]      = pForest->TreeList[itree]->Y;
      specID[itree] = pForest->TreeList[itree]->SpecID;
   }

   for (int ispec = 0; ispec < nspec; ++ispec)
      abund[ispec] = pForest->SpecAbund[ispec];


   NumericVector PCF(pForest->nBins1);
   NumericVector PropCon(pForest->nBins1);

   NumericVector SARq(pForest->SARq_n);
   NumericVector Area(pForest->SARq_n);

   for (int ir=0; ir<pForest->nBins1; ++ir){
      PCF[ir] = pForest->PCF_all[ir];
      PropCon[ir] = pForest->PropCon[ir];
   }

   for (int i=0; i < pForest->SARq_n; ++i){
		SARq[i] = pForest->SARq_m[i];
      Area[i] = pForest->SARq_scales[i];
   }

   pForest->clearSpecies();
   pForest->clearTrees();

   delete pForest;
   delete pPara;
   delete pSettings;

   List output = List::create(_["seed"]     = seed,
                              _["X"]        = x,
                              _["Y"]        = y,
                              _["SpecID"]   = specID,
                              _["Abund"]    = abund,
                              _["nSpecies"] = nspec,
                              _["Shannon"]  = shannon,
                              _["Simpson"]  = 1 - simpson,
                              _["SAD"]      = SAD,
                              _["Area_m2"]  = Area,
                              _["SAR"]      = SARq,
                              _["Fr"]       = PropCon,
                              _["PCF"]      = PCF);

   return(output);
}

