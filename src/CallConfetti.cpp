//---------------------------------------------------------------------------

#include "Forest.h"
#include <Rcpp.h>
using namespace Rcpp;


//' @useDynLib confettiRbasic
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List EvalConfetti(NumericVector pars,
                  int ngen       = 100,
                  int nsteps_out = 1,
                  int ntrees     = 10000,
                  double xext    = 500,
                  double yext    = 500,
                  double rmax    = 100,
                  double bw1     = 1.0,
                  int metaSAD    = 0
                  )
{
	CModelSettings* pSettings = new CModelSettings(ngen,
                                                  false,
                                                  5,
                                                  ntrees*100,
                                                  metaSAD,
                                                  ntrees,
                                                  xext,
                                                  yext,
                                                  rmax,
                                                  bw1
                                                  );
	CPara* pPara = new CPara(pars[0],
                            pars[1],
                            pars[2],
                            pars[3],
                            pars[4],
                            pars[5],
                            pars[6],
                            pars[7],
                            pars[8],
                            pars[9],
                            pars[10],
                            pars[11],
                            pars[12],
                            pars[13],
                            pars[14],
                            pars[15],
                            pars[16],
                            pars[17],
                            pars[18],
                            pars[19],
                            pars[20],
                            pars[21],
                            pars[22]
                            );

	int seed = as<int>(runif(1,0,999999));

   CForest* pForest = new CForest(seed, pSettings);
   pForest->pPars = pPara;

   //init species and trees
   pForest->initSpecies();
   pForest->initTrees();

   //prepare output data
   int nout{0};
   if (nsteps_out > 1) nout = nsteps_out + 1;
   else                nout = 1;

   IntegerVector Time(nout);

   IntegerVector nSpecies(nout);
   NumericVector Shannon(nout);
   NumericVector Simpson(nout);

   IntegerMatrix Abundance(nout, pForest->SpecMax);
   IntegerMatrix SAD(nout, pForest->MaxSAD);

   NumericMatrix PCF(nout, pForest->nBins1);
   NumericMatrix PropCon(nout, pForest->nBins1);
   NumericMatrix SARq(nout, pForest->SARq_n);

   NumericVector r(pForest->nBins1);
   for (int ir=0; ir < pForest->nBins1; ++ir)
      r[ir] = pForest->rvec1[ir];

   NumericVector Area(pForest->SARq_n);
   for (int i=0; i < pForest->SARq_n; ++i)
      Area[i] = pForest->SARq_scales[i];

   int output_interval = ngen/nsteps_out;

   int nspec{0};
   double shannon{0.0}, simpson{0.0};

   //save initial condition
   if (nsteps_out > 1){

      Time[0] = 0;

      pForest->GetDiversity(nspec, shannon, simpson);

      nSpecies[0] = nspec;
      Shannon[0] = shannon;
      Simpson[0] = 1.0 - simpson;

      for (unsigned int ispec = 0; ispec < pForest->SpecAbund.size(); ++ispec)
         Abundance(0,ispec) = pForest->SpecAbund[ispec];

      pForest->GetSAD();
      for (int iclass = 0; iclass < pForest->MaxSAD; ++iclass)
         SAD(0,iclass) = pForest->SAD[iclass];

      pForest->GetPPA();
      for (int ir = 0; ir < pForest->nBins1; ++ir){
         PCF(0,ir)     = pForest->PCF_all[ir];
         PropCon(0,ir) = pForest->PropCon[ir];
      }

      pForest->GetSARq();
      for (int i = 0; i < pForest->SARq_n; ++i)
         SARq(0,i) = pForest->SARq_m[i];
   }

   int64_t BD_max = ntrees * ngen;

   bool stoprun = false;
   int istep {0};

   while (pForest->BD_total < BD_max){

      for (int i = 0; i < ntrees; ++i) {
         stoprun = pForest->BirthDeathAsync();
         if (stoprun == true){
            Rcout<<"STOP"<<std::endl;
            break;
         }
      }

      if (stoprun == true){
         break;
      }

      ++istep;

      //Rcout<<pForest->BD_total<<std::endl;

      if (istep % output_interval == 0) {

         int iout = istep/output_interval;
         if (nsteps_out == 1) iout = 0;
         Rcout<<"Time step: "<<istep<<std::endl;
         //Rcout<<"iout" <<iout<<std::endl;

         Time[iout] = istep;

         pForest->GetDiversity(nspec, shannon, simpson);

         nSpecies[iout] = nspec;
         Shannon[iout] = shannon;
         Simpson[iout] = 1.0 - simpson;

         for (unsigned int ispec = 0; ispec < pForest->SpecAbund.size(); ++ispec)
            Abundance(iout,ispec) = pForest->SpecAbund[ispec];

         pForest->GetSAD();
         for (int iclass = 0; iclass < pForest->MaxSAD; ++iclass)
            SAD(iout,iclass) = pForest->SAD[iclass];

         pForest->GetPPA();
         for (int ir = 0; ir < pForest->nBins1; ++ir){
            PCF(iout,ir)     = pForest->PCF_all[ir];
            PropCon(iout,ir) = pForest->PropCon[ir];
         }

         pForest->GetSARq();
         for (int i = 0; i < pForest->SARq_n; ++i)
            SARq(iout,i) = pForest->SARq_m[i];
      }
   }  // while BD_total < BD_max

   NumericVector tree_x(ntrees);
   NumericVector tree_y(ntrees);
   IntegerVector tree_specID(ntrees);

   for (int itree = 0; itree < ntrees; ++itree){
      tree_x[itree]      = pForest->TreeList[itree]->X;
      tree_y[itree]      = pForest->TreeList[itree]->Y;
      tree_specID[itree] = pForest->TreeList[itree]->SpecID;
   }

   //Species level output
   IntegerVector speciesID(pForest->SpecMax);
   NumericVector metaRelAbund(pForest->SpecMax);
   NumericVector mDisp(pForest->SpecMax);
   NumericVector CNDD(pForest->SpecMax);
   NumericVector pRec(pForest->SpecMax);

   for (int ispec = 0; ispec < pForest->SpecMax; ++ispec){

      speciesID[ispec] = ispec;

      if (ispec==0) metaRelAbund[ispec] = pForest->CumRelAbundMeta[ispec];
      else          metaRelAbund[ispec] = pForest->CumRelAbundMeta[ispec] -  pForest->CumRelAbundMeta[ispec-1];

      mDisp[ispec] = pForest->SpecPars[ispec].meanDisp;
      CNDD[ispec]  = pForest->InteractMat[ispec][ispec];
      pRec[ispec]  = pForest->SpecPars[ispec].probRec;
   }

   pForest->clearSpecies();
   pForest->clearTrees();

   delete pForest;
   delete pPara;
   delete pSettings;

   DataFrame trees = DataFrame::create(_["X"]         = tree_x,
                                       _["Y"]         = tree_y,
                                       _["SpecID"]    = tree_specID);

   DataFrame species = DataFrame::create(_["speciesID"] = speciesID,
                                         _["metaRelAbund"] = metaRelAbund,
                                         _["meanDisp"] = mDisp,
                                         _["CNDD"] = CNDD,
                                         _["pRec"] = pRec
                                         );

   List output = List::create(_["Generations"] = Time,
                              _["nSpecies"]  = nSpecies,
                              _["Shannon"]   = Shannon,
                              _["Simpson"]   = Simpson,
                              _["Abundance"] = Abundance,
                              _["SAD"]       = SAD,
                              _["Area_m2"]   = Area,
                              _["radius"]    = r,
                              _["SAR"]       = SARq,
                              _["Fr"]        = PropCon,
                              _["PCF"]       = PCF,
                              _["Species"]   = species,
                              _["Trees"]     = trees
                              );
   return(output);
}

