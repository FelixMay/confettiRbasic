// ---------------------------------------------------------------------------
#include "ParaSettings.h"
#include <fstream>
#include <iostream>
#include <string>
using std::string;

// ---------------------------------------------------------------------------
CModelSettings::CModelSettings() :
   m_nGen{100},
   m_stepsOut{0},
   m_cellSize{5.0},
   m_Jm{1000000},
   m_metaSAD{0},
   m_nTrees{10000},
   m_Xext{500.0},
   m_Yext{500.0},
   m_rmax{100.0},
   m_bw1{1.0}
{}

CModelSettings::CModelSettings(
   int nGen,
   bool stepsOut,
   double cellSize,
   int Jm,
   int metaSAD,
   int nTrees,
   double Xext,
   double Yext,
   double rmax,
   double bw1
   ) :
   m_nGen{nGen},
   m_stepsOut{stepsOut},
   m_cellSize{cellSize},
   m_Jm{Jm},
   m_metaSAD{metaSAD},
   m_nTrees{nTrees},
   m_Xext{Xext},
   m_Yext{Yext},
   m_rmax{rmax},
   m_bw1{bw1}
{}

// ---------------------------------------------------------------------------
CModelSettings::~CModelSettings(){};

// ---------------------------------------------------------------------------
CPara::CPara(int    metaSR1,
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
             int trade_disp_pRec1,
             double a_disp_pRec1,
             double b_disp_pRec1,
             double c_disp_pRec1
             ) :
   metaSR{metaSR1},
	metaCV{metaCV1},
   m{m1},
   r_max{r_max1},
   aRec{aRec1},
   disp_mean{disp_m1},
   disp_sd{disp_m1* disp_cv1},
   CNDD_mean{CNDD_m1},
   CNDD_sd{CNDD_m1 * CNDD_cv1},
   pRec_mean{pRec_m1},
   pRec_sd{pRec_m1 * pRec_cv1},
   trade1_CNDD_pRec{trade_CNDD_pRec1},
   a_CNDD_pRec{a_CNDD_pRec1},
   b_CNDD_pRec{b_CNDD_pRec1},
   c_CNDD_pRec{c_CNDD_pRec1},
   trade2_CNDD_abund{trade_CNDD_abund1},
   a_CNDD_abund{a_CNDD_abund1},
   b_CNDD_abund{b_CNDD_abund1},
   trade3_disp_pRec{trade_disp_pRec1},
   a_disp_pRec{a_disp_pRec1},
   b_disp_pRec{b_disp_pRec1},
   c_disp_pRec{c_disp_pRec1}
{}
// ---------------------------------------------------------------------------

