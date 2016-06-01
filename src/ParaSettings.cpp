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
             double disp1,
             double m_CNDD1,
             double cv_CNDD1
             ) :
   metaSR{metaSR1},
	metaCV{metaCV1},
   m{m1},
   r_max{r_max1},
   aRec{aRec1},
   disp_mean{disp1},
   m_CNDDspec{m_CNDD1},
   cv_CNDDspec{cv_CNDD1},
   sd_CNDDspec{m_CNDDspec * cv_CNDDspec}
{}
// ---------------------------------------------------------------------------

