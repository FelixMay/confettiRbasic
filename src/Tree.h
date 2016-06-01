//---------------------------------------------------------------------------

#ifndef TreeH
#define TreeH

#include <cmath>
#include <vector>
#include <map>
#include "Cell.h"

//---------------------------------------------------------------------------
class CTree
{
public:
	unsigned int TreeID;
	double X;
	double Y;
	unsigned int SpecID;

	double R; //ZOI Radius
	double NCI; // neighbourhood crowding index


	CTree(unsigned int id, double x, double y, int spec, double r)
	{
		TreeID = id;
		X = x; Y = y;
		SpecID = spec;
		R = r;
	};

	~CTree()
	{};
};


typedef std::vector<CTree*>::iterator TreeIterV;
typedef std::list<CTree*>::iterator TreeIterL;

//---------------------------------------------------------------------------
#endif
