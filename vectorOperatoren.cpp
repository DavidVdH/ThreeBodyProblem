#include <vector>
#include "vec.h"

using namespace std;

/*
These functions allow to do elementwise operations in STL vectors. The * operator was overloaded twice in order to be able to multiply a vector by just one scalar or each element with a different scalar. Testing is done in main(). 
*/


inline vector<vec> operator +(vector<vec> vecLHS, vector<vec> vecRHS)
{
	for(int i = 0; i < vecLHS.size(); i++)
	{
		vecLHS[i] += vecRHS[i]; 
	}
	
	return vecLHS;
}


inline vector<vec> operator *(vector<double> scalars, vector<vec> vecRHS)
{
	for(int i = 0; i < vecRHS.size(); i++)
	{
		vecRHS[i] *= scalars[i]; 
	}

	return vecRHS;
}

inline vector<vec> operator *(double scalar, vector<vec> vecRHS)
{
	for(int i = 0; i < vecRHS.size(); i++)
	{
		vecRHS[i] *= scalar; 
	}

	return vecRHS;
}


inline vec sum(vector<vec> Vec)
{
	vec vecSum;
	for(int i = 0; i < Vec.size(); i++)
	{
		vecSum += Vec[i]; 
	}

	return vecSum;
}


/*

MAIN IS ENKEL OM TE TESTEN OF ALLES WERKT

int main()
{
	vector<vec> container1, container2;
	vector<double> mass;
	vec vec1(2,3);

	for(int i = 0; i<3; i++)
	{
		container1.push_back(vec(i,i+1)); //[ (0,1) ; (1,2) ; (2;3) ]	
		container2.push_back(vec(i,i));   //[ (0,0) ; (1,1) ; (2,2) ]
		mass.push_back( i*i );		  //[ 0 ; 1 ; 4 ]
	}

	vector<vec> scal = 5 * container1;
	
	for(int i = 0; i<scal.size(); i++)
	{scal[i].print();}

	vector<vec> som = container1 + container2;

	for(int i = 0; i<som.size(); i++)
	{som[i].print();}
	
	vector<vec> force = mass * container1;

	for(int i = 0; i < force.size(); i++)
	{force[i].print();}

	sum(container1).print();

	vector<vec> relDist = container1 - vec1;

	for(int i = 0; i < relDist.size(); i++)
	{relDist[i].print();}

	return 0;
}


*/
