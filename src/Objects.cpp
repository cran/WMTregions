/*
 * Objects.cpp
 *
 *  Created on: Jun 9, 2010
 *      Author: chetnik
 */

// Here we describe different objects which are used in the algorithm
// P.Bazovkin, 2009.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <string.h>
#include <bitset>
//#include "R/R.h"
//#include "R/Rmath.h"

#define NDEBUG
#define BOOST_NO_MEMBER_TEMPLATE_FRIENDS

#include <boost/dynamic_bitset.hpp>
#include <boost/cstdint.hpp>
//#include <boost/filesystem.hpp>

#ifdef INCL_OBJECTS

using namespace std;

// ******************************* Timer routine **************************************************

   clock_t starttimer, finishtimer;
   list<double>    timestamps;
   vector<double>  summtimes;
   double cumulcount;

   void StartTiming()
   {
	   timestamps.clear();
	   starttimer = clock();
	   cumulcount = 0.0f;
   }

   void RecordTime()
   {
		finishtimer = clock();
		double durationtimer = (double)(finishtimer - starttimer) / CLOCKS_PER_SEC;

		timestamps.push_front(cumulcount);
		timestamps.push_front(durationtimer);
	}

   clock_t resumecount, stopcount;

   void ResumeCumulTime()
   {
	   resumecount = clock();
   }

   void StopCumulTime()
   {
	   stopcount   = clock();
	   cumulcount += (double)(stopcount - resumecount) / CLOCKS_PER_SEC;
   }

class SpaceObject
{
protected:

	int dim;           // Dimension

public:

	SpaceObject()
	{
		dim = 3;       // Default setting to 3d
	}

	SpaceObject(int _dim)
	{
		dim = _dim;
	}

	int Dim()
	{
		return dim;
	}
};

class Point : public SpaceObject
{

public:

	// Coordinates
	vector<float> coord;

public:

	Point()
	{
		vector<float> invec(dim, 0);
		coord = invec;
	}

	Point(int _dim):
				SpaceObject(_dim)
	{
		vector<float> invec(dim, 0);
		coord = invec;
	}

	// For 3d-case for convenience
	//property float x
	//{
	//	float get()
	//	{
	//		return coord[0];
	//	}

	//	void set(float value)
	//	{
	//		coord[0] = value;
	//	}
	//}

	//property float y
	//{
	//	float get()
	//	{
	//		return coord[1];
	//	}

	//	void set(float value)
	//	{
	//		coord[1] = value;
	//	}
	//}

	//property float z
	//{
	//	float get()
	//	{
	//		return coord[2];
	//	}

	//	void set(float value)
	//	{
	//		coord[2] = value;
	//	}
	//}
};

class Vector : public SpaceObject
{
public:

	// Coordinates
	vector<float> coord;

	// "Combinatorial" coordinates - numbers of starting and ending point in a permutation
	int           comb[2];

	Vector()
	{
		for(int i=0; i<dim; i++)
		{
			//float refpoint = 0;
			float refpoint = 1;
			coord.push_back(refpoint);
		}

		coord[0] = 1;

	}

	Vector(int _dim):
				SpaceObject(_dim)
	{
		for(int i=0; i<dim; i++)
		{
			//float refpoint = 0;
			float refpoint = 1;
			coord.push_back(refpoint);
		}

		coord[0] = 1;

	}

	// Create a vector to the point
	void ToPoint(Point _point)
	{
		if (this->Dim() == _point.Dim())
		{
			this->coord.assign(_point.coord.begin(), _point.coord.end());
		}
		else if (this->Dim() > _point.Dim())
		{
			//vector<float> padding;
			//padding.assign( this->Dim() - _point.Dim(), 0 );

			this->coord.assign(_point.coord.begin(), _point.coord.end());
			this->coord.insert(this->coord.end()+1, this->Dim() - _point.Dim(), 0);      // Padding set to 0
		}
		else
		{
			// Error
		}
	}

	static Vector CreateVector(Point point)
	{
		Vector vect(point.Dim());
		vect.ToPoint(point);

		return vect;
	}

	static float ScalarMultiply(Vector v1, Vector v2)
	{
		float res = 0;
		for(int i=0; i<v1.Dim() && i<v2.Dim(); i++)
		{
			res = res + v1.coord[i] * v2.coord[i];
		}

		return res;

	}

	bool PointOrder( Point point1, Point point2)
	{
		return  ScalarMultiply(*this, CreateVector(point1)) > ScalarMultiply(*this, CreateVector(point2));
	}

	// For default vector
	static bool InitialPointOrder( Point point1, Point point2)
	{
		Vector vect(point1.Dim());
		return  ScalarMultiply(vect, CreateVector(point1)) > ScalarMultiply(vect, CreateVector(point2));
	}

	// Determining the current vector through end and start points
	void FromPointToPoint(Point _start, Point _end)
	{
		for(int i = 0; i < dim; i++)
		{
			this->coord[i] = _end.coord[i] - _start.coord[i];
		}
	}

	// Reverse the vector's direction
	void Reverse()
	{
		for(int i = 0; i < dim; i++)
		{
			this->coord[i] = (-1) * this->coord[i];
		}
	}

	float GetLength()
	{
		return sqrt(Vector::ScalarMultiply(*this, *this));
	}

	// Returns cosine of the angle between vectors
	static float AngleCos(Vector vect1, Vector vect2)
	{
		return Vector::ScalarMultiply(vect1, vect2) / (vect1.GetLength() * vect2.GetLength());
	}

	// Return a 3d-vector, orthogonal to vect1 and vect2 and "codirected" with dir
	static Vector FindOrthogonal3d(Vector vect1, Vector vect2, Vector _dir)
	{
		Vector vec(3);

		vec.coord[0] = 1;
		vec.coord[2] = (vect2.coord[1]*vect1.coord[0] - vect2.coord[0]*vect1.coord[1])/(vect2.coord[2]*vect1.coord[1] - vect2.coord[1]*vect1.coord[2]);
		vec.coord[1] = (-vect1.coord[0] - vect1.coord[2]*vec.coord[2]) / vect1.coord[1];

		// Specifying direction
		if(Vector::ScalarMultiply(vec, _dir) < 0) vec.Reverse();

		return vec;
	}

	// Returns the sum of elements of the vector
	float SumOfElements()
	{
		float res = 0;
		for(int j = 0; j < this->dim; j++)
			res += this->coord[j];

		return res;
	}

	// Vector addition
	void Add(Vector _vec)
	{
		for(int i = 0; i < this->dim; i++)
		{
			this->coord[i] += _vec.coord[i];
		}
	}

	// Multiplying by a constant
	void Scale(float _lambda)
	{
		for(int i = 0; i < this->dim; i++)
		{
			this->coord[i] = this->coord[i] * _lambda;
		}
	}

	// Setting length to 1
	void Normalize()
	{
		float divratio = 1.0f / this->GetLength();

		this->Scale(divratio);
	}

};

//#ifdef incl_objects
Vector stat_direction;               // !!!Warning!!! (Static member)
//#endif

class Permutation
{
public:

	vector<Point> points;

	///*static*/ Vector stat_direction;               // !!!Warning!!! (Static member)

	Permutation()
	{
	}

	// For statically defined vector - don't use in general tasks!
	static bool StaticPointOrder( Point point1, Point point2)
	{
		return  Vector::ScalarMultiply(stat_direction, Vector::CreateVector(point1)) < Vector::ScalarMultiply(stat_direction, Vector::CreateVector(point2));
		//return true;
	}

	// Generalized support function based reordering
	void Support(Vector _p)
	{
		//bool (Vector::*pointOrder)(Point , Point ) = &Vector::PointOrder;
		//sort(points.begin(), points.end(), (_p.*pointOrder));

		stat_direction = _p;      // Passing argument to static vector
		sort(points.begin(), points.end(), &Permutation::StaticPointOrder);

	}


	// Support function based reordering only for initial vectors (defaults)
	void InitialSupport()
	{
		sort(points.begin(), points.end(), &Vector::InitialPointOrder);
	}
};

class ExtremePoint : public Point
{
public:

	vector<int> index_perm;                        // Index permutation

	boost::dynamic_bitset<unsigned int> neighs;                // Neighbours of the current extreme point (basing on the index_perm)

	int _border;                                   // Index of a border element

	ExtremePoint(int _dim):
					Point(_dim)
	{
	}

	ExtremePoint(int _dim, int __border):
					Point(_dim)
	{
		_border = __border;
	}

	// Returns the index of a border element in the initial permutation
	int GetInitBorder()
	{
		return index_perm[_border];
	}

	// Method to return ordered map of n points (basing of index_perm), where their position is displayed (1 in array - before border point; 0 - after)
	vector<char> GetMap()
	{
		vector<char> elem_map(this->index_perm.size());

		for(int i = 0; i < this->index_perm.size(); i++)
		{
			// For all points before the border
			if(i < _border)
			{
				elem_map[index_perm[i]] = 1;
			}
			// For the other points
			else
			{
				elem_map[index_perm[i]] = 0;
			}

		}

		// Mark the border point in another way
		elem_map[index_perm[_border]] = -1;

		return elem_map;
	}

	// Method to return ordered map of n points (basing of index_perm), where their position is displayed (1 in array - before border point; 0 - after)
	boost::dynamic_bitset<unsigned int> GetSmallPartialMap()
	{
		int amount = this->index_perm.size();

		boost::dynamic_bitset<unsigned int> elem_map(amount/* + 32*/);

		for(int i = 0; i < amount; i++)
		{
			// For all points before the border
			if(i < _border)
			{
				elem_map[index_perm[i]] = 1;
			}
			// For the other points
			else
			{
				elem_map[index_perm[i]] = 0;
			}

		}

		//// Mark the border point in another way
		//int log_length = log(amount) / log(2) + 1;
		//bitset<log_length> border_bits(_border);

		//for(int i = 0; i < log_length; i++)
		//{
		//	elem_map[amount+i] = border_bits[i];
		//}

		return elem_map;
	}

	// The partial globally but full locally order of the points in a facet
	static bool PartialEPointOrder(ExtremePoint ep1, ExtremePoint ep2)
	{
		return ep1.index_perm[ep1._border] > ep2.index_perm[ep2._border];
		//return ep1._border > ep2._border;
	}

	bool Equal(ExtremePoint ep)
	{
		return this->index_perm[this->_border] == ep.index_perm[ep._border];
	}

};

class DataVector : public Vector
{
public:

	int start;                   // Index of the starting data point
	int end;                     // Index of the ending data point

	bool start_typeI;                  // Associating cluster type
	bool end_typeI;                  // Associating cluster type

	DataVector()
	{
		start = -3;
		end   = -3;
	}

	DataVector(int _dim, int _start, int _end):
											Vector(_dim)
	{
		start = _start;
		end   = _end;
	}

	// If generating data points are defined
	bool defined()
	{
		return this->start >= 0 && this->end >= 0;
	}

};


class kFace : public SpaceObject
{
protected:
	
	int real_dim;     // Real dimension in d-space

public:

	Vector normalvec;                      // Some normal to the k-face
	
	// Combinatoric stuff
	list<int>   anchors;                // Array of anchors for defining subsets
	list<int>   cardinals;              // Cardinals of defining subsets in correspondence with anchors

	vector<int> index_perm;             // $\pi$ function for the facet (more precisely, it should be $\Pi$ class of $\pi$-s )


	kFace(int _dim, int _K):
			SpaceObject(_dim),
			normalvec(_dim)
	{
		real_dim = _K;
	}

	int K()
	{
		return real_dim;
	}

	// Increasing the real dimension by 1
	void CatchOneDoF()
	{
		if(this->K() < this->Dim()-1) this->real_dim++;
	}

	// Reducing the real dimension by 1
	void RelaxOneDoF()
	{
		if(this->K() > 0) this->real_dim--;
	}

	vector<float> SolveLinearSystem(vector<vector<float> > A, vector<float> b)
	{
		// The first stage of Gaussian elimination
		for(int i = 0; i < dim; i++)
		{
			// Interchanging rows if the main element is 0
			//if(A[i][i] == 0)
			if(fabs(A[i][i]) < 0.001f)
			{
				for(int j = i+1; j < dim; j++)
				{
					if(A[j][i] != 0)
					{
						vector<float> tempvc = A[i];
						A[i] = A[j];
						A[j] = tempvc;
						float tempb = b[i];
						b[i] = b[j];
						b[j] = tempb;
						break;
					}
				}
			}

			// Normalizing the current row
			float coeff = A[i][i];
			for(int k = i; k < dim; k++)
			{
				A[i][k] = A[i][k] / coeff; 
			}
			b[i] = b[i] / coeff;

			// Setting to 0 the lower subcolumn
			for(int j = i+1; j < dim; j++)
			{
				coeff = A[j][i];
				
				for(int k = i; k < dim; k++)
				{
					A[j][k] = A[j][k] - A[i][k] * coeff; 
				}
				b[j] = b[j] - b[i] * coeff; 
			}
		}

		// The second stage of Gaussian elimination
		for(int i = dim-1; i >= 0; i--)
		{
			// Setting to 0 the upper subcolumn
			for(int j = i-1; j >= 0; j--)
			{
				float coeff = A[j][i];
				
				A[j][i] = 0; 
				b[j] = b[j] - b[i] * coeff; 
			}
		}

		return b;
	}
	
};


class Facet : public kFace
{
public:

	list<ExtremePoint> nodes;  

	//Vector direction;					   // Orientation of the facet

	vector<int> basic_perm;				   // The index array of the point adding which we have got this facet

	int zone;                              // Indicates zone in permutation, in which facet's nodes have differences: 0 - before border_index, 1 - after it. 

	bool truncated;                        // If this facet is truncated

	list<ExtremePoint> extra_nodes;        // For truncated one

	bool marked;                           // Whether all ridges are already marked

	Facet(int _dim):
				kFace(_dim,_dim-1)
	{
		truncated = false;
		marked    = false;

	}

	Facet(int _dim, bool _trunc):
				kFace(_dim,_dim-1)
	{
		truncated = _trunc;
		marked    = false;

	}

	Facet(kFace _face):
		kFace(_face)
	{
		marked    = false;
	}


	// Gives the hash code for the facet IF nodes are ordered (!!!Warning!!!)
	// Returns (n + d * 32)-bit hash (n for the basic map and 32 for each epoint) 
	boost::dynamic_bitset<unsigned int> GetHashMap()
	{
		int amount = nodes.begin()->index_perm.size();

		int volume = amount + dim * 32;
		
		boost::dynamic_bitset<unsigned int> hash(volume);

		// Implementing basic map into the hash
		boost::dynamic_bitset<unsigned int> basic_map = nodes.begin()->GetSmallPartialMap();
		for(int j = 0; j < amount; j++)
		{
			hash[j] = basic_map[j];
		}

		int counter = amount;

		// Implementing individual treats of the epoints into the hash
		list<ExtremePoint>::iterator nodeit;                   // Iterator for the nodes of facets
		for( nodeit = nodes.begin() ; nodeit != nodes.end() ; nodeit++)
		{
			bitset<32> border_bits(nodeit->GetInitBorder());
			for(int i = 0; i < 32; i++)
			{
				hash[counter+i] = border_bits[i];
			}

			counter += 32;

			// !Can be not only for truncated (now is not max efficient)
			if(this->truncated)
			{
				hash[nodeit->index_perm[nodeit->_border]] = false;
			}
		}

		return hash;

	}

	// The _index_of_deleted-th node is replaced by the _new_node (saving the order)
	void ReplaceNode(ExtremePoint _new_node, ExtremePoint _old_node)
	{

		bool erased = false;
		bool inserted = false;
		bool comeback = false;

		// Warning: not optimal insertion algorithm!
		list<ExtremePoint>::iterator nodeit;                   // Iterator for the nodes of facets
		
		for(nodeit = nodes.begin(); nodeit != nodes.end(); nodeit++)
		{

			if(!erased && nodeit->Equal(_old_node))
			{
				list<ExtremePoint>::iterator nextit = nodeit;
				nextit++;
				nodes.erase(nodeit);
				nodeit = nextit;
				erased = true;
				//comeback = true;
				if(nodeit == nodes.end())
				{
					break;
				}
			}

			// Searching a place for insertion with saving the order
			if(!inserted && ExtremePoint::PartialEPointOrder(*nodeit, _new_node))
			{
				nodes.insert(nodeit, _new_node);
				inserted = true;
			}

			// If both operations are completed
			if(erased && inserted)
			{
				return;
			}
		}

		if(!inserted)
		{
			nodes.insert(nodes.end(), _new_node);
		}
	}

	// Insertion of a node is made according to the order
	void InsertNode(ExtremePoint _new_node)
	{

		// Warning: not optimal insertion algorithm!
		list<ExtremePoint>::iterator nodeit;                   // Iterator for the nodes of facets
		
		for(nodeit = nodes.begin(); nodeit != nodes.end(); nodeit++)
		{

			// Searching a place for insertion with saving the order
			if(ExtremePoint::PartialEPointOrder(*nodeit, _new_node))
			{
				nodes.insert(nodeit, _new_node);
				return;
			}

		}

		nodes.insert(nodes.end(), _new_node);

	}

	void AdjustNodes()
	{
		//sort<ExtremePoint>(this->nodes.begin(), this->nodes.end(), &ExtremePoint::PartialEPointOrder);
	}

	// *** To hyperplane approach ***

	boost::dynamic_bitset<unsigned int> basic_scheme;  // The basic map for this facet (except d points with variable positions defining the facet) 
	
	list<int>     def_set;       // A set of $d$ points from the data cloud defining a position of the facet (without the absolute member)

	int           old_point;     // The neighbour from which a current facet is obtained
	Vector        old_normal;    // The neighbour's normal from which a current facet is obtained  
	Vector		  old_vector;    // A vector tossed off when generating the facet  

	int           got_point;     // An obtained point while generating a current facet

	int			  anchor;        // $k$ from (8) the number of points lower than the basis 

	float         abs_member;    // The absolute member in the facet's hyperplane equation 

	bool          pos_ready;     // Whether a position is calculated 
	bool          doubled;       // If it is the second facet for the current def_set
	


public:

	// Solving a linear system $\bmA \cdot \bmr = \bmb$ to get $\bmr$ = hypercoord
	void CalculateEquationH(vector<vector<float> > A, vector<float> b)
	{
		// Augmenting matrix and vector (to get non-vanishing solution)
		b[dim-1] = 1;
		vector<float> sinvec(dim, 1);
		A[dim-1] = sinvec;

		// The solution is obtained
		normalvec.coord = SolveLinearSystem(A, b);

		// Set a direction of the facet (relat. 0-point)

		pos_ready = true;    // Position calculated
	}

	vector<float> GiveVecBasis(vector<Vector> _defvecs, int _replace)
	{
		vector<float> b(dim, 0);

		_defvecs[_replace] = this->old_normal;

		// Augmenting matrix and vector (to get non-vanishing solution)
		b[dim-1] = 1;

		vector<vector<float> > A(dim);

		for(int i = 0; i < _defvecs.size(); i++)
		{
			A[i] = _defvecs[i].coord;
		}
		vector<float> sinvec(dim, 1);
		A[dim-1] = sinvec;

		// The solution is obtained
		return SolveLinearSystem(A, b);
	}

	// _replace - the number of a point in def. set to be replaced
	vector<float> GiveGenBasis(vector<Vector> _defvecs, int _replace)
	{
		/*
		int vecind = _replace;
		int grnum  = -1;
		do
		{
			grnum++;
			vecind -= this->cardinals[grnum] + 1;
		}while(vecind >= 0);

		vecind += this->cardinals[grnum] + 1;

		if(vecind == 0)
		{
			// Remove one vector
			_defvecs[_replace - grnum] = this->old_normal;
		}
		else if(vecind == this->cardinals[grnum])
		{
			// Remove one vector
			_defvecs[_replace - grnum - 1] = this->old_normal;
		}
		else
		{
			// Get sum of two vectors
			_defvecs[_replace - grnum - 1].Add(_defvecs[_replace - grnum]);

			// Replace a vector
			_defvecs[_replace - grnum] = this->old_normal;
		}
		*/

		vector<float> b(dim, 0);

		// Augmenting matrix and vector (to get non-vanishing solution)
		b[dim-1] = 1;

		vector<vector<float> > A(dim);

		for(int i = 0; i < _defvecs.size(); i++)
		{
			A[i] = _defvecs[i].coord;
		}
		vector<float> sinvec(dim, 1);
		A[dim-1] = sinvec;

		// The solution is obtained
		return SolveLinearSystem(A, b);
	}

	void CalculateAbsoluteMemberH(Point _node)
	{
		Vector vec(this->dim);
		vec.coord = _node.coord;

		this->abs_member = - Vector::ScalarMultiply(this->normalvec, vec);

	}

	float ClassifyPoint(Point _pnt)
	{
		Vector vec(this->dim);
		vec.coord = _pnt.coord;

		return Vector::ScalarMultiply(this->normalvec, vec) + this->abs_member;
	}

	void SetDirection()
	{
		Point nullp(dim);

		for(int i = 0; i < dim; i++)
		{
			nullp.coord[i] = 0;
		}

		// Selecting a proper halfspace

		if(ClassifyPoint(nullp) > 0)
		{
			this->normalvec.Reverse();
			this->abs_member = - this->abs_member;
		}

	}

	boost::dynamic_bitset<unsigned int> GetHashMapH(int _num)
	{

		boost::dynamic_bitset<unsigned int> hmap(_num);
		hmap.reset();

		list<int>::iterator wit;
		for(wit = def_set.begin(); wit != def_set.end(); wit++)
		{
			hmap[*wit] = 1;
		}

		return hmap;

	}

	// Get an index permutation of an arbitrary node of the facet
	vector<int> GetIndexPerm(int _num)
	{

		int num = _num;
		vector<int> arbind(num);

		boost::dynamic_bitset<unsigned int> hmap = GetHashMapH(num);

		int li = 0;
		int bi = this->anchor;
		int hi = this->anchor + this->dim;

		// Using info: in basic scheme "0" - for lower and basic; "1" - for higher. In hash-map  "1" - for basic; "0" - for others
		for(int i = 0; i < num; i++)
		{
			if(this->basic_scheme[i])
			{
				// Highers

				arbind[hi] = i;
				hi++;
			}
			else if(hmap[i])
			{
				// Basics

				arbind[bi] = i;
				bi++;
			}
			else
			{
				// Lowers

				arbind[li] = i;
				li++;
			}

		}

		return arbind;
	}


};

/*
class TruncatedFacet : public Facet
{
public:

	list<ExtremePoint> extra_nodes;

	TruncatedFacet(int _dim)
		: Facet(_dim)
	{
	}
};
*/

class HashTable
{
protected:

	int expbase;                               // Actually, number of points in a data cloud

	unsigned int ceiling;

	boost::dynamic_bitset<unsigned int> mass;

	// According to the indentification theorem, it should uniquely represent ${\cal X}_F$
	unsigned int HashCode(Facet _currfacet)
	{
		// Coding each ${\cal A}_l$ by the product and the sum of its points indices requires 2 int32 blocks
		boost::dynamic_bitset<unsigned int> hmap(_currfacet.anchors.size() * 32 * 3 + 32);
		hmap.reset();

		// Mapping info about all ${\cal A}_l$ sequentially
		// (this sequence encodes an order of ${\cal A}_l$s)
		// There can be many other ways of mapping!
		list<int>::iterator itset, itcard;
		int tmpbase = 0;
		for(itset = _currfacet.anchors.begin(), itcard = _currfacet.cardinals.begin();
			itset != _currfacet.anchors.end();
			itset++, itcard++, tmpbase += 3)
		{
			int currsum = 0, currprod = 1, currprod2 = 1;
			for(int i = *itset; i < *itset + *itcard; i++)
			{
				currsum += _currfacet.index_perm[i];
				currprod = currprod * _currfacet.index_perm[i];
				currprod2 = currprod2 * (_currfacet.index_perm[i]+1);
			}

			hmap.m_bits[tmpbase]   = currsum;
			hmap.m_bits[tmpbase+1] = currprod;
			hmap.m_bits[tmpbase+2] = currprod2;
		}

		// To tract a "doubled" case (when there is only 1 ${\cal A}$)
		if(_currfacet.anchors.size() == 1)
		{
			hmap.m_bits[tmpbase] = (_currfacet.normalvec.coord[0] > 0);
		}

		// Transforming the map into a number in the int32 interval
		// It is a general stuff (independent on the mapping method)
		// To understand the procedure see theory!

		const int prime = 2147483647;  // Mersenne prime number $2^31-1$

		boost::uint64_t psum = 0;
		for(int i = 0; i < hmap.m_bits.size(); i++)
		{
			boost::uint64_t toadd = ((boost::uint64_t)(hmap.m_bits[i]) << (i%31)) % prime;

			psum += toadd;
		}

		psum = psum % prime;

		// Now, psum == (hmap % prime). The result is strongly exact!

		// Narrowing interval for hashing (to reduce a hash table)

		return (unsigned int)psum % this->ceiling;

	}

	// According to the indentification theorem, it should uniquely represent ${\cal X}_F$
	// The same as HashCode(), just with a different mapping method (the 1st part)
	unsigned int HashCode_full(Facet _currfacet)
	{
		boost::dynamic_bitset<unsigned int> hmap(expbase * _currfacet.anchors.size() + 1);
		hmap.reset();

		// Mapping info about all ${\cal A}_l$ sequentially
		// (this sequence encodes an order of ${\cal A}_l$s)
		list<int>::iterator itset, itcard;
		int tmpbase = 0;
		for(itset = _currfacet.anchors.begin(), itcard = _currfacet.cardinals.begin();
			itset != _currfacet.anchors.end();
			itset++, itcard++, tmpbase += expbase)
		{
			int _strt = *itset;
			int _stp  = *itset + *itcard;
			for(int i = _strt; i < _stp; i++)
			{
				hmap[tmpbase + _currfacet.index_perm[i]] = 1;
			}
		}

		// To tract a "doubled" case (when there is only 1 ${\cal A}$)
		if(_currfacet.anchors.size() == 1)
		{
			//hmap[tmpbase] = (_currfacet.anchors.front() > (expbase - _currfacet.anchors.front() - _currfacet.cardinals.front()));
			hmap[tmpbase] = (_currfacet.normalvec.coord[0] > 0);
		}

		// Transforming the map into a number in the int32 interval
		// It is a general stuff (independent on the mapping method)
		// To understand the procedure see theory!

		const int prime = 2147483647;  // Mersenne prime number $2^31-1$

		boost::uint64_t psum = 0;
		for(int i = 0; i < hmap.m_bits.size(); i++)
		{
			boost::uint64_t toadd = ((boost::uint64_t)(hmap.m_bits[i]) << (i%31)) % prime;

			psum += toadd;
		}

		psum = psum % prime;

		// Now, psum == (hmap % prime). The result is strongly exact!

		// Narrowing interval for hashing (to reduce a hash table)

		return (unsigned int)psum % this->ceiling;

	}

public:
	vector<float> weight;

protected:

	vector<int> wstairs;

	unsigned int HashCode_ridgeCombin(kFace _ridge)
	{
		//::ResumeCumulTime();

		// Define how many indices are packed into one integer
		int pack = 2; int indspace = 32 / pack;
		boost::dynamic_bitset<unsigned int> hmap((expbase+pack-(expbase%pack)) * indspace);
		//boost::dynamic_bitset<unsigned int> hmap(expbase * _ridge.anchors.size() + (expbase - _ridge.Dim())*32);
		hmap.reset();

		// General positioning info
		for(int i = 0; i < this->expbase; i++)
		{
			hmap.m_bits[_ridge.index_perm[i] / pack] += (wstairs[i] << indspace * (_ridge.index_perm[i] % pack));
		}

		// Mapping info about all ${\cal A}_l$ 
		int alcount = -1;
		list<int>::iterator itset, itcard;
		for(itset = _ridge.anchors.begin(), itcard = _ridge.cardinals.begin();
			itset != _ridge.anchors.end();
			itset++, itcard++, alcount--)
		{
			int _strt = *itset;
			int _stp  = *itset + *itcard;
			for(int i = _strt; i < _stp; i++)
			{
				//hmap.m_bits[_ridge.index_perm[i]] = alcount;
				int toadd = ((alcount-wstairs[i]) << indspace * (_ridge.index_perm[i] % pack));
				hmap.m_bits[_ridge.index_perm[i] / pack] += toadd;
			}		
		}


		const int prime = 2147483647;  // Mersenne prime number $2^31-1$
		
		boost::uint64_t psum = 0;
		for(int i = 0; i < hmap.m_bits.size(); i++)
		{		
			// Scramble the hash
			hmap.m_bits[i] ^= -i;

			boost::uint64_t toadd = ((boost::uint64_t)(hmap.m_bits[i]) << (i%31)) % prime;

			psum += toadd;
		}

		psum = psum % prime;

		// Now, psum == (hmap % prime). The result is strongly exact!

		// Narrowing interval for hashing (to reduce a hash table)

		//::StopCumulTime();
		return (unsigned int)psum % this->ceiling;
	}

	unsigned int HashCode_ridge(kFace _ridge, int _unique)
	{
		//::ResumeCumulTime();
		boost::dynamic_bitset<unsigned int> hmap(expbase * _ridge.anchors.size() + 32);
		//boost::dynamic_bitset<unsigned int> hmap(expbase * _ridge.anchors.size() + (expbase - _ridge.Dim())*32);
		hmap.reset();

		// Mapping info about all ${\cal A}_l$ sequentially
		// (this sequence encodes an order of ${\cal A}_l$s)
		list<int>::iterator itset, itcard;
		int tmpbase = 32;
		for(itset = _ridge.anchors.begin(), itcard = _ridge.cardinals.begin();
			itset != _ridge.anchors.end();
			itset++, itcard++, tmpbase += expbase)
		{
			int _strt = *itset;
			int _stp  = *itset + *itcard;
			for(int i = _strt; i < _stp; i++)
			{
				hmap[tmpbase + _ridge.index_perm[i]] = 1;
			}
		}

		// Adding unique information for the given ridge
		hmap.m_bits[0] = _unique;

		/*
		int curind = tmpbase/32 + 1;

		list<int>::iterator iti, itic;
		list<int>::iterator itass;
		list<bool>::iterator itasst;

		int  afterI = 0;
		int elem;
		for(iti = _ridge.anchors.begin(), itic = _ridge.cardinals.begin(); iti != _ridge.anchors.end(); iti++, itic++)
		{
			// Gathering all type II clusters before current type I
			int typeIoccur = *iti;
			if(typeIoccur > afterI)
			{
				elem = _ridge.index_perm[afterI];

				for(int i = afterI+1; i < typeIoccur; i++)
				{
					if( this->weight[i] > this->weight[i-1] )
					{
						hmap.m_bits[curind] = elem;
						curind++;
						elem = _ridge.index_perm[i];
					}
					else
					{
						elem += _ridge.index_perm[i];
					}
				}

			}

			afterI = typeIoccur + *itic;
		}

		// Gathering all type II clusters after the last type I
		if(afterI < expbase)
		{
			elem = _ridge.index_perm[afterI];
		}
		for(int i = afterI + 1; i < expbase; i++)
		{
			if( this->weight[i] > this->weight[i-1] )
			{
				hmap.m_bits[curind] = elem;
				curind++;
				elem = _ridge.index_perm[i];
			}
			else
			{
				elem += _ridge.index_perm[i];
			}
		}
		*/

		// Transforming the map into a number in the int32 interval
		// It is a general stuff (independent on the mapping method)
		// To understand the procedure see theory!

		const int prime = 2147483647;  // Mersenne prime number $2^31-1$

		boost::uint64_t psum = 0;
		for(int i = 0; i < hmap.m_bits.size(); i++)
		{
			boost::uint64_t toadd = ((boost::uint64_t)(hmap.m_bits[i]) << (i%31)) % prime;

			psum += toadd;
		}

		psum = psum % prime;

		// Now, psum == (hmap % prime). The result is strongly exact!

		// Narrowing interval for hashing (to reduce a hash table)

		//::StopCumulTime();
		return (unsigned int)psum % this->ceiling;

	}

	unsigned int HashCode2(Facet _currfacet)
	{
		boost::dynamic_bitset<unsigned int> hmap(expbase);
		hmap.reset();
		list<int>::iterator wit;

		for(wit = _currfacet.def_set.begin(); wit != _currfacet.def_set.end(); wit++)
		{
			hmap[*wit] = 1;
		}

		boost::uint64_t hcode = 0;

		for(int i = hmap.m_bits.size()-1; i >= 0; i--)
		{
			hcode = (hcode << 32);                    // Shifting by 4 bytes (int32)
			hcode += (int)hmap.m_bits[i];
		}

		return hcode % this->ceiling;
	}

	boost::dynamic_bitset<unsigned int> mass_ridge1;
	boost::dynamic_bitset<unsigned int> mass_ridge2;

public:

	void Reset()
	{
		mass.reset();
		mass_ridge1.reset();
		mass_ridge2.reset();
	}

protected:
	void InitializeHash()
	{
		ceiling = mass_ridge1.size();
		Reset();
	}

	HashTable()
		:mass(1),
		mass_ridge1(83234597),
		mass_ridge2(83234597)
	{
		this->expbase = 0;
		InitializeHash();
	}

public:

	HashTable(int _expbase)
		:mass(1),
		mass_ridge1(83234597),
		mass_ridge2(83234597),
		wstairs(_expbase)
	{
		this->expbase = _expbase;
		InitializeHash();
	}

	void StartHashTable(vector<float> _weight)
	{
		this->weight = _weight;

		int staircount = 1;
		wstairs[0] = staircount;
		for(int i = 1; i < this->expbase; i++)
		{
			if(weight[i] > weight[i-1]) staircount++;
			wstairs[i] = staircount;
		}
	}

	void Mark(Facet _new_fac)
	{
		this->mass[HashCode_full(_new_fac)] = 1;
	}

	bool Check(Facet _currfacet)
	{
		return this->mass[HashCode_full(_currfacet)];
	}

	void MarkRidge(kFace _ridge/*, int _unique*/)
	{
		//unsigned int hscode = HashCode_ridge(_ridge, _unique);
		unsigned int hscode = HashCode_ridgeCombin(_ridge);
		
		if(this->mass_ridge1[hscode])
			this->mass_ridge2[hscode] = 1;
		else
			this->mass_ridge1[hscode] = 1;

	}

	bool CheckRidge(kFace _ridge/*, int _unique*/)
	{
		//return this->mass_ridge2[HashCode_ridge(_ridge, _unique)]; 
		return this->mass_ridge2[HashCode_ridgeCombin(_ridge)]; 
	}
};

#undef INCL_OBJECTS
#endif /*INCL_OBJECTS*/
