// Here we describe the main agents which realize the algorithm
// P.Bazovkin, 2009.

//#include <windows.h>		// Header File For Windows
#include <time.h>
//#include "SDL/SDL.h"
//#include <GL/gl.h>			// Header File For The OpenGL32 Library
//#include <GL/glu.h>			// Header File For The GLu32 Library
//#include <GL/glext.h>		// Header File For The Glaux Library

// Routine with the static elements
#define INCL_OBJECTS
#include "Objects.cpp"		// Header File For my private objects Library
//#undef INCL_OBJECTS

// General interface for algorithmic agents
class Agent
{
};


// ********************** Class for displaying the result *****************************************

	bool	keys[256];		// Array Used For The Keyboard Routine
	bool	active;			// Window Active Flag Set To TRUE By Default


class ResultAg : public Agent
{
public:

	bool	fullscreen;	    // Fullscreen Flag Set To Fullscreen Mode By Default

	list<Facet> trimmed_region;
	vector<Point> data_cloud;

	Point coord_center;         // Coordinates of the centroid

	void PrintResults(const char* _dir, const int _bound)
	{
		//ofstream os("TRegion.dat");
		
		char pathtofile[100]; 
		FILE *stream, *stream_ub;
		strcpy(pathtofile, _dir);
		stream    = fopen( strcat(pathtofile, "TRegion_vertices.dat"), "w" );
		strcpy(pathtofile, _dir);
		stream_ub = fopen( strcat(pathtofile, "TRegion_vertices_bound.dat"), "w" );
		
		//os << "The calculated trimmed region: (" << this->trimmed_region.size() << " facets) \n";
		//fprintf( stream, "The calculated trimmed region: (%d facets) \n", this->trimmed_region.size() );

		list<Facet>::iterator facetit;

		// For each facet
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
			// For printing only a boundary facets: the upper one if _bound == 1, and the lower one if _bound == -1. Don't write anything if _bound == 0
			bool writebound = 
				((_bound == 1) && !(facetit->normalvec.coord[0] < 0 || facetit->normalvec.coord[1] < 0 || facetit->normalvec.coord[2] < 0)) ||
				((_bound == -1) && !(facetit->normalvec.coord[0] > 0 || facetit->normalvec.coord[1] > 0 || facetit->normalvec.coord[2] > 0));
			
			//os << "( ";
			fprintf( stream,"( " );
			if(writebound) fprintf( stream_ub,"( " );

			list<ExtremePoint>::iterator expointit;

			// For each node of the facet
			for(expointit = facetit->nodes.begin(); expointit != facetit->nodes.end(); expointit++)
			{
				//os << "(";
				fprintf( stream,"(" );
				if(writebound) fprintf( stream_ub,"(" );

				for(int i = 0; i < expointit->coord.size(); i++)
				{
					//os << expointit->coord[i] + this->coord_center.coord[i] << ";";
					fprintf( stream,"%.3f;", expointit->coord[i]/* + this->coord_center.coord[i]*/ );
					if(writebound) fprintf( stream_ub,"%.3f;", expointit->coord[i]/* + this->coord_center.coord[i]*/ );
				}

				//os << ") ";
				fprintf( stream,") " );
				if(writebound) fprintf( stream_ub,") " );
			
			}
			
			//os << " )\n\n";
			fprintf( stream," )\n\n" );
			if(writebound) fprintf( stream_ub," )\n\n" );

		}

		//os.close();
		fclose( stream );
		fclose( stream_ub );

		return;

		//if(!boost::filesystem::create_directory("WMTR-temp")) return;

		ofstream tos("WMTR-temp/Timing.dat");
		tos << "Collected timestamps: \n";
		
		list<double>::iterator tit;
		for(tit = ::timestamps.begin(); tit != ::timestamps.end(); tit++)
		{
			tos << *tit << "\n";
		}

		tos.close();
	}

	void PrintResults(int _d, int _n, int _ind, int _indcount)
	{
		char filenm[200];
		sprintf(filenm, "Trimmed_region_%d_%d_%d.dat", _d, _n, _ind);

		ofstream os(filenm);
		
		os << "The calculated trimmed region: (" << this->trimmed_region.size() << " facets) \n";

		list<Facet>::iterator facetit;

		// For each facet
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
			os << "( ";

			list<ExtremePoint>::iterator expointit;

			// For each node of the facet
			for(expointit = facetit->nodes.begin(); expointit != facetit->nodes.end(); expointit++)
			{
				os << "(";

				for(int i = 0; i < expointit->coord.size(); i++)
				{
					os << expointit->coord[i]/* + this->coord_center.coord[i]*/ << ";";
				}

				os << ") ";
			}
			
			os << " )\n\n";
		}

		os.close();

		FILE* tos;
		tos = fopen( "Timing.dat", "a");
		fprintf(tos, "%d, %d: ", _d, _n);
		
		list<double>::iterator tit;
		for(tit = ::timestamps.begin(); tit != ::timestamps.end(); tit++)
		{
			fprintf(tos, "%f\t", *tit);
		}
		double tperfacet = ::timestamps.front() / this->trimmed_region.size();
		fprintf(tos, "- %f per facet (%d facets)\n", tperfacet, this->trimmed_region.size());

		if(_ind == 1)
		{
			::summtimes.push_back(::timestamps.front());
			::summtimes.push_back(tperfacet);
			::summtimes.push_back(this->trimmed_region.size());
		}
		else 
		{
			::summtimes[0] += ::timestamps.front();
			::summtimes[1] += tperfacet;
			::summtimes[2] += this->trimmed_region.size();
		}
		
		// Average results for the combination 
		if(_ind == _indcount)
		{
			::summtimes[0] = ::summtimes[0] / _indcount;
			::summtimes[1] = ::summtimes[1] / _indcount;
			::summtimes[2] = ::summtimes[2] / _indcount;
			
			fprintf(tos, "Average times for %d, %d: %f - %f per facet (%d facets)\n",  _d, _n, ::summtimes[0], ::summtimes[1], ::summtimes[2]);
			
			::summtimes.clear();

		}

		fclose(tos);
	}

	void PrintResultsHyperplanes(const char* _dir, const int _bound)
	{

		char pathtofile[100]; 
		strcpy(pathtofile, _dir);
		ofstream osv    ( strcat(pathtofile,"tmp_vrtheap.dat") );
		strcpy(pathtofile, _dir);
		ofstream ostr   ( strcat(pathtofile,"TRegion.dat") );
		strcpy(pathtofile, _dir);
		ofstream ostrbnd( strcat(pathtofile,"TRegion_bound.dat") );
		//ofstream os("Trimmed_region.dat");
		
		//os << "The calculated trimmed region: (" << this->trimmed_region.size() << " facets) \n";

		list<Facet>::iterator facetit;

		// For each facet
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
			// For printing only a boundary facets: the upper one if _bound == 1, and the lower one if _bound == -1. Don't write anything if _bound == 0
			bool writebound = 
				((_bound == 1) && !(facetit->normalvec.coord[0] < 0 || facetit->normalvec.coord[1] < 0 || facetit->normalvec.coord[2] < 0)) ||
				((_bound == -1) && !(facetit->normalvec.coord[0] > 0 || facetit->normalvec.coord[1] > 0 || facetit->normalvec.coord[2] > 0));
			
			list<ExtremePoint>::iterator expointit;

			// For each node of the facet
			for(expointit = facetit->nodes.begin(); expointit != facetit->nodes.end(); expointit++)
			{
				for(int i = 0; i < expointit->coord.size(); i++)
				{
					osv << expointit->coord[i]/* + this->coord_center.coord[i]*/ << " ";
				}

				osv << "\n";
			}

  			Vector shift;
			shift.coord = this->coord_center.coord;
			  // Writing hyperplanes
				for(int i = 0; i < facetit->normalvec.coord.size(); i++)
				{
					ostr << facetit->normalvec.coord[i] << " ";
					if(writebound) ostrbnd << facetit->normalvec.coord[i] << " ";
				}

				ostr << facetit->abs_member/*- Vector::ScalarMultiply(facetit->normalvec, shift)*/ << endl;
				if(writebound) ostrbnd << facetit->abs_member/*- Vector::ScalarMultiply(facetit->normalvec, shift)*/ << endl;
			
		}

		osv.close();
		ostr.close();
		ostrbnd.close();

	}

	void PrintResultsHyperplanes(int _d, int _n, int _ind, int _indcount, char* _type)
	{
		char filenm[200];
		sprintf(filenm, "Trimmed_region_%d_%d_%s_%d.dat", _d, _n, _type, _ind);

		ofstream os(filenm);
		
		os << "The calculated trimmed region: (" << this->trimmed_region.size() << " facets) \n";

		Vector shift;
		shift.coord = this->coord_center.coord;

		list<Facet>::iterator facetit;

		// For each facet
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
			os << "( ";

			for(int i = 0; i < facetit->normalvec.coord.size(); i++)
			{
				os << facetit->normalvec.coord[i] << ", ";
			}

			
			os << " ), abs.dist. = " << -(facetit->abs_member/* - Vector::ScalarMultiply(facetit->normalvec, shift)*/) << "\n";
		}

		os.close();

		FILE* tos; 
		sprintf(filenm, "Timing_%s.dat", _type);
		tos = fopen( filenm, "a");
		fprintf(tos, "%d, %d: ", _d, _n);
		
		list<double>::iterator tit;
		for(tit = ::timestamps.begin(); tit != ::timestamps.end(); tit++)
		{
			fprintf(tos, "%f\t", *tit);
		}
		double tperfacet = ::timestamps.front() / this->trimmed_region.size();
		fprintf(tos, "- %f per facet (%d facets)\n", tperfacet, this->trimmed_region.size());

		if(_ind == 1)
		{
			::summtimes.push_back(::timestamps.front());
			::summtimes.push_back(tperfacet);
			::summtimes.push_back(this->trimmed_region.size());
		}
		else 
		{
			::summtimes[0] += ::timestamps.front();
			::summtimes[1] += tperfacet;
			::summtimes[2] += this->trimmed_region.size();
		}
		
		// Average results for the combination 
		if(_ind == _indcount)
		{
			::summtimes[0] = ::summtimes[0] / _indcount;
			::summtimes[1] = ::summtimes[1] / _indcount;
			::summtimes[2] = ::summtimes[2] / _indcount;
			
			fprintf(tos, "Average times for %d, %d: %f - %f per facet (%d facets)\n",  _d, _n, ::summtimes[0], ::summtimes[1], ::summtimes[2]);
			
			::summtimes.clear();

		}

		fclose(tos);
	}

	vector<float> pmin;
	vector<float> pmax;
	vector<float> pdelta;

	// Only for dim = 3
	void FindBorders()
	{
		pmin.push_back(100000);
		pmin.push_back(100000);
		pmin.push_back(100000);

		pmax.push_back(-100000);
		pmax.push_back(-100000);
		pmax.push_back(-100000);

		list<Facet>::iterator facetit;

		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
			list<ExtremePoint>::iterator expointit;
			
			expointit = facetit->nodes.begin();

			for(int j = 0; j < 3; j++)
			{
				if(expointit->coord[0] < pmin[0])
				{
					pmin[0] = expointit->coord[0];
				}
				if(expointit->coord[1] < pmin[1])
				{
					pmin[1] = expointit->coord[1];
				}
				if(expointit->coord[2] < pmin[2])
				{
					pmin[2] = expointit->coord[2];
				}

				if(expointit->coord[0] > pmax[0])
				{
					pmax[0] = expointit->coord[0];
				}
				if(expointit->coord[1] > pmax[1])
				{
					pmax[1] = expointit->coord[1];
				}
				if(expointit->coord[2] > pmax[2])
				{
					pmax[2] = expointit->coord[2];
				}

				expointit++;
			}
		}

		pdelta.push_back(pmax[0] - pmin[0]);
		pdelta.push_back(pmax[1] - pmin[1]);
		pdelta.push_back(pmax[2] - pmin[2]);
	}
	
	ResultAg()
	{

	}

	//LRESULT	CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);	// Declaration For WndProc

	
private:
	// technical variables
	int flborder;
	bool markfirst;

public:


	void GLSceneToExport(const char* _dir, const int _bound)									// Here's Where We Do All The Drawing
	{
		// Only for d=3
		if(this->coord_center.coord.size() != 3) return;

		int currcol[3];

		//// ********************* Display points from the data cloud **************************

		
		char pathtofile[100]; 
		strcpy(pathtofile, _dir);
		ofstream osvisd( strcat(pathtofile, "tmp_vis_dt.dat"));

		currcol[0] = 255;
		currcol[1] = 0;
		currcol[2] = 0;

		vector<Point>::iterator pointit;

		float fat = 0.1f;

		// Typing coordinates of the mean at the first place
		//osvisd << this->coord_center.coord[0] << " " << this->coord_center.coord[1] << " " << this->coord_center.coord[2] << endl;
		osvisd << 0 << " " << 0 << " " << 0 << endl;

		for(pointit = this->data_cloud.begin(); pointit != this->data_cloud.end(); pointit++)
		{
			// Displaying as spheres	
			osvisd << pointit->coord[0] << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;

			// Displaying as pyramides	
			/* osvisd << pointit->coord[0] + fat << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] + fat << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] - fat << endl;					

				osvisd << pointit->coord[0] + fat << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] + fat << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] + fat << endl;					
				
				osvisd << pointit->coord[0] + fat << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] - fat << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] + fat << endl;					

				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] + fat << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] - fat << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] + fat << endl; */					
		}

		osvisd.close();


		// ***********************************************************************************


		// ********************* Display the trimmed region **********************************

		strcpy(pathtofile, _dir);
		ofstream osvisf( strcat(pathtofile, "tmp_vis_fc.dat"));

		list<Facet>::iterator facetit;

		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{
			// For printing only a boundary facets: the upper one if _bound == 1, and the lower one if _bound == -1. Don't write anything if _bound == 0
			if(((_bound == 1) && (facetit->normalvec.coord[0] < 0 || facetit->normalvec.coord[1] < 0 || facetit->normalvec.coord[2] < 0)) ||
				((_bound == -1) && (facetit->normalvec.coord[0] > 0 || facetit->normalvec.coord[1] > 0 || facetit->normalvec.coord[2] > 0)))
				continue;

			if( !facetit->truncated /*&&  ++stop<flborder*/  )
			{

		
				float peak = max(fabs(facetit->normalvec.coord[0]), max(fabs(facetit->normalvec.coord[1]), fabs(facetit->normalvec.coord[2])));
				float gr_sigma = 2;
				float gr_mu    = 0.5f;
				float colx = -facetit->normalvec.coord[0] / (gr_sigma * peak) + gr_mu;
				float coly = -facetit->normalvec.coord[1] / (gr_sigma * peak) + gr_mu;
				float colz = -facetit->normalvec.coord[2] / (gr_sigma * peak) + gr_mu;
				//glColor3f(colx, coly, colz);

				list<ExtremePoint>::iterator expointit;
				
				expointit = facetit->nodes.begin();
				osvisf << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;		
				
				expointit++;
				osvisf << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;					
				
				expointit++;
				osvisf << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;
			}
		}

	

		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
			// For printing only a boundary facets: the upper one if _bound == 1, and the lower one if _bound == -1. Don't write anything if _bound == 0
			if(((_bound == 1) && (facetit->normalvec.coord[0] < 0 || facetit->normalvec.coord[1] < 0 || facetit->normalvec.coord[2] < 0)) ||
				((_bound == -1) && (facetit->normalvec.coord[0] > 0 || facetit->normalvec.coord[1] > 0 || facetit->normalvec.coord[2] > 0)))
				continue;


			if( facetit->truncated /*&& ++stop<flborder*/ )
			{

				float peak = max(fabs(facetit->normalvec.coord[0]), max(fabs(facetit->normalvec.coord[1]), fabs(facetit->normalvec.coord[2])));
				float gr_sigma = 2;
				float gr_mu    = 0.5f;
				float colx = -facetit->normalvec.coord[0] / (gr_sigma * peak) + gr_mu;
				float coly = -facetit->normalvec.coord[1] / (gr_sigma * peak) + gr_mu;
				float colz = -facetit->normalvec.coord[2] / (gr_sigma * peak) + gr_mu;

				list<ExtremePoint>::iterator expointit;
				list<ExtremePoint>::iterator expointitback;
				
				//glBegin(GL_POLYGON);

				//glColor3f(colx, coly, colz);

				vector<vector<float> > cvec(facetit->nodes.size());
				int j;
				for(j = 0, expointit = facetit->nodes.begin(); expointit != facetit->nodes.end(); expointit++, j++)
				{
					cvec[j] = expointit->coord;

				}

				// Hexagonal
				if(cvec.size() == 6)
				{
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[1])[0] << " " << (cvec[1])[1] << " " << (cvec[1])[2] << endl;
					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;

					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;
					osvisf << (cvec[3])[0] << " " << (cvec[3])[1] << " " << (cvec[3])[2] << endl;
					osvisf << (cvec[4])[0] << " " << (cvec[4])[1] << " " << (cvec[4])[2] << endl;

					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[5])[0] << " " << (cvec[5])[1] << " " << (cvec[5])[2] << endl;
					osvisf << (cvec[4])[0] << " " << (cvec[4])[1] << " " << (cvec[4])[2] << endl;

					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[4])[0] << " " << (cvec[4])[1] << " " << (cvec[4])[2] << endl;
				}
				// Quad
				else
				{
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[1])[0] << " " << (cvec[1])[1] << " " << (cvec[1])[2] << endl;
					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;

					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;
					osvisf << (cvec[3])[0] << " " << (cvec[3])[1] << " " << (cvec[3])[2] << endl;
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
				}

				
			}
		}

		osvisf.close();
		

		// ***********************************************************************************
		
		
		// ********************* Display edges of the trimmed region *************************

		strcpy(pathtofile, _dir);
		ofstream osvisr( strcat(pathtofile, "tmp_vis_rg.dat") );
		
		//glColor3f(0.4f,0.4f,0.4f);						
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			

			list<ExtremePoint>::iterator expointit, expointitn;
			expointit = facetit->nodes.begin();
			expointitn = expointit;
			expointitn++;
			for(; expointitn != facetit->nodes.end(); expointit++, expointitn++)
			{
				osvisr << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;		
				osvisr << expointitn->coord[0] << " " << expointitn->coord[1] << " " << expointitn->coord[2] << endl;		
			}
			expointitn = facetit->nodes.begin();
			osvisr << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;		
			osvisr << expointitn->coord[0] << " " << expointitn->coord[1] << " " << expointitn->coord[2] << endl;		
				
		}

		osvisr.close();

		// ***********************************************************************************


	}




};
//*************************************************************************************************


// ********************** Class for building a trimmed region *************************************

class ProcessAg : public Agent
{
public:

	Permutation perm;             // Current permutation

	Point initial_center;         // Coordinates of the centroid

	string WMTD_type;              // Notion of a WMTR

protected:

	list<Facet> trimmed_region; // Trimmed region: set of the facets

	float alpha;                  // Zonoid depth

	list<Facet> queue;            // Queue of the not processed facets

	int   dim  ;                  // Dimension of the cloud
	int   num  ;                  // Number of points

	Vector p;                     // Current direction for support function

	int border_index;             // Index of the border point (affected by the depth)
	vector<float> weight;         // Weigts vector (non-decreasing function on the permutation)

	HashTable ohash_table;        // "Optimistic" hash table
	
	void NormalizeWeights()
	{
		// Smoothing too small weights
		for(int i = 0; i < num; i++)
		{
			if(this->weight[i] < 1.0f / (this->num * 100))
				this->weight[i] = 0.0f;
		}

		float sum_of_weights = 0;
		for(int i = 0; i < num; i++)
		{
			sum_of_weights += this->weight[i];
		}

		// If there is an incorrect alpha parametere
		if(sum_of_weights < 0.01f)
		{
			 this->weight[num-1] = 1;
			 sum_of_weights = 1;
		}

		for(int i = 0; i < num; i++)
		{
			this->weight[i] = this->weight[i] / sum_of_weights;
		}
	}

public:

	// Generates weight vector basing on the type of a WMTR
	void WeightsGenerator()
	{
		if(this->WMTD_type == "zonoid")
		{
			int zonbord = num - floor(alpha * num) - 1;
			this->border_index = zonbord;

			for(int i = 0; i < this->num; i++)
			{
				if( i<zonbord )
				{
					this->weight[i] = 0;
				}
				else if( i==zonbord)
				{
					this->weight[i] = (alpha * num - floor(alpha * num)) / (alpha * num);
				}
				else
				{
					this->weight[i] = 1.0f / (alpha * num);
				}
			}
		}
		else if(this->WMTD_type == "ECH") 
		{
			int beta = 1.0f / this->alpha;

			for(int j = num-1; j >= 0; j--)
			{
				if(j < beta-1)
				{
					this->weight[j]= 0.0f;
				}
				else if(j == num-1)
				{
					//this->weight[i]= 1.0f / Comb(num, beta);
					this->weight[j] = ((float)beta)/num;
				}
				else
				{
					this->weight[j]= weight[j+1] * (j+2-beta) / (j+1);
				}
			}
		}
		else if(this->WMTD_type == "ECH*") 
		{
			float beta = 1.0f / this->alpha;

			for(int i = 0; i < num; i++)
			{
				this->weight[i]= ( pow(i+1, beta) - pow(i, beta) ) / pow(num, beta);
			}
		}
		else if(this->WMTD_type == "geometrical") 
		{
			float gmult = (1-alpha) / (1-pow(alpha, num));

			for(int i = 0; i < num; i++)
			{
				this->weight[i]= pow(alpha, num-i-1) * gmult;
			}
		}
		else if(this->WMTD_type == "general") 
		{
		}
		else
		{
			int zonbord = num - floor(alpha * num) - 1;
			this->border_index = zonbord;

			for(int i = 0; i < this->num; i++)
			{
				if( i<zonbord )
				{
					this->weight[i] = 0;
				}
				else if( i==zonbord)
				{
					this->weight[i] = (alpha * num - floor(alpha * num)) / (alpha * num);
				}
				else
				{
					this->weight[i] = 1.0f / (alpha * num);
				}
			}
		}
	}

	// Use only this constructor for this class
	ProcessAg(string _type, float _depth, int _dim, int _num, vector<Point> _cloud, vector<float> _manweights):
																p(_dim),
																initial_center(_dim),
																ohash_table(_num),
																weight(_num)
	{
		this->WMTD_type = _type;
		
		dim   = _dim;
		num   = _num;

		alpha = _depth;

		// ************** Normalization of the data ************

		vector<Point>::iterator it;

		// Calculate centroid
		for(it = _cloud.begin(); it != _cloud.end(); it++)
		{
			for(int j=0; j<dim; j++)
			{
				initial_center.coord[j] += it->coord[j];

			}

		}

		for(int j=0; j<dim; j++)
		{
			initial_center.coord[j] = initial_center.coord[j] / num;

		}

		//// Set centroid to 0
		//for(it = _cloud.begin(); it != _cloud.end(); it++)
		//{
		//	for(int j=0; j<dim; j++)
		//	{
		//		it->coord[j] = it->coord[j] - initial_center.coord[j];

		//	}

		//}
		// ************************************************************


		perm.points.assign(_cloud.begin(), _cloud.end());

		// Ordering points according to DEFAULT support vector
		perm.Support(p);

		if(this->WMTD_type == "general")
			weight = _manweights;
		else
			WeightsGenerator();

		NormalizeWeights();

		this->ohash_table.StartHashTable(this->weight);

	}

protected:

	// Working with the hash-table
	
	list<boost::dynamic_bitset<unsigned int> > hash_table;
	
	static bool HashOrder(boost::dynamic_bitset<unsigned int> bs1, boost::dynamic_bitset<unsigned int> bs2)
	{
		//::ResumeCumulTime();
		// Warning: not optimal procedure
		int bsrange = bs1.size();
		for(int i = 0; i < bsrange; i++)
		{
			if(bs1[i] > bs2[i])
			{
				//::StopCumulTime();
				return true;
			}
			else if(bs1[i] < bs2[i])
			{
				//::StopCumulTime();
				return false;
			}
		}
		//::StopCumulTime();
		return true;

		/*int byterange = bs1.m_bits.size();
		for(int i = 0; i < byterange; i++)
		{
			if((int)bs1.m_bits[i] > (int)bs2.m_bits[i])
				return true;
			else if((int)bs1.m_bits[i] < (int)bs2.m_bits[i])
				return false;
		}
		return true;*/
	}
	
	static bool ReversedHashOrder(boost::dynamic_bitset<unsigned int> bs1, boost::dynamic_bitset<unsigned int> bs2)
	{
		return !HashOrder(bs1, bs2);
	}

	void MarkByHash(Facet _new_fac)
	{
		// Warning: not optimal insertion algorithm!
		list<boost::dynamic_bitset<unsigned int> >::iterator hashit;

		bool inserted = false;
		
		boost::dynamic_bitset<unsigned int> curr_hash = _new_fac.GetHashMap();

		if(this->hash_table.size() == 0)
		{
			this->hash_table.push_back(curr_hash);
		}
		else
		{
			for(hashit = this->hash_table.begin(); hashit != this->hash_table.end(); hashit++)
			{
				// Searching a place for insertion with saving the order
				if(ProcessAg::HashOrder(*hashit, curr_hash))
				{
					this->hash_table.insert(hashit, curr_hash);
					inserted = true;
					return;
				}
			}
			
			if(!inserted)
			{
				this->hash_table.insert(hash_table.end(), curr_hash);
			}
		}
	}

	// NOT key function (see ContainsFacet) 
	bool CheckInHash(Facet _new_fac)
	{
		return binary_search(this->hash_table.begin(), this->hash_table.end(), _new_fac.GetHashMap(), &ProcessAg::ReversedHashOrder);
	}

private:
	
	// Factorial
	boost::uint64_t Fact(int _n)
	{
		if (_n <= 1)
		{
			return 1;
		}
		else
		{
			return _n*Fact(_n-1);
		}
	}

	// Number of combinations _n by _b
	int Comb(int _n, int _b)
	{
		return Fact(_n)/(Fact(_n-_b) * Fact(_b)); 
	}

protected:

	// Argument: index permutation of extr point
	// Result  : coordinates of this extr point
	vector<float> CurrentExtremePoint(vector<int> _index_perm)
	{
		Point extreme(this->dim);

		for(int j = 0; j < this->num; j++)
		{
			for(int i = 0; i<this->dim; i++)
			{
				extreme.coord[i] += this->perm.points[_index_perm[j]].coord[i] * this->weight[j];
			}

		}

		return extreme.coord;


	}

	vector<float> CurrentExtremePoint_gener(vector<int> _index_perm, Permutation _perm)
	{
		Point extreme(this->dim);

		for(int j = 0; j < this->num; j++)
		{
			for(int i = 0; i<this->dim; i++)
			{
				extreme.coord[i] += _perm.points[_index_perm[j]].coord[i] * this->weight[j];
			}

		}

		return extreme.coord;

	}

	// Neighbouring vertices
	boost::dynamic_bitset<unsigned int> FindNeighbours(ExtremePoint _epoint)
	{
		// The resulting array - corresponding to index permutation of the extr point
		boost::dynamic_bitset<unsigned int> res(this->num);

		// Creating the bundle of vectors to each point in the data cloud from the boundary point
		vector<Vector> bundle(this->num - 1);
			
		//_epoint.coord = CurrentExtremePoint(_epoint.index_perm);

		// Creating bundle of (n-1) vectors
		int vect_counter = 0;
		for(int i = 0; i < this->num; i++)
		{
			Vector vec(this->dim);
			
			// !!!!!!!! Warning!!! Index changed (after border_index)
			if(i < this->border_index)
			{
				vec.FromPointToPoint( perm.points[_epoint.index_perm[i]] , perm.points[_epoint.index_perm[this->border_index]]);
				bundle[vect_counter] = vec;
				vect_counter++ ;
			}
			else if(i > this->border_index)
			{
				vec.FromPointToPoint( perm.points[_epoint.index_perm[this->border_index]] , perm.points[_epoint.index_perm[i]]);
				bundle[vect_counter] = vec;
				vect_counter++ ;
			}

		}

//#define MY_TESTING
#ifdef MY_TESTING
		
		Vector yy(this->dim);
		
		yy.coord[0] = 1;
		yy.coord[1] = 2;
		yy.coord[2] = -1;
		bundle[0] = yy;

		yy.coord[0] = 0;
		yy.coord[1] = 1;
		yy.coord[2] = 1;
		bundle[1] = yy;

		yy.coord[0] = 2;
		yy.coord[1] = -1;
		yy.coord[2] = 2;
		bundle[2] = yy;

		yy.coord[0] = 1;
		yy.coord[1] = 2;
		yy.coord[2] = 1;
		bundle[3] = yy;

#endif

		// Now we have a bundle of (n-1) vectors codirected with extrpointvector (order according to index_perm)

		float epsilon = 0.001;                  // The precision for the floating point operations

		// Procedure of testing neighbourhood for each point
		for(int i = 0; i < this->num - 1; i++)
		{
			bool criterion      = false;        // Whether current point is a neighbour
			bool second_attempt = false;

			do
			{
			
				// Forming the matrix for simplex method
				vector<Vector> curr_bundle(this->num-1);
				curr_bundle = bundle;


				/* The method with reversing the current element 
				curr_bundle[i].Reverse();
				curr_bundle.erase(curr_bundle.begin() + this->border_index);     // !!Warning: index changes
				*/

				// ********* The way based on the simplex method ***********

				// Searching the first basic feasible solution basing on ex
				
				float M = 128;
				
				// Initializing
				Vector c(this->dim);
				c = curr_bundle[i];

				// Mirroring some dimensions of the vectors to transform the vector c to a nonnegative one
				for(int k = 0; k < dim; k++)
				{
					if(c.coord[k] < 0)
					{
						for(int q = 0; q < num-1; q++)
						{
							curr_bundle[q].coord[k] = -curr_bundle[q].coord[k];

						}

						c.coord[k] = -c.coord[k]; 
					}
				}

				// Creating vector of 1
				vector<float> ones(dim,1);
				Vector univect(dim);
				univect.coord = ones;
				if(second_attempt) univect.Reverse();

				curr_bundle.push_back(univect);

				// Memorizing the checksum for using in the criterion
				float checksum = curr_bundle[i].SumOfElements() * M;

				// Transforming matrix to the universal form
				for(int q = 0; q < dim; q++)
				{
					vector<float> tmpvec(dim, 0);
					tmpvec[q] = 1;
					
					Vector zusaetz_vektor(dim);
					zusaetz_vektor.coord = tmpvec;
					curr_bundle.push_back(zusaetz_vektor);
				}
				
				// Forming vector b according to the system (8)
				vector<float> b(num + dim + 1,0);	        // (n+1) element is the objective value in the dual task			
				
				for(int q = 0; q < num; q++)
				{
					b[q] = -curr_bundle[q].SumOfElements() * M;

				}
				b[num-1] += -1;
				
				// Iterations of the simplex method for the dual task 

				int main_x;

				// Do while the minus in the objective line is not found
				do
				{

					// The horizontal coordinate of the column of the selected non-basic variable (if remains =-1, then there is no such column)
					main_x = -1;

					for(int k = 0; k < num + dim; k++)
					{
						if( b[k] < -epsilon )                              // Checking if b[k] < 0
						{
							main_x = k;
							break; // Is found 
						}

					}

					if(main_x != -1)
					{

						// Searching coordinate of the basic variable to take of the basis
						int main_y = -1;

						float min_proportion = 100000;
						for(int k = 0; k < dim; k++)
						{
							float proportion = min_proportion;
							if(curr_bundle[main_x].coord[k] > 0)
								proportion = c.coord[k] / curr_bundle[main_x].coord[k];
							
							// Warning! - endless loop is possible
							if(proportion < min_proportion && proportion >= -epsilon)        // ">= -epsilon" is equivalent to ">= 0" 
							{
								min_proportion = proportion;
								main_y         = k;
							}

						}

	#ifdef MY_TESTING

	#endif

						// The solution is not limited!!!
						if(main_y < 0)
						{
							// If the dual task has unlimited solution -> the initial task has no feasible solution
							// Such result is possible only in one from 2 attempts

							//MessageBox(NULL,"Potential error the simplex method - solution is not limited", "Would you like to continue?",MB_YESNO|MB_ICONQUESTION);
							criterion = false;
							b[num + dim] = checksum + 1;   // Warning! It is made to fix (b[num + dim] != checksum)
							break;
							
							// return;
						}

						// Coordinates of the key element for the current simplex method loop are found - (main_y, main_x)

						for(int k = 0; k < dim; k++)
						{
							// For not key rows
							if(k != main_y)
							{
								float tmp_coeff = curr_bundle[main_x].coord[k] / curr_bundle[main_x].coord[main_y];
								
								for(int q = 0; q < num + dim; q++)
								{
									curr_bundle[q].coord[k] -= tmp_coeff *  curr_bundle[q].coord[main_y];
								}
								
								c.coord[k] -= tmp_coeff *  c.coord[main_y];
							}
							// For the key row
							else
							{
								float tmp_coeff = 1 / curr_bundle[main_x].coord[main_y];

								for(int q = 0; q < num + dim; q++)
								{
									curr_bundle[q].coord[main_y] = tmp_coeff * curr_bundle[q].coord[main_y];
								}
								
								c.coord[main_y] = tmp_coeff * c.coord[main_y];

								// Correcting objective line coeffs
								tmp_coeff = b[main_x];
								
								for(int q = 0; q < num + dim; q++)
								{
									b[q] -= tmp_coeff * curr_bundle[q].coord[main_y];
								}

								b[num + dim] -= tmp_coeff * c.coord[main_y];
							}

						}
					}

				}while(main_x != -1);

				// Having the value of the objective function in the optimal solution, check the criterion
				if( fabs(b[num + dim] - checksum) < epsilon )         // Checking equality
				{
					criterion      = true;
					second_attempt = false;
				}
				else
				{
					criterion      = false;
					second_attempt = !second_attempt;
				}

			}while(second_attempt);

			if(i < this->border_index)
			{
				res[i] = criterion;
			}
			else
			{
				res[i+1] = criterion;
			}


			// We've got a set of vectors curr_bundle, which enable us to test neighbourhood: if all vectors leave codirected
			
			// !!!! Warning !!!! Now we are using a stub - only for 3d
			// !!!!!!!!!! Start stub
			bool flag = 0;
			// !!!!!!!!!! End stub
		}

		// !!! Warning - senseless code
		res[this->border_index] = false;


		return res;
	}

private:

	// Returns whether current facet has already been included into the trimmed region 
	bool ContainsFacet(Facet _currfacet)
	{
		return binary_search(this->hash_table.begin(), this->hash_table.end(), _currfacet.GetHashMap(), &ProcessAg::ReversedHashOrder);
	}

protected:

	// Gets an extrpoint and finds a facet with it, searching among neighbours from the specified _part (left or right from the border_index)
	Facet FindFirstFacetInBlock(ExtremePoint _start_node, int _zone)
	{
		// We create d points for the first facet. Then find their neighbours
		Facet ffacet(this->dim);
		_start_node.neighs = FindNeighbours(_start_node);
		ffacet.nodes.push_back(_start_node);

		boost::dynamic_bitset<unsigned int> fneigh = _start_node.neighs;

		int search_start;
		int search_stop;

		// If neighs are being searched in the first part of the permutation
		if(_zone == 0)
		{
			search_start = 0;
			search_stop  = this->border_index;
		}
		// ... in the second part of ...
		else
		{
			search_start = this->border_index;
			search_stop  = this->num;
		}

		int v = search_start;
		int found_neighs   = 0;
		bool all_traversed = false;

		// Searching (dim-1) neighbours 
		while((found_neighs < this->dim - 1) && !all_traversed)
		{
			// Finding (d-1) neighs of the first point and adding to ffacet
			while( fneigh[v] == 0 )
			{
				v++;                  // Skip all not neigh's
				
				if (v == search_stop)
				{
					all_traversed = true;
					break;
				}
			}

			if(all_traversed)
			{
				continue;
			}

			// Creating new extr point
			ExtremePoint node(this->dim, this->border_index);
			
			// Swapping indices of two elements
			node.index_perm = _start_node.index_perm;
			int temp                            = node.index_perm[v];
			node.index_perm[v]                  = node.index_perm[this->border_index];
			node.index_perm[this->border_index] = temp;

			node.coord = CurrentExtremePoint(node.index_perm);

			node.neighs = FindNeighbours(node);
			//ffacet.nodes.push_back(node);
			ffacet.InsertNode(node);
			found_neighs++;

			// Obtain common neighbours with the new node
			for(int q = search_start; q < search_stop; q++)
			{
				if(q == v)
				{
					fneigh[this->border_index] &= node.neighs[q];    // Warning: operator &=
				}
				else if(q == this->border_index)
				{
					fneigh[v] &= node.neighs[q];    // Warning: operator &=
				}
				else
				{
					fneigh[q] &= node.neighs[q];    // Warning: operator &=
				}
			}
			
			v++;
			
			if (v == search_stop)
			{
				all_traversed = true;
				continue;
			}
		}

		if(found_neighs < this->dim - 1)
		{
			Facet empty_fac(this->dim);
			return empty_fac;           // First facet is not found!
		}

		// Sorting nodes in the first facet according to the extreme point partial order
		//sort(ffacet.nodes.begin(), ffacet.nodes.end(), &ExtremePoint::PartialEPointOrder);

		// Assigning zone to the found facet
		ffacet.zone = _zone;
		
		return ffacet;
	}

public:

	list<Facet> Compute()
	{
		
		// ************************** Create the first extr point *****************
		
		vector<int> init_index_perm(this->num);
		for(int i=0; i<this->num; i++)
		{
			init_index_perm[i] = i;
		}

		ExtremePoint fnode(this->dim, this->border_index);
		fnode.index_perm = init_index_perm;
		fnode.coord      = CurrentExtremePoint(fnode.index_perm);

		// ************************************************************************

		// ************************** Create the first facet **********************

		// We create d points for the first facet. Then find their neighbours
		Facet ffacet = FindFirstFacetInBlock(fnode, 0);

		if(ffacet.nodes.empty())
		{
			ffacet = FindFirstFacetInBlock(fnode, 1);
		}

		// Mark facet by putting its hash-code to the hash-table
		this->MarkByHash(ffacet);

		// Placing the facet to the algorithms queue
		this->queue.push_back(ffacet);
		
		// ************************************************************************

		// ************************** Create all the facets ***********************
				
		do
		{
			// Take the upper facet
			Facet currfacet = queue.front();
			
			if( !currfacet.truncated )
			{

				list<ExtremePoint>::iterator nodeit;                   // Iterator for the nodes of facets

				int search_start;
				int search_stop;

				// If neighs are being searched in the first part of the permutation
				if(currfacet.zone == 0)
				{
					search_start = 0;
					search_stop  = this->border_index;
				}
				// ... in the second part of ...
				else
				{
					search_start = this->border_index;
					search_stop  = this->num-1;
				}

				//// Using lists of neighs for the different nodes of the facet obtain common neighs

				// Swapping every point sequently, except basic point (not to return to the parent-facet)
				//for(int i = 0; i < this->dim; i++)
				for(nodeit = currfacet.nodes.begin(); nodeit != currfacet.nodes.end(); nodeit++)
				{
					// Don't go back: generating (d-1) ridges
					// !!!Warning!!! for the first facet this condition is ever met (all is OK) but inexplicitly - error possible 
					if( nodeit->index_perm != currfacet.basic_perm )
					{
						boost::dynamic_bitset<unsigned int> common_neigh(num);

						// We create a mask for conjunction of neigh sets 
						//(1 - for all points in the zone (except main point of nodeit -> hence we need to find only the single 1) and 0 - for others)
						common_neigh.reset();
						for(int k = search_start; k <= search_stop; k++)
						{
							common_neigh[nodeit->index_perm[k]] = true;
						}
						common_neigh[nodeit->index_perm[this->border_index]] = false;

						list<ExtremePoint>::iterator tmpnodeit;
						//for(int j = 0; j < this->dim; j++)
						for(tmpnodeit = currfacet.nodes.begin(); tmpnodeit != currfacet.nodes.end(); tmpnodeit++)
						{
							// Skip erased point
							if( tmpnodeit != nodeit )
							{
								// Operations for every from (d-1) points of the ridge

								for(int k = search_start; k <= search_stop; k++)
								{
									common_neigh[tmpnodeit->index_perm[k]] &= tmpnodeit->neighs[k];
								}
							}

						}

						// Now in vector common_neigh we have a single 1 - that is our sought-for point

						// Inserting a new neighbour extr point to the facet
						
						int t = common_neigh.find_first();

						if(t == -1)
						{
							// If there is no neighbour, i.e. this ridge is a border with a truncated facet

							Facet truncfacet(this->dim, true);

							// Searching neighbours in the opposite zone

							truncfacet.zone = currfacet.zone;

							int trunc_start;
							int trunc_stop;

							// If neighs are being searched in the first part of the permutation
							if( (1 - currfacet.zone) == 0 )      // Zone is opposite
							{
								trunc_start = 0;
								trunc_stop  = this->border_index;
							}
							// ... in the second part of ...
							else
							{
								trunc_start = this->border_index;
								trunc_stop  = this->num-1;
							}

							boost::dynamic_bitset<unsigned int> trunc_neigh(num);

							// We create a mask for conjunction of neigh sets 
							trunc_neigh.reset();
							for(int k = trunc_start; k <= trunc_stop; k++)
							{
								trunc_neigh[nodeit->index_perm[k]] = true;
							}
							trunc_neigh[nodeit->index_perm[this->border_index]] = false;

							for(tmpnodeit = currfacet.nodes.begin(); tmpnodeit != currfacet.nodes.end(); tmpnodeit++)
							{
								// Skip erased point
								if( tmpnodeit != nodeit )
								{
									// Operations for every from (d-1) points of the ridge

									for(int k = trunc_start; k <= trunc_stop; k++)
									{
										trunc_neigh[tmpnodeit->index_perm[k]] &= tmpnodeit->neighs[k];
									}
								}

							}

							int trn = trunc_neigh.find_first();

							if(trn == -1)
							{
								// !!! Error - impossible!
							}

							// Searching extra nodes

							// !!! Stub !!!

							if(true)
							{
								//truncfacet.basic_scheme.reset();

								//for(tmpnodeit = currfacet.nodes.begin(); tmpnodeit != currfacet.nodes.end(); tmpnodeit++)
								//{
								//	// Skip erased point
								//	if( tmpnodeit != nodeit )
								//	{
								//		truncfacet.basic_scheme[tmpnodeit->index_perm[this->border_index]] = true;
								//	}
								//}

								//truncfacet.basic_scheme[trn] = true;

								truncfacet.nodes = currfacet.nodes;

								// Find index of the point trn in any point of the current ridge
								if(nodeit != currfacet.nodes.begin())
								{
									tmpnodeit = currfacet.nodes.begin();
								}
								else
								{
									tmpnodeit = --(currfacet.nodes.end());
								}

								int r   = 0;
								while(tmpnodeit->index_perm[r] != trn)
								{
									r++;
								}

								ExtremePoint neighbour = *tmpnodeit;

								// Swapping elements in neighbour
								int temp                                 = neighbour.index_perm[r];
								neighbour.index_perm[r]                  = neighbour.index_perm[this->border_index];
								neighbour.index_perm[this->border_index] = temp;

								// Compute coordinates of the new extr point
								neighbour.coord = CurrentExtremePoint(neighbour.index_perm);
								
								neighbour.neighs = FindNeighbours(neighbour);
								
								truncfacet.ReplaceNode(neighbour, *nodeit);     // The i-th node is replaced by the neighbour (saving the order)
								truncfacet.basic_perm = neighbour.index_perm;

								vector<int> protoperm(this->num - this->dim);
								vector<int> circleperm(this->dim);
								ExtremePoint protonode = neighbour;

								boost::dynamic_bitset<unsigned int> ggmap = neighbour.GetSmallPartialMap();

								int suchecircle = 0;
								int suchel = 0;
								int sucher = protoperm.size() - 1;
								tmpnodeit = truncfacet.nodes.begin();
								for(int gg = 0; gg < this->num; gg++)
								{
									if(tmpnodeit != truncfacet.nodes.end() && tmpnodeit->index_perm[this->border_index] == gg)
									{
										circleperm[suchecircle] = gg;
										tmpnodeit++;
										suchecircle++;
									}
									else
									{
										if( ggmap[gg] == 1)
										{
											protoperm[suchel] = gg;
											suchel++;
										}
										else
										{
											protoperm[sucher] = gg;
											sucher--;
										}
									}
								}

								int ins_place = suchel;

								// Checking in the hash-table : If current facet is new 
								if( ! this->ContainsFacet(truncfacet) )
								{
									// Mark facet by putting its hash-code to the hash-table
									this->MarkByHash(truncfacet);

									Facet extfacet = truncfacet;
									extfacet.nodes.clear();

									// d(d-1) vertices
									for(int uu = 0; uu < this->dim; uu++)
										for(int pp = 0; pp < this->dim; pp++)
										{
											vector<int> tcperm = circleperm;
											if(truncfacet.zone == 1)
											{
												if(pp != 1)
												{
													tcperm[1] = circleperm[uu];
													tcperm[uu] = circleperm[1];

													int temp1 = tcperm[0];
													tcperm[0] = tcperm[pp];
													tcperm[pp] = temp1;

													vector<int> truncindperm = protoperm;
													truncindperm.insert(truncindperm.begin() + ins_place, tcperm.begin(), tcperm.end());

													neighbour.index_perm = truncindperm;
													neighbour.coord = CurrentExtremePoint(neighbour.index_perm);
													extfacet.nodes.push_back(neighbour);
												}
											}
											else
											{
												if(pp != this->dim - 2)
												{
													tcperm[this->dim - 2] = circleperm[uu];
													tcperm[uu] = circleperm[this->dim - 2];

													int temp1 = tcperm[this->dim - 1];
													tcperm[this->dim - 1] = tcperm[pp];
													tcperm[pp] = temp1;

													vector<int> truncindperm = protoperm;
													truncindperm.insert(truncindperm.begin() + ins_place, tcperm.begin(), tcperm.end());

													neighbour.index_perm = truncindperm;
													neighbour.coord = CurrentExtremePoint(neighbour.index_perm);
													extfacet.nodes.push_back(neighbour);
												}
											}


										}


									// Push new facet to the algorithm queue
									this->trimmed_region.push_back(extfacet);
								}
							}

							// !!! End Stub !!!

						}
						
						else
						{
						
							// We've found a true neighbour - point t
							// !!! Warning !!! If currfacet is a cut, then there are more than 2 1s in common_neigh -> facets can be lost

							// Find index of the point t in any point of the current ridge
							if(nodeit != currfacet.nodes.begin())
							{
								tmpnodeit = currfacet.nodes.begin();
							}
							else
							{
								tmpnodeit = --(currfacet.nodes.end());
							}

							int r   = 0;
							while(tmpnodeit->index_perm[r] != t)
							{
								r++;
							}

							ExtremePoint neighbour = *tmpnodeit;

							// Swapping elements in neighbour
							int temp                                 = neighbour.index_perm[r];
							neighbour.index_perm[r]                  = neighbour.index_perm[this->border_index];
							neighbour.index_perm[this->border_index] = temp;

							// Compute coordinates of the new extr point
							neighbour.coord = CurrentExtremePoint(neighbour.index_perm);
							
							// Add neighbour to the current ridge and obtain new facet
							Facet partfacet = currfacet;              // Current ridge
							
							neighbour.neighs = FindNeighbours(neighbour);
							
							partfacet.ReplaceNode(neighbour, *nodeit);     // The i-th node is replaced by the neighbour (saving the order)
							partfacet.basic_perm = neighbour.index_perm;

							// Checking in the hash-table : If current facet is new 
							if( ! this->ContainsFacet(partfacet) )
							{
								// Mark facet by putting its hash-code to the hash-table
								this->MarkByHash(partfacet);

								// Push new facet to the algorithm queue
								this->queue.push_back(partfacet);
							}

						}

					}

				}

				// *** Create first facets from another blocks ***
				{

					// Changing zone for creating first facets
					int new_zone = 1 - currfacet.zone;

					int far_start;
					int far_stop;

					// If neighs are being searched in the first part of the permutation
					if(new_zone == 0)
					{
						far_start = 0;
						far_stop  = this->border_index;
					}
					// ... in the second part of ...
					else
					{
						far_start = this->border_index;
						far_stop  = this->num-1;
					}
			
					// Now we create all neighbours from new_zone of the nodes

					list<ExtremePoint>::iterator tmpnodeit;
					for(tmpnodeit = currfacet.nodes.begin(); tmpnodeit != currfacet.nodes.end(); tmpnodeit++)
					{

						Facet farfacet(this->dim);

						for(int k = far_start; k <= far_stop; k++)
						{
							// Defining a far neighbour and creating a first facet basing on it

							if(tmpnodeit->neighs[k])
							{
								ExtremePoint farneighbour = *tmpnodeit;

								// Swapping elements in farneighbour
								int temp                                    = farneighbour.index_perm[k];
								farneighbour.index_perm[k]                  = farneighbour.index_perm[this->border_index];
								farneighbour.index_perm[this->border_index] = temp;

								// Compute coordinates of the new extr point
								farneighbour.coord = CurrentExtremePoint(farneighbour.index_perm);										
								
								farneighbour.neighs = FindNeighbours(farneighbour);
								
								// Create a new first facet
								farfacet = FindFirstFacetInBlock(farneighbour, 1 - new_zone);              

								// Checking in the hash-table : If this first facet is new 
								if(!(farfacet.nodes.empty()) && ! this->ContainsFacet(farfacet) )
								{
									// Mark facet by putting its hash-code to the hash-table
									this->MarkByHash(farfacet);

									// Push new facet to the algorithm queue
									this->queue.push_back(farfacet);
								}
							}
						}

						farfacet = FindFirstFacetInBlock(*tmpnodeit, new_zone);

						// Checking in the hash-table : If this first facet is new 
						if(!(farfacet.nodes.empty()) && ! this->ContainsFacet(farfacet) )
						{
							// Mark facet by putting its hash-code to the hash-table
							this->MarkByHash(farfacet);

							// Push new facet to the algorithm queue
							this->queue.push_back(farfacet);
						}
					}
				}
			}
			
			// If currfacet is TRUNCATED
			else
			{
				
			}
			
			// Compute facet's coordinates and transfer it to the result list

			this->queue.pop_front();

			this->trimmed_region.push_back(currfacet);


		}while(!queue.empty());


		return this->trimmed_region;
	}

private:

	// Obtaining $d-1$ independent vectors from a set of $d$ points (in general position)
	vector<vector<float> > GenIndepVectors(list<int> _def_set)
	{
		vector<vector<float> > ind_bundle(this->dim);

		list<int>::iterator iter;
		list<int>::iterator pivot_iter = _def_set.begin();

		int j = 0;
		for(iter = _def_set.begin(); iter != _def_set.end(); iter++)
		{
			if(iter != pivot_iter)
			{
				Vector vec(this->dim);
				vec.FromPointToPoint(perm.points[*pivot_iter], perm.points[*iter]);
				ind_bundle[j] = vec.coord;
				j++;
			}
		}

		return ind_bundle;
	}

	// Generate combination (curr_comb.count num) by increasing pos element index
	vector<int> NextCombination(vector<int> curr_comb, int n)
	{
		if(curr_comb.size() == 0) return curr_comb;

		vector<int> new_comb = curr_comb;
		int c = new_comb.size();
		
		int pos = c-1;
		bool found = false;
		do
		{
			if(new_comb[pos] < n-c+pos)
			{
				new_comb[pos]++;
				found = true;
				for(int j = pos+1; j < c; j++)
				{
					new_comb[j] = new_comb[pos] + j - pos;
				}
			}

			pos--;
		}while( pos >= 0 && found == false);

		if(!found)
		{
			new_comb[0] = -1;
		}

		return new_comb;
	}

protected:

	
	// Generate set of points
	list<int> NextSet(list<int> curr_comb, int n)
	{
		list<int>::iterator lit;
		lit = curr_comb.begin();

		vector<int> vcomb(dim);

		for(int u = 0; u < dim; u++)
		{
			vcomb[u] = *lit;
			lit++;
		}

		vcomb = NextCombination(vcomb, n);
		
		lit = curr_comb.begin();
		for(int u = 0; u < dim; u++)
		{
			*lit = vcomb[u];
			lit++;
		}

		return curr_comb;
	}


	void MarkByHashH(Facet _new_fac)
	{
		// Warning: not optimal insertion algorithm!
		list<boost::dynamic_bitset<unsigned int> >::iterator hashit;

		bool inserted = false;
		
		boost::dynamic_bitset<unsigned int> curr_hash = _new_fac.GetHashMapH(this->num);

		if(this->hash_table.size() == 0)
		{
			this->hash_table.push_back(curr_hash);
		}
		else
		{
			for(hashit = this->hash_table.begin(); hashit != this->hash_table.end(); hashit++)
			{
				// Searching a place for insertion with saving the order
				if(ProcessAg::HashOrder(*hashit, curr_hash))
				{
					this->hash_table.insert(hashit, curr_hash);
					inserted = true;
					//::StopCumulTime();
					return;
				}
			}
			
			if(!inserted)
			{
				this->hash_table.insert(hash_table.end(), curr_hash);
			}
		}
	}

	// Creates a matrix (d x d) with first K rows - defining vectors of _fac
	vector<vector<float> > GetDefVectors_kFace(kFace _fac)
	{
		vector<vector<float> > def_vecs(this->dim);
		int plz = 0;

		list<int>::iterator itset, itcard;
		for(itset = _fac.anchors.begin(), itcard = _fac.cardinals.begin();
			itset != _fac.anchors.end(); 
			itset++, itcard++)
		{
			for(int y = 0; y < *itcard - 1; y++)
			{
				Vector dvec(this->dim);
				dvec.FromPointToPoint(perm.points[_fac.index_perm[*itset + y]], perm.points[_fac.index_perm[*itset + y + 1]]);
				
				def_vecs[plz] = dvec.coord;
				plz++;
			}
		}
	
		return def_vecs;
	}

	// Returns whether current facet has already been included into the trimmed region 
	bool ContainsFacetH(Facet _currfacet)
	{
		bool fres = binary_search(this->hash_table.begin(), this->hash_table.end(), _currfacet.GetHashMapH(this->num), &ProcessAg::ReversedHashOrder);
		return fres;
	}

	list<DataVector> GetCriticalVectors(Facet _currfac)
	{
		// !!! Attention: _currfac.anchors should modified (with removed vector)

		// The first step - forming associating clusters

		list<int>     assoc_clust;            // Association clusters starting points in the pi()
		list<bool>    clustyp;                // Types of the association clusters (true - Type I; false - Type II)

		list<int>::iterator iti, itic;
		list<int>::iterator itass;
		list<bool>::iterator itasst;

		int  afterI = 0;
		for(iti = _currfac.anchors.begin(), itic = _currfac.cardinals.begin(); iti != _currfac.anchors.end(); iti++, itic++)
		{
			// Gathering all type II clusters before current type I
			int typeIoccur = *iti;
			if(typeIoccur > afterI)
			{
				assoc_clust.push_back(afterI);
				clustyp.push_back( false );

				for(int i = afterI+1; i < typeIoccur; i++)
				{
					if( this->weight[i] > this->weight[i-1] )
					{
						assoc_clust.push_back(i);
						clustyp.push_back( false );
					}
				}

			}

			// Current type I cluster
			assoc_clust.push_back(typeIoccur);
			clustyp.push_back( true );

			afterI = typeIoccur + *itic;
		}

		// Gathering all type II clusters after the last type I
		if(afterI < this->num)
		{
			assoc_clust.push_back(afterI);
			clustyp.push_back( false );
		}
		for(int i = afterI + 1; i < this->num; i++)
		{
			if( this->weight[i] > this->weight[i-1] )
			{
				assoc_clust.push_back(i);
				clustyp.push_back( false );
			}
		}

		// The second step - Cartesian multiplication

		list<DataVector> found_vecs;   // ${\cal V}_{F*}$

		// Forming "associated sets" $X^j_{assoc}$
		list<list<int> > ass_sets;
		for(itass = assoc_clust.begin(), itasst = clustyp.begin(); itass != assoc_clust.end(); itass++, itasst++)
		{
			list<int> X_j_assoc; 
			
			if(*itasst)   // If type I
			{
				X_j_assoc.push_back(*itass);
			}
			else          // If type II
			{
				list<int>::iterator itassn = itass;
				itassn++; 

				if(itassn != assoc_clust.end())
				{
					for(int jj = *itass; jj < *itassn; jj++)
					{
						X_j_assoc.push_back(jj);
					}
				}
				else
				{
					for(int jj = *itass; jj < this->num; jj++)
					{
						X_j_assoc.push_back(jj);
					}
				}
			}

			ass_sets.push_back(X_j_assoc); 
		}

		list<list<int> >::iterator itx;
		for(itx = ass_sets.begin(), itasst = clustyp.begin(); itx != ass_sets.end(); itx++, itasst++)
		{
			list<list<int> >::iterator itxn    = itx;
			list<bool>::iterator      itasstn = itasst;
			itxn++, itasstn++;

			// *itx = $X^j_assoc$, *itxn = $X^{j+1}_assoc$, *itasst and *itasstn - corresponding cluster types

			// Cartesian multiplication followed by the union
			if(itxn != ass_sets.end()/* && !(*itasst && *itasstn)*/)
			{
				// Cartesian
				list<Vector> part_found;

				list<int> X_j_assoc  = *itx;
				list<int> X_j1_assoc = *itxn;

				for(iti = X_j_assoc.begin(); iti != X_j_assoc.end(); iti++)
				{
					for(itic = X_j1_assoc.begin(); itic != X_j1_assoc.end(); itic++)
					{
						DataVector crit_vec(this->dim, *iti, *itic);

						crit_vec.start_typeI = *itasst;
						crit_vec.end_typeI   = *itasstn;

						crit_vec.FromPointToPoint(this->perm.points[_currfac.index_perm[*iti]], this->perm.points[_currfac.index_perm[*itic]]);
						
						// Union
						found_vecs.push_back(crit_vec);
					}
				}
			}
		}

		return found_vecs; 

	}



public:

	// To a given k-dimensional face find some orthogonal vector that is also orthogonal to its normal
	Vector FindRotDirection(kFace _facet)
	{
		// Defining some rotation vector 
		vector<vector<float> > def_vectors = GetDefVectors_kFace(_facet);

		// Augmenting the matrix:
		// 1. Orthogonal to normvec
		def_vectors[_facet.K()] = _facet.normalvec.coord;
		//def_vectors[j] = normvec.coord;
		// 2. Last coordinates = 0 except of the last = 1
		for(int k = _facet.K()+1; k < dim; k++)
		{
			vector<float> onevec(this->dim, 0);
			onevec[k] = 1;

			def_vectors[k] = onevec;
		}

		Vector rotvec(this->dim);

		vector<float> bf(this->dim, 0);
		bf[this->dim-1] = 1;
		rotvec.coord = _facet.SolveLinearSystem(def_vectors, bf);

		rotvec.Normalize();

		return rotvec;

	}

	// Find the neighboring (k+1)-face in a given direction
	kFace RemoveOneDoF(kFace _currface, Vector _rotvec)
	{
		kFace newface  = _currface;
		Vector normvec = _currface.normalvec;

		list<DataVector> crit_vectors = this->GetCriticalVectors(newface);

		list<DataVector>::iterator itcrit;
		float minplambda = 1000001;
		float minnlambda = 1;
		Vector deltavec(this->dim);

		DataVector minp, minn;

		// Searching for a minimal rotation on an angle < 90 degrees, if there is
		for(itcrit = crit_vectors.begin(); itcrit != crit_vectors.end(); itcrit++)
		{
			deltavec.coord = itcrit->coord;
			
			float lambda = - Vector::ScalarMultiply(deltavec, normvec) / Vector::ScalarMultiply(deltavec, _rotvec);

			if(lambda > 0.0000001f)
			{
				if (lambda < minplambda)
				{
					minplambda = lambda;
					minp       = *itcrit;
				}								
			}
		}

		// If the rotation should done on the angle > 90 degrees
		if( !minp.defined() )
		{
			for(itcrit = crit_vectors.begin(); itcrit != crit_vectors.end(); itcrit++)
			{
				deltavec.coord = itcrit->coord;
				
				float lambda = - Vector::ScalarMultiply(deltavec, normvec) / Vector::ScalarMultiply(deltavec, _rotvec);
				
				if(lambda <= -0.0000001f)
				{
					if (lambda < minnlambda)
					{
						minnlambda = lambda;
						minn       = *itcrit;
					}
				}
			}
		}

		DataVector neighn;      // Found neighbor

		if( minp.defined() )
		{
			neighn  = minp;

			Vector vecp = normvec;
			Vector deltap = _rotvec;
			deltap.Scale(minplambda);
			vecp.Add(deltap);

			newface.normalvec = vecp;
		}
		else
		{
			neighn = minn;

			Vector vecn = normvec;
			vecn.Reverse();
			Vector deltan = _rotvec;
			deltan.Scale(-minnlambda);
			vecn.Add(deltan);

			newface.normalvec = vecn;
		}


		// Normalizing the normalvec (normvec & _rotvec must have length=1)
		float tonormal;
		if( minp.defined() )
			tonormal = 1.0f / (sqrt(1.0f + minplambda * minplambda));
		else
			tonormal = 1.0f / (sqrt(1.0f + minnlambda * minnlambda));

		newface.normalvec.Scale(tonormal);

		// If 2 Type II clusters do associate, create a new ${\cal{A}}_q$ 

		// Modifying ${\cal X}_F$ by adding the found vector (to some existing or new ${\cal{A}}_q$)

		// Editing index_perm
		// Warning: let neighn - the found vector to be intoduced

		if(neighn.start > neighn.end)
		{
			swap(neighn.start, neighn.end);

			swap(neighn.start_typeI, neighn.end_typeI);
		}

		list<int>::iterator srcit1, srcit2; 
		
		// Type I - Type I
		if(neighn.start_typeI && neighn.end_typeI)
		{
			// newface.index_perm stays unchanged

			// Locating corresponding ${\cal{A}}$s
			for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
				srcit1 != newface.anchors.end(); 
				srcit1++, srcit2++)
			{
				if(*srcit1 > neighn.start) break;
			}
			
			newface.anchors.erase(srcit1); 
			
			int added = *srcit2;
			srcit2--;
			*srcit2 += added;
			srcit2++;
			newface.cardinals.erase(srcit2);
			
		}
		// Type II - Type II
		else if(!neighn.start_typeI && !neighn.end_typeI)
		{
			// Shifting start point to the right border of the type II cluster
			int qq = neighn.start;
			while(this->weight[qq+1] == this->weight[qq])
			{
				qq++;
			}

			if(qq > neighn.start)
			{
				swap(newface.index_perm[qq], newface.index_perm[neighn.start]);
				neighn.start = qq;
			}

			// Shifting end point to the left border of the type II cluster
			qq = neighn.end;
			while(this->weight[qq-1] == this->weight[qq])
			{
				qq--;
			}

			if(qq < neighn.end)
			{
				swap(newface.index_perm[qq], newface.index_perm[neighn.end]);
				neighn.end = qq;
			}

			// newface.index_perm is ready

			// Locating place for a new ${\cal{A}}$
			for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
				srcit1 != newface.anchors.end(); 
				srcit1++, srcit2++)
			{
				if(*srcit1 > neighn.start) break;
			}

			newface.anchors.insert(srcit1, neighn.start);
			newface.cardinals.insert(srcit2, 2);
		}
		else
		{
			// Type I - Type II
			if(neighn.start_typeI)
			{
				//// Shifting end point to the left border of the type II cluster
				//int qq = neighn.end;
				//while(this->weight[qq-1] == this->weight[qq])
				//{
				//	qq--;
				//}

				// newface.index_perm is ready

				// Locating the ${\cal{A}}$
				for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
					srcit1 != newface.anchors.end(); 
					srcit1++, srcit2++)
				{
					if(*srcit1 > neighn.start) break;
				}

				srcit2--;
				*srcit2 += 1;
				srcit1--;

				if(*srcit1 + *srcit2 - 1 < neighn.end)
				{
					swap(newface.index_perm[*srcit1 + *srcit2 - 1], newface.index_perm[neighn.end]);
					neighn.end = *srcit1 + *srcit2;
				}
			}
			// Type II - Type I
			else
			{
				//// Shifting start point to the right border of the type II cluster
				//int qq = neighn.start;
				//while(this->weight[qq+1] == this->weight[qq])
				//{
				//	qq++;
				//}

				// newface.index_perm is ready

				// Locating the ${\cal{A}}$
				for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
					srcit1 != newface.anchors.end(); 
					srcit1++, srcit2++)
				{
					if(*srcit1 > neighn.start) break;
				}

				*srcit1 -= 1;
				*srcit2 += 1;

				if(*srcit1 > neighn.start)
				{
					swap(newface.index_perm[*srcit1], newface.index_perm[neighn.start]);
					neighn.start = *srcit1;
				}
			}
		}

		newface.CatchOneDoF();
		return newface;
	}

public:

	void ComputeHyperplanes()
	{
		// ************************* Smart first facet ****************************


		// Create an "empty" first facet
		kFace ffac(this->dim, 0);

		// Initial $\pi(i)$ - trivial
		vector<int>  initial_ind(num);
		for(int fu = 0; fu < num; fu++)
			initial_ind[fu] = fu;

		ffac.index_perm = initial_ind;

		// Initial support vector determining the direction for the ffacet search
		Vector initvec(this->dim);

		// Getting the first permutation - only once for the first facet!
		this->perm.Support(initvec);

		// Now, only anchors and cardinals must be constructed (also normalvec)

		ffac.normalvec = initvec;
		ffac.normalvec.Normalize();

		// Grasping d-1 degrees of freedom
		for(int j = 0; j < this->dim - 1; j++)
		{
			// Create the direction to shift the normal vector
			Vector rotvec = this->FindRotDirection(ffac);

			// Moving normvec by means of rotvec, which is orthogonal to normvec

			ffac = this->RemoveOneDoF(ffac, rotvec);
	
		}
		
		//ffac.normalvec.Normalize();
		Facet found_facet(ffac);
		found_facet.doubled   = false;

		Point ffnode(this->dim);
		ffnode.coord = CurrentExtremePoint(found_facet.index_perm);
		found_facet.CalculateAbsoluteMemberH(ffnode);
		
		// Push new facet to the algorithm queue
		this->queue.push_front(found_facet);

		// ************************************************************************


	
		/*
		// ************************** Create the first facet **********************

		list<int>    fdefset(dim);
		int fu = 0;
		for(list<int>::iterator fit = fdefset.begin(); fit != fdefset.end(); fit++)
		{
			*fit = fu;

			fu++;
		}

		bool ffound = false;

		// Create an "empty" first facet
		Facet ffacet(this->dim);

		do
		{
			// Try combination
			
			//ffacet.def_set = fdefset;

			// Constructing normal-vector
			ffacet.CalculateEquationH(GenIndepVectors(fdefset), b);

			Permutation fperm = this->perm;
			fperm.Support(ffacet.normalvec);
			
			// Checking the conditions (3)
			Vector deltavec(this->dim);

			for(int j = 0; j < this->num - 1; j++)
			{
				// Searching for group of defining points
				deltavec.FromPointToPoint(fperm.points[j], fperm.points[j+1]);
				if( fabs(Vector::ScalarMultiply(ffacet.normalvec, deltavec)) < 0.001f )
				{
					if(this->weight[j] < this->weight[j+dim-1])
					{
						ffound = true;

						this->perm = fperm;
						this->p    = ffacet.normalvec;

						ffacet.anchors.clear();
						ffacet.cardinals.clear();
						
						ffacet.anchors.push_back(j);
						ffacet.cardinals.push_back(dim);
					}

					break;
				}
			}

			fdefset = NextSet(fdefset, this->num);

		}while(!ffound);

		// Final processing
		{
			ffacet.doubled   = false;

			// Initial $\pi(i)$ - trivial
			vector<int>  initial_ind(num);
			for(fu = 0; fu < num; fu++)
				initial_ind[fu] = fu;

			ffacet.index_perm = initial_ind;

			ffacet.normalvec.Normalize();

			Point fnode(this->dim);
			fnode.coord = CurrentExtremePoint(ffacet.index_perm);
			ffacet.CalculateAbsoluteMemberH(fnode);
			
			// Push new facet to the algorithm queue
			this->queue.push_front(ffacet);

		}
		*/

		::RecordTime();
		
		// ************************************************************************

		// ************************** Create all the facets ***********************
				
		do
		{
			// Take the upper facet $F$
			Facet currfacet = queue.front();

			this->queue.pop_front();

			//if(currfacet.doubled) 
			//{
			//	queue.pop_front();
			//	if(queue.empty()) break;
			//	continue;
			//}

			// Removing one vector thus enabling rotating

			// Iterating through ${\cal{A}}_l \in {\cal X}_F$
			list<int>::iterator iterset, itercard;
			for(iterset = currfacet.anchors.begin(), itercard = currfacet.cardinals.begin(); 
				iterset != currfacet.anchors.end(); 
				iterset++, itercard++)
			{
				// "Cracking" each ${\cal{A}}_l$ into 2 sets. There are $2^{\|{\cal{A}}\|}-2$

				list<boost::dynamic_bitset<unsigned int> > crackmasks;

				for(int ccc = 1; ccc < pow(2.0, *itercard) - 1.5; ccc++ )
				{
					boost::dynamic_bitset<unsigned int> crackmask(*itercard, ccc);

					crackmasks.push_back(crackmask);
				}

				// Filtering crackmaps according to the condition (3)

				list<boost::dynamic_bitset<unsigned int> >::iterator itcrack;

				//for(int k = *iterset; k < *iterset+*itercard; k++)
				for(itcrack = crackmasks.begin(); itcrack != crackmasks.end(); itcrack++)
				{
					int bulk = itcrack->count();

					// Moving in the permutation according to the cracking map
					if( (this->weight[*iterset+*itercard-1] > this->weight[*iterset+bulk]) && ( this->weight[*iterset+bulk-1] > this->weight[*iterset]) ||
						(bulk == 1)           && (this->weight[*iterset+*itercard-1] > this->weight[*iterset+bulk]) ||
						(bulk == *itercard-1) && (this->weight[*iterset+bulk-1]      > this->weight[*iterset]) ||
						*itercard == 2)
					{
						// Splitting the ${\cal{A}}_l$ according to the cracking map

						//::ResumeCumulTime();
						// Be cautious: currfacet must be recovered
						currfacet.anchors.insert(iterset, *iterset);
						int origset = *iterset;
						*iterset  = *iterset + bulk;
						list<int>::iterator itsetnew  = iterset;
						iterset--;

						currfacet.cardinals.insert(itercard, bulk);
						int origcard = *itercard;
						*itercard = *itercard - bulk;
						list<int>::iterator itcardnew = itercard;
						itercard--;

						// When removing only 1 point, or eliminating ${\cal{A}}_l$ containing only 2 points
						if(bulk == 1)
						{
							list<int>::iterator empset = iterset;
							list<int>::iterator empcard = itercard;
							iterset++; iterset++;      // 2-step move
							itercard++; itercard++;      // 2-step move
							currfacet.anchors.erase(empset);
							currfacet.cardinals.erase(empcard);
						}

						if(bulk == itcrack->size()-1)
						{
							currfacet.anchors.erase(itsetnew);
							currfacet.cardinals.erase(itcardnew);
						}
						
						kFace newridge = currfacet;
						newridge.RelaxOneDoF();
						
						// Modifying index_perm
						int ww1 = origset, ww2 = origset + origcard - 1;
						for(int ww = 0; ww < itcrack->size(); ww++)
						{
							if((*itcrack)[ww])
							{
								newridge.index_perm[ww1] = currfacet.index_perm[origset + ww];
								ww1++;
							}
							else
							{
								newridge.index_perm[ww2] = currfacet.index_perm[origset + ww];
								ww2--;
							}
						}

						// Saving removed vector (one of the possible)
						//newridge.old_vector.FromPointToPoint(this->perm.points[newridge.index_perm[origset]], this->perm.points[newridge.index_perm[origset + bulk]]);
						
						// Recovering currfacet
						if(bulk < itcrack->size()-1)
						{
							currfacet.anchors.erase(itsetnew);
							currfacet.cardinals.erase(itcardnew);
						}

						if(bulk == 1)
						{
							currfacet.anchors.insert(iterset, origset);
							currfacet.cardinals.insert(itercard, origcard);
							iterset--;
							itercard--;
						}

						*iterset  = origset;
						*itercard = origcard;

						// currfacet recovered

						// Either only marking or processing the currfacet
						if(!currfacet.marked)
						{
							this->ohash_table.MarkRidge(newridge);
							continue;
						}

						// If the ridge has been double blocked - don't process it 
						if(this->ohash_table.CheckRidge(newridge))
							continue;

						// Taking a vector to detect direction of an infinitesimal rotation (either clockwise or counterclockwise)
						Vector infinitesimalvec(this->dim);
						infinitesimalvec.FromPointToPoint(perm.points[newridge.index_perm[origset]], perm.points[newridge.index_perm[origset + origcard - 1]]);

						Vector rotvec = this->FindRotDirection(newridge);

						if( Vector::ScalarMultiply(infinitesimalvec, rotvec) <  -0.00001f)
							rotvec.Reverse();

						// Now rotvec is orthogonal to currfacet.normalvec. <rotvec;currfacet.normalvec> is a B2 basis
						// Direction of rotation (clcws. or counterclcws.) is shown by infinitesimalvec 

						/*
						// Calculating unique number for the newridge in the combinatorially defined class of ridges
						Vector univec(this->dim);
						univec.coord = CurrentExtremePoint(newridge.index_perm);

						Vector proj1(this->dim), proj2(this->dim), proj(this->dim);
						proj1 = currfacet.normalvec; proj1.Scale(Vector::ScalarMultiply(univec, currfacet.normalvec));
						proj2 = rotvec;				 proj2.Scale(Vector::ScalarMultiply(univec, rotvec));
						proj = proj1; proj.Add(proj2);

						float unifl  = proj.coord[0] / proj.GetLength();
						// Managing 'float comparison' problem
						int   uniint = (int)(unifl * 111512);
						*/

						//::StopCumulTime();


						Facet neigh_facet(this->RemoveOneDoF(newridge, rotvec));
						neigh_facet.marked = false;

						this->queue.push_front(neigh_facet);




					} // end processing a feasible crackmask
				}  // end of loop for a current crackmask

		}

			// Transfer currfacet from queue to zonoid region

			if(!currfacet.marked)
			{
				currfacet.marked = true;
				this->queue.push_back(currfacet);
			}
			else
			{

				// Calculating absolute member to accomplish identification
				Point idnode(this->dim);
				idnode.coord = CurrentExtremePoint(currfacet.index_perm);
				currfacet.CalculateAbsoluteMemberH(idnode);

				//::RecordTime();
				this->trimmed_region.push_back(currfacet);
			}


		}while(!queue.empty());

		//return this->trimmed_region;
	}

private:


	// Returns a set of splittings of the _sublist into a selected combination and the rest 
	template<class ArbClass> list<vector<vector<ArbClass> > > SelectCombins(vector<ArbClass> _vec, int _num)
	{
		int veclen = _vec.size(); 
		int numofcombs = Comb(veclen, _num);

		list<vector<vector<ArbClass> > > reslist;

		vector<int> comb(_num), contracomb(veclen - _num);
		for(int i=0; i < _num; i++)
		{
			comb[i]       = i;
		}
		for(int i=_num; i < veclen; i++)
		{
			contracomb[i-_num] = i;
		}

		vector<vector<ArbClass> > splitt(2);
		vector<ArbClass> empt1(_num);          splitt[0] = empt1;
		vector<ArbClass> empt2(veclen - _num); splitt[1] = empt2;

		for(int w = 0; w < _num; w++)
			(splitt[0])[w] = _vec[comb[w]]; 
		for(int w = 0; w < veclen - _num; w++)
			(splitt[1])[w] = _vec[contracomb[w]];
		//splitt[0] = comb; splitt[1] = contracomb;

		reslist.push_back(splitt);
		
		for(int q=0; q < numofcombs-1; q++)
		{
			comb = NextCombination(comb, veclen);
			
			int hund = 0;
			for(int i=0; i < veclen; i++)
			{
				if(hund < _num && comb[hund] == i )
					hund++;
				else contracomb[i-hund] = i; 
			} 
		
			for(int w = 0; w < _num; w++)
				(splitt[0])[w] = _vec[comb[w]]; 
			for(int w = 0; w < veclen - _num; w++)
				(splitt[1])[w] = _vec[contracomb[w]];
			//splitt[0] = comb; splitt[1] = contracomb;
		
			reslist.push_back(splitt);
		}

		return reslist;
	}

public:
	
	void CalculateAllVertices()
	{
		list<Facet>::iterator facetit;

		// For each facet
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{
			list<int>::iterator iterit, itercard;

			// Combinations within ${\cal A}_l$ are independent of ones in ${\cal A}_k, k \ne l$
			vector<list<vector<vector<int> > > > als;
			for(iterit = facetit->anchors.begin(), itercard = facetit->cardinals.begin();
				iterit != facetit->anchors.end();
				iterit++, itercard++)
			{
				// Create a list from $A_l$ and weights homogeneity map
				vector<int> al, alordered;
				list<int> alwmap;
				int stopal = *iterit+*itercard;
					
				int lastjump = *iterit;
				for(int i = *iterit; i < stopal; i++)
				{
					//al.push_back(facetit->index_perm[i]);
						
					if( (i > *iterit) && (weight[i] > weight[i-1]) )
					{
						alwmap.push_back(i - lastjump);
						lastjump = i; 
					}
				}
				alwmap.push_back(stopal - lastjump);

				al.assign(facetit->index_perm.begin() + *iterit, facetit->index_perm.begin() + stopal);

				list<vector<vector<int> > > alcombs; // List of 2-element vectors - splittings of A_l into permuted and unpermuted parts
				vector<vector<int> > defcomb(2);
				defcomb[0] = alordered; defcomb[1] = al;
				alcombs.push_back(defcomb);

				list<vector<vector<int> > >::iterator itspl,itpartspl;
				do
				{
					list<vector<vector<int> > > nextalcombs;

					for(itspl = alcombs.begin(); itspl != alcombs.end(); itspl++)
					{
						list<vector<vector<int> > > partcombs = SelectCombins((*itspl)[1], alwmap.front());
						
						// Augmenting partial combination with an already ordered part
						for(itpartspl = partcombs.begin(); itpartspl != partcombs.end(); itpartspl++)
						{
							(*itpartspl)[0].insert((*itpartspl)[0].begin(), (*itspl)[0].begin(), (*itspl)[0].end());
						}

						nextalcombs.splice(nextalcombs.begin(), partcombs);
					}

					alcombs = nextalcombs;

					alwmap.pop_front();
				}while(!alwmap.empty());

				als.push_back(alcombs);
			}

			// Generating final permutations (combining independent variants from all $A_l$)
			list<vector<int> > ready_perms;
			vector<list<vector<vector<int> > > >::iterator iterals;

			ready_perms.push_back(facetit->index_perm);

			for(iterit = facetit->anchors.begin(), itercard = facetit->cardinals.begin(), iterals = als.begin();
				iterit != facetit->anchors.end();
				iterit++, itercard++, iterals++)
			{
				list<vector<int> > step_ready_perms;
				int stopal = *iterit+*itercard;

				for(list<vector<int> >::iterator  itsub = ready_perms.begin(); itsub != ready_perms.end(); itsub++)
				{
					for(list<vector<vector<int> > >::iterator tmpal = iterals->begin(); tmpal != iterals->end(); tmpal++)
					{
						vector<int> augmperm = *itsub;
						
						for(int j = *iterit; j < stopal; j++)
						{
							augmperm[j] = ((*tmpal)[0])[j-*iterit];
						}

						step_ready_perms.push_back(augmperm);
					}
				}

				ready_perms = step_ready_perms;
			}

			// Final genearing of vertices

			facetit->nodes.clear();

			for(list<vector<int> >::iterator  itsub = ready_perms.begin(); itsub != ready_perms.end(); itsub++)
			{
				ExtremePoint newnode(this->dim);

				newnode.coord = CurrentExtremePoint(*itsub);
				facetit->nodes.push_back(newnode);
			}
			
		}
	}
	
	void CalculateAllVertices3d()
	{
		if(dim == 3)
		{
			list<Facet>::iterator facetit;

			// For each facet
			for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
			{
				// If there are 2 Al's
				if(facetit->anchors.size() == 2)
				{
					// Faces with 4 vertices are generated
					facetit->truncated = true;

					vector<int> piperm = facetit->index_perm;

					ExtremePoint idnode(this->dim);

					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

					swap(piperm[facetit->anchors.front()], piperm[facetit->anchors.front()+1]);
					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

					swap(piperm[facetit->anchors.back()], piperm[facetit->anchors.back()+1]);
					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

					swap(piperm[facetit->anchors.front()], piperm[facetit->anchors.front()+1]);
					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

				}
				// If there is only 1 Al
				else
				{
					int unianc = facetit->anchors.front();

					// If not all weights are different - 3 vertices
					if( this->weight[unianc] == this->weight[unianc + 1])
					{
						facetit->truncated = false;

						vector<int> piperm = facetit->index_perm;

						ExtremePoint idnode(this->dim);

						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc + 1], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
					}
					else if( this->weight[unianc+1] == this->weight[unianc + 2] )
					{
						facetit->truncated = false;

						vector<int> piperm = facetit->index_perm;

						ExtremePoint idnode(this->dim);

						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
					}
					// Else - 6 vertices
					else
					{
						facetit->truncated = true;

						vector<int> piperm = facetit->index_perm;

						ExtremePoint idnode(this->dim);

						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc + 1], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc + 1], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
					}
					
					/*list<ExtremePoint>::iterator iterex;
					for(iterex = facetit->nodes.begin(); iterex != facetit->nodes.end(); iterex++)
					{
						Vector vec(this->dim);
						vec.coord = iterex->coord;
						float chdelta = fabs(Vector::ScalarMultiply(vec, facetit->normalvec) + facetit->abs_member);
						if(chdelta > 0.0001f)
						{
							vec.Reverse(); // ERROR
						}
					}*/
				}
			}
		}

	}


	list<Facet> ReadyTR()
	{
		return this->trimmed_region;
	}
};


// ******************************* Class for receiving original data cloud **************************

class InputAg : public Agent
{
public:

	int   dim  ;                  // Dimension of the cloud
	vector<Point> cloud;
	Point centroid;

private:

	string WMTD_type;
	float depth;                  // Parameter for computing trimmed region
	int   num  ;                  // Number of points
	vector<float> manweights;

public:

	InputAg()
	{
		dim   = 3;
		num   = 0;

		depth = 1;
	}

	// Receive cloud of points from the specified text file
	int Receive(const char* _source,const char* _dir)
	{
		try
		{
			char pathtofile[100]; strcpy(pathtofile, _dir);
			ifstream is(strcat(pathtofile,_source));     // Open file with data

			is >> WMTD_type;

			is >> depth;
			
			is >> dim;
			is >> num;

			// If the type of the trimmed region is given by manual specifying of the weight vector
			if(this->WMTD_type == "general")
			{
				vector<float> tmpwe(num);
				manweights = tmpwe;
				for(int i = 0; i < num; i++)
				{
					is >> manweights[i];

				}
			}

			// Reading coordinates to objects...
			for(int i=0; i<num; i++)
			{
				Point p(dim);
				for(int j=0; j<dim; j++)
				{
					is >> p.coord[j];
				}
				
				cloud.push_back(p);       // Adding a point to cloud

			}
			
			is.close();

			return 0;
		}
		catch(std::exception ex)
		{
			return -1;
		}
		
		return 0;
	}

	// Receive cloud of points from the specified text file
	int Receive(int _d, int _n, int _ind, char* _type, char* _dir)
	{
		cloud.clear();

		try
		{
			char filenm[200];
			sprintf(filenm, "Cloud_%d_%d_%d.dat", _d, _n, _ind);
			ifstream is(strcat(_dir,filenm));     // Open file with data

			//is >> WMTD_type;
			WMTD_type = _type;

			is >> depth;
			
			is >> dim;
			is >> num;

			// If the type of the trimmed region is given by manual specifying of the weight vector
			if(this->WMTD_type == "general")
			{
				vector<float> tmpwe(num);
				manweights = tmpwe;
				for(int i = 0; i < num; i++)
				{
					is >> manweights[i];

				}
			}

			// Reading coordinates to objects...
			for(int i=0; i<num; i++)
			{
				Point p(dim);
				for(int j=0; j<dim; j++)
				{
					is >> p.coord[j];
				}
				
				cloud.push_back(p);       // Adding a point to cloud

			}
			
			is.close();

			return 0;
		}
		catch(std::exception ex)
		{
			return -1;
		}
		
		return 0;
	}

	/*list<Facet>*/void StartComputing(ResultAg* _result_agent)
	{
		// Creating main processor for realizing algorithm 
		ProcessAg* processor = new ProcessAg(WMTD_type, depth, dim, num, cloud, manweights);

		this->cloud    = processor->perm.points;
		this->centroid = processor->initial_center;
		
		::StartTiming();
		//return processor->Compute();
		processor->ComputeHyperplanes();
		::RecordTime();

		if(this->dim == 3)
			processor->CalculateAllVertices3d();
		else
			processor->CalculateAllVertices();

		//return processor->ReadyTR();
		_result_agent->trimmed_region = processor->ReadyTR();

	}

};

// ************************************************************************************************
// *********************** Main Function **********************************************************

extern "C"{

void  ComputeWMTR(char** nameofsource, char** _wdir, int* _bound)			// Window Show State
{
	
	bool single_mode = true;

	ResultAg* result_agent = new ResultAg();
	InputAg*  input_agent  = new InputAg();
	
	if(single_mode)
	{
		int res = input_agent->Receive(*nameofsource, *_wdir);

		/*result_agent->trimmed_region = */input_agent->StartComputing(result_agent);
		result_agent->data_cloud     = input_agent->cloud;
		result_agent->coord_center   = input_agent->centroid;

		result_agent->PrintResultsHyperplanes(*_wdir, *_bound);
		result_agent->PrintResults(*_wdir, *_bound);
		result_agent->GLSceneToExport(*_wdir, *_bound);
	}
	else
	// Group mode
	{
		for(int dvar = 4; dvar <= 6; dvar ++)
			for(int nvar = 10/*dvar + 1*/; nvar <= 25; nvar += 5)
			{
				for(int tp = 0; tp < 4; tp++)
				{
					char* trtype;
					
					switch(tp)
					{
					case 0:
						{
							trtype = "ECH";
							break;
						}
					case 1:
						{
							trtype = "contECH";
							break;
						}
					case 2:
						{
							trtype = "geometrical";
							break;
						}
					case 3:
						{
							trtype = "zonoid";
							break;
						}
					}

					int indcount = 2;
					for(int ivar = 1; ivar <= indcount; ivar ++)
					{
						if(input_agent->Receive(dvar, nvar, ivar, trtype, *_wdir) == 0)
						{
							/*result_agent->trimmed_region = */input_agent->StartComputing(result_agent);
							result_agent->data_cloud     = input_agent->cloud;
							result_agent->coord_center   = input_agent->centroid;

							result_agent->PrintResultsHyperplanes(dvar, nvar, ivar, indcount, trtype);
							
							// Free memory
							result_agent->trimmed_region.clear();
						}
					}
				}
			}
	}


}

}
