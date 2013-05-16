#ifndef NN_AUX_H
#define NN_AUX_H

#include <cmath>
#include <stack>
#include <list>
#include <vector>
#include <algorithm>
#include <queue>

using namespace std;

// INFINITY is needed as an upper bound for every
// distance that might be encountered during search
// So be carefull : if the data set has a very unusual scaling,
// actual distances might exceed this value, which would lead
// to wrong results

#ifndef DBL_MAX
#define DBL_MAX 1E+38
#endif

#ifndef INFINITY
#define INFINITY DBL_MAX
#endif

class neighbor	// a "neighbor" is described by its index (into the set of points) and it's distance (to the query point)
{
	protected:
		long i; 	  // index of neighbor
		double d;	  // distance of neighbor to query point	  
	public:
		neighbor() {};
		neighbor(const long I, const double D) : i(I), d(D) {}; 
		inline long index() const { return i; };
		inline double dist() const { return d; };
		inline long& index() { return i; }
		inline double& dist() { return d; }
		inline bool operator< (const neighbor& x) const { return d < x.d; }            		
};

class neighborCompare : public binary_function<neighbor, neighbor, bool> {  
	public:
  		inline bool operator()(const neighbor& x, const neighbor& y) const { return x.dist() < y.dist(); } 
};

class SortedNeighborTable
{
	protected:
		long NNR;		// number of neighbors to be searched
  priority_queue<neighbor, vector<neighbor> , neighborCompare > pq;	
		double hd;		// cache highest distance
	
	public:
		SortedNeighborTable() : NNR(1), hd(INFINITY) {};
		SortedNeighborTable(const long nnr) : NNR(nnr), hd(INFINITY) {};
		~SortedNeighborTable() {};
			
		inline double highdist() const { return hd; };
		void insert(const neighbor& x);
		
		inline void init_search(const long nnr) { NNR = nnr; hd = INFINITY; }
		long finish_search(vector<neighbor>& v);
};

class cluster
{
	private:
#ifdef USE_OWN_CLUSTER_MEMORY_HANDLER	
		static cluster* headOfFreeList;	
#endif
	
	public:
		long center;		// index of center point for this cluster (points direcly into the point set 
		double Rmax;		// if Rmax <= 0 we have a terminal node (so we have to use fabs(Rmax) to get the true value for Rmax)							
		double g_min;	

		union {
			cluster* left;	// used in case of a non-terminal node
			long start;		// used in case of a terminal node
		};
		union {
			cluster* right;		// used in case of a non-terminal node
			long length;		// used in case of a terminal node
		};

		// this class-specific constants specifies how
	  	// many cluster objects fit into a big memory block;
	  	
		static long OLD_BLOCK_SIZE;
	  	static long BLOCK_SIZE;

		cluster() : center(0), Rmax(INFINITY), g_min(0), start(0), length(0)  {};
		cluster(const long c) : center(c), Rmax(INFINITY), g_min(0), start(0), length(0) {};
		cluster(const long s, const long l) : center(0), Rmax(INFINITY), g_min(0), start(s), length(l)  {};
		cluster(const long s, const long l, const long c) : center(c), Rmax(INFINITY), g_min(0), start(s), length(l)  {};		

		~cluster() {};
	
		inline int is_terminal() const { return (Rmax <= 0); }; 
		inline double R_max() const { return fabs(Rmax); };

#ifdef USE_OWN_CLUSTER_MEMORY_HANDLER
	  	static void* operator new(size_t size);
 		static void operator delete(void *deadObject, size_t size);
#endif		
};

// During k-nearest neighbor search, clusters are treated as searchitems
// These searchitems are inserted into a priority queue
class searchitem
{
	protected:	
		const cluster* c; 	// pointer to cluster object
		double d;			// distance from query point to the cluster's center
		double dbrother;	// distance from query point to brother cluster's center
		
		double dmin;		// miniumum distance from  query point to any point in cluster, accumulated through all tree levels
	 	double dmax;		// maximal distance from  query point to any point in cluster, accumulated through all tree levels	
	public:
		searchitem() {};
			
		inline searchitem(const cluster* C, const double D) 	
		: c(C), d(D), dbrother(INFINITY), dmin(D - C->R_max()), dmax(D+C->R_max()) {};	
			
		inline searchitem(const cluster* C, const double D, const double Dbrother, const searchitem& parent)
		 :  c(C), d(D), dbrother(Dbrother),	 
		 	dmin(std::max(std::max(0.0, 0.5*(D-Dbrother+c->g_min)), std::max(D - C->R_max(), parent.dmin))),
			dmax(std::min(parent.dmax, D + C->R_max())) {}; 	
	
		inline const cluster* clusterp() const { return c; };
		inline double dist() const { return d; };	
		inline double dist_brother() const { return dbrother; };	
		
		inline double d_min() const { return dmin; }; 		// accumulated value of d_min for this cluster
		inline double d_max() const { return dmax; }; 		// accumulated value of d_max for this cluster
};

class searchitemCompare : public binary_function<searchitem, searchitem, bool> {  
	public:		
		inline bool operator()(const searchitem& x, const searchitem& y) const {				
			if (x.d_min() == y.d_min()) 		
				return x.d_max() > y.d_max();	
			else
				return x.d_min() > y.d_min(); 			
		}; 
};

typedef cluster* cluster_pointer;							
typedef vector<cluster_pointer> cluster_pointer_vector; 	

#endif
