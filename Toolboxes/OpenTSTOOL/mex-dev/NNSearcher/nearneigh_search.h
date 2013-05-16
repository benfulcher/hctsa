#ifndef NEARNEIGH_SEARCH_H
#define NEARNEIGH_SEARCH_H

/* This file contains a class collection for nearest neighbor (NN) searching */
/* It is possible to search k (k=1..N) nearest neighbors to a given query point ("fixed-mass" approach) */
/* or to search all neighbors within a distance r (range) around the query point ("fixes size" approach) */
  
/* Several different algorithms for NN search are available. */
/* These algorithms are encapsulated in different classes. */
/* All algorithms should work with any kind of metric, not only the euclidian metric */

/* These algorithms are available : */

/* Brute : O(N*N) A class that implements a brute force approach. The advantage of */
/* brute force is that no preprocessing is required. If only a few query */
/* points are given and the number of points in the data set is small, brute force */
/* can be a good choice, but in general there are faster algorithms */

/* ATRIA : A class that implements an advanced triangle inequality algorithm */
/* During preprocessing, a search tree is constructed by dividing the set of data points */
/* in two (sub)clusters. Each cluster is than subdivided until a minimum number of points is reached */
/* During search, the triangle inequality is used to exclude cluster from further search */
/* ATRIA might be a good choice for unevenly distributed points in very high dimensional spaces. */

/* Christian Merkwirth, DPI Goettingen, 1998 - 2000 */

/* Here are some parameters that might need some tuning for your special case, */
/* but for general purposes these are OK */

/* Parameters for the ATRIA nearest neighbor search */
/* A cluster will not be further subdivided if it contains less */
/* than ATRIAMINPOINTS points */
/* A smaller value might accelerate actual search */
/* but increase pre-processing time */
/* Memory consumption will not change very much when choosing a smaller */
/* value for ATRIAMINPOINTS */

#define ATRIAMINPOINTS 64

// header files beloning to this package
#include "mextools/Utilities.h"
#include "nn_aux.h"

// base class for the nearest neihbor search which defines a common interface
template<class POINT_SET>
class nearneigh_searcher : protected My_Utilities	
{
	protected :
		int err;  					// error state, == 0 means OK, every other value is a failure
	
		const POINT_SET& points;
#ifdef MATLAB_MEX_FILE
		long Nused;
#else		
		const long Nused; 			// number of points of the data set actually used
#endif

#ifdef PROFILE
		long points_searched;
		long number_of_queries;
#endif	

		SortedNeighborTable table;	

		// test point number #index of points 
		
		template<class ForwardIterator>
		inline void test(const long index, ForwardIterator qp, const double thresh) {
#ifdef PARTIAL_SEARCH		
			const double d = points.distance(index, qp, thresh);
#else
			const double d = points.distance(index, qp);
#endif			
			if (d < thresh) 
				table.insert(neighbor(index,d));
#ifdef PROFILE
			points_searched++;
#endif
		} 
				
	public :
		typedef POINT_SET point_set;
	
		// prepare searching for a point set points
		// excl gives the number of samples from the end of points that should be omitted
		// from searching
		
		nearneigh_searcher(const POINT_SET& p, const long excl = 0);
		~nearneigh_searcher();
		
		// return the error state, = 0 means OK
		// an error will only occur when the object was not initialized correctly
		inline int geterr() const { return err; } 
		
		long Number_of_Points() const { return Nused; }
		const POINT_SET& get_point_set() const { return points; }
					
		double search_efficiency() const { 
#ifdef PROFILE
			return ((double)points_searched)/((double)Nused*number_of_queries);
#else
			return -1;		// search_efficiency not measured
#endif
		}			 		
};

// Brute Force search algorithm
template<class POINT_SET>
class  Brute : public nearneigh_searcher<POINT_SET>
{
     using nearneigh_searcher<POINT_SET>::Nused;
     using nearneigh_searcher<POINT_SET>::table;
     using nearneigh_searcher<POINT_SET>::points;
     using nearneigh_searcher<POINT_SET>::err;

	protected:
		template<class ForwardIterator>		
		void search(ForwardIterator query_point, const long first, const long last, const double epsilon);	
			
	public:
		Brute(const POINT_SET& p, const long excl = 0);
		~Brute() {};

		template<class ForwardIterator>
		long search_k_neighbors(vector<neighbor>& v, const long k, ForwardIterator query_point, const long first = -1, const long last = -1, const double epsilon = 0);	
		
		// search (and count) number of points within distance 'radius' from the query point
		// an unsorted vector v of neigbors is returned	
		template<class ForwardIterator>
		long search_range(vector<neighbor>& v, const double radius, ForwardIterator query_point, const long first = -1, const long last = -1);
						
		// count number of points within distance 'radius' from the query point ("correlation sum")
		template<class ForwardIterator>
		long count_range(const double radius, ForwardIterator query_point, const long first = -1, const long last = -1);

		template<class ForwardIterator>
		void count_range(long& count1, long& count2, const double radius1, const double radius2, ForwardIterator query_point, const long first, const long last);
};

// Advanced triangle inequaltity algorithm
template<class POINT_SET>
class ATRIA : public nearneigh_searcher<POINT_SET>
{
     using nearneigh_searcher<POINT_SET>::Nused;
     using nearneigh_searcher<POINT_SET>::table;
     using nearneigh_searcher<POINT_SET>::points;
     using nearneigh_searcher<POINT_SET>::err;

	protected :	
#ifdef MATLAB_MEX_FILE	
		long MINPOINTS;
#else
		const long MINPOINTS;
#endif	
		cluster root;

#ifdef MATLAB_MEX_FILE	
		neighbor* permutation_table;
#else
		neighbor* const permutation_table;
#endif			
	
		typedef typename POINT_SET::Metric METRIC;
		typedef searchitem SearchItem;
			
		priority_queue<SearchItem, vector<SearchItem>, searchitemCompare> search_queue;
		stack<SearchItem, vector<SearchItem> > SearchStack;		// used for range searches/counts 

		long total_clusters;
		long terminal_nodes;
		long total_points_in_terminal_node;

#ifdef PROFILE			
		unsigned long terminal_cluster_searched;
#endif 	
		void create_tree();
		void destroy_tree();
				
		pair<long, long> find_child_cluster_centers(const cluster* const c, neighbor* const Section, const long c_length);				
				
		long assign_points_to_centers(neighbor* const Section, const long c_length, pair<cluster*, cluster*> childs);	
					
		template<class ForwardIterator>		
		void search(ForwardIterator query_point, const long first, const long last, const double epsilon);	
		
	public:
		ATRIA(const POINT_SET& p, const long excl = 0, const long minpts = ATRIAMINPOINTS);
		~ATRIA();
			
		// search for k nearest neighbors of the point query_point, excluding all points
		// with indices between first and last from the search
		// returns number of nearest neighbor found and a sorted vector of neighbors (by reference)
		// an error in search_k_neighbors() will not result in an errorstate for the searcher ( see geterr() ) 
		template<class ForwardIterator>
		long search_k_neighbors(vector<neighbor>& v, const long k, ForwardIterator query_point, const long first = -1, const long last = -1, const double epsilon = 0);	
						
		// count number of points within distance 'radius' from the query point ("correlation sum")
		template<class ForwardIterator>
		long count_range(const double radius, ForwardIterator query_point, const long first = -1, const long last = -1);
			
		// search (and count) number of points within distance 'radius' from the query point
		// an unsorted vector v of neigbors is returned
		template<class ForwardIterator>
		long search_range(vector<neighbor>& v, const double radius, ForwardIterator query_point, const long first = -1, const long last = -1);
				
		inline double data_set_radius() const { return root.Rmax; };
		inline long total_tree_nodes() const { return total_clusters; };
		
#ifdef MATLAB_MEX_FILE
		// In case we have a matlab mex-file, offer a pair of functions to
		// store/retrieve an object of type ATRIA 
		ATRIA(const POINT_SET& p, const mxArray* inStructArr);	// create an ATRIA object from data stored in struct array inStructArr
		mxArray* store();		// store an ATRIA object into an MATLAB struct array
#endif		
};


template<class POINT_SET>
nearneigh_searcher<POINT_SET>::nearneigh_searcher(const POINT_SET& p, const long excl)
: err(0), points(p), Nused(points.size() - excl)
#ifdef PROFILE
	, points_searched(0), number_of_queries(0)
#endif
{	
	if ((Nused < 1) || (excl < 0))
	{
		cerr << "Wrong parameters for nearest neighbour search" << endl;
		cerr << "Nused : "  << Nused << endl;
		err = 1;
		return;
	}
}

template<class POINT_SET>
nearneigh_searcher<POINT_SET>::~nearneigh_searcher()
{
#ifdef PROFILE
	if (number_of_queries == 0)
		cout << "No queries were done" << endl;
	else
		cout << "Average percentage of points searched " 
		     << (100.0*(double)points_searched)/((double)Nused*number_of_queries)
			 << "% (" << ((long)ceil(((double)points_searched)/number_of_queries)) << ")" << endl;
#endif
}


template<class POINT_SET>
Brute<POINT_SET>::Brute(const POINT_SET& p, const long excl) 
:  nearneigh_searcher<POINT_SET>(p, excl)
{            
     if (this->err) {
	  cerr << "Brute : Error initializing parent object" << endl;
	  return;
     }
}

template<class POINT_SET> template<class ForwardIterator>
void Brute<POINT_SET>::search(ForwardIterator query_point, const long first, const long last, const double epsilon) 
{

// 	for (long j=0; j < Nused; j++) {
// 		if ((j < first) || (j > last)) test(j, query_point, table.highdist());
// 		
// 	}	
	long j;
	
	for (j=0; j <= first-1; j++) {
		test(j, query_point, table.highdist());
	}
	for (j=std::max(last+1, first); j < Nused; j++) {
		test(j, query_point, table.highdist());
	}		
}

template<class POINT_SET> template<class ForwardIterator>
long Brute<POINT_SET>::search_k_neighbors(vector<neighbor>& v, const long k, ForwardIterator query_point, const long first, const long last, const double epsilon) 
{
	// FIXME : make shure table is empty 12.Okt.1998 cmerk
	
#ifdef PROFILE
	number_of_queries++;
#endif

	table.init_search(k);
	search(query_point, first, last, epsilon);
	
	return table.finish_search(v);	// append table items to v, v should be empty, afterwards table is empty
}


template<class POINT_SET> template<class ForwardIterator>
long Brute<POINT_SET>::search_range(vector<neighbor>& v, const double radius, ForwardIterator query_point, const long first, const long last) 
{
	long count = 0;

#ifdef PROFILE
	number_of_queries++;
#endif

	for (long j=0; j < Nused; j++)
		if ((j < first) || (j > last)) {
			const double d = points.distance(j, query_point);
				
			if (d<=radius) { 
				v.push_back(neighbor(j,d));
				count++;	
			}	
#ifdef PROFILE
			points_searched++;
#endif			
		}
	
	return count;
}

template<class POINT_SET> template<class ForwardIterator>
long Brute<POINT_SET>::count_range(const double radius,ForwardIterator query_point, const long first, const long last) 
{
	long count = 0;

#ifdef PROFILE
	number_of_queries++;
#endif

	for (register long j=0; j < Nused; j++)
		if ((j < first) || (j > last)) {	
			if (points.distance(j, query_point) <= radius) count++;	
#ifdef PROFILE
			points_searched++;
#endif			
		}
	
	return count;
}

template<class POINT_SET> template<class ForwardIterator>
void Brute<POINT_SET>::count_range(long& count1, long& count2, const double radius1, const double radius2, ForwardIterator query_point, const long first, const long last) 
{
	count1 = 0;
	count2 = 0;
	
#ifdef PROFILE
	number_of_queries++;
#endif

	for (register long j=0; j < Nused; j++)
		if ((j < first) || (j > last)) {	
			const double d = points.distance(j, query_point);
			if (d <= radius1) count1++;
			if (d <= radius2) count2++;	
#ifdef PROFILE
			points_searched++;
#endif			
		}
}

template<class POINT_SET>
ATRIA<POINT_SET>::ATRIA(const POINT_SET& p, const long excl, const long minpts) 
	:	nearneigh_searcher<POINT_SET>(p, excl), root(1,Nused-1), MINPOINTS(minpts), 
		permutation_table(new neighbor[Nused]), total_clusters(1),
		total_points_in_terminal_node(0), terminal_nodes(0)
{
#ifdef VERBOSE
	cout << "ATRIA Constructor" << endl;
	cout << "Size of point set : " << p.size() << "  points of dimension " << p.dimension() << endl;
	cout << "Number of points used : " << Number_of_Points() << endl;
	cout << "MINPOINTS : " << MINPOINTS << endl;
#endif
#ifdef PROFILE			
	terminal_cluster_searched = 0;
#endif 	    
	        
	if (err) {
		cerr << "Error initializing parent object" << endl;
		return;
	}
	
	if (permutation_table == 0) { 
		cerr << "Out of memory" << endl; 
		err = 1; 
		return;
	}
				
	create_tree();
		
#ifdef VERBOSE
	cout << "Created tree structure for ATRIA searcher" << endl;
#endif
}

template<class POINT_SET>
ATRIA<POINT_SET>::~ATRIA()
{
#ifdef VERBOSE
	cout << "ATRIA Destructor" << endl;
#endif

#ifdef PROFILE
	cout << "Total_clusters : " << total_clusters << endl;
	cout << "Total number of points in terminal nodes : " << total_points_in_terminal_node << endl;
	cout << "Average number of points in a terminal node : " << ((double)total_points_in_terminal_node)/terminal_nodes << endl;
	cout << "Average number of terminal nodes visited : " << ((double)terminal_cluster_searched)/number_of_queries << endl;
#endif	

	destroy_tree();

	delete[] permutation_table;
}

template<class POINT_SET>
pair<long, long> ATRIA<POINT_SET>::find_child_cluster_centers(const cluster* const c, neighbor* const Section, const long length)
{
	pair<long, long> centers(-1, -1);
	
	if (c->Rmax == 0) { 		// if all data points seem to be identical
#ifdef VERBOSE		
		cout << "ATRIA : Data seem to be singular, search may be very inefficient" << endl;
#endif
#ifdef PROFILE
		terminal_nodes++;
		total_points_in_terminal_node += length;
#endif   

		return centers; // indicate that there's no need to further divide this data set
	}

	long index = 0;
	long center_right = Section[index].index(); 	
	double dist = Section[index].dist(); // points.distance(c->center, Section[index].index()); 
			
	// Compute right center, the point that is farthest away from the c->center
	
	for (long i=1; i < length; i++) {
		const double d = Section[i].dist();
 		if (d>dist) {
 			dist = d;
 			center_right = Section[i].index();
			index = i;
		}
 	}
	
	centers.second = center_right;
	
	// move this center the the last (rightmost) element of this Section
	My_Utilities::swap(Section, index, length-1);
	
	// compute left center, the point that is farthest away from the center_right 

	index = 0;
	long center_left = Section[index].index();	
	dist = points.distance(center_right, Section[index].index());
	Section[index].dist() = dist;
	
 	for (long i=1; i < length-1; i++) {			
 		const double d = points.distance(center_right, Section[i].index());
		Section[i].dist() = d;
 		if (d > dist) {
			dist = d;
 			center_left = Section[i].index();
			index = i;
 		}
 	}
	
	// move this center the the first (leftmost) element of this Section
	My_Utilities::swap(Section, index, 0);
	
	centers.first = center_left;
	
	return centers;
}

// assign each point to the nearest center, using a kind of quicksort like sorting procedure
template<class POINT_SET>
long ATRIA<POINT_SET>::assign_points_to_centers(neighbor* const Section, const long c_length, pair<cluster*, cluster*> childs)	
{
	const long center_left  = childs.first->center;

	register long i = 0; 
	register long j = c_length-1;
	
	// maximal distance fron one cluster's center to points belonging to this cluster
	double Rmax_left = 0;
	double Rmax_right = 0;
	
	// minimal gap in distances to both cluster centers
	double g_min_left = INFINITY;
	double g_min_right = INFINITY;
	
	while(1) {
		short i_belongs_to_left = 1;
		short j_belongs_to_right = 1;

		while(i+1 < j) {
			i++;
			const double dl = points.distance(center_left, Section[i].index());
			// ruse information instead of calculating dr = points.distance(center_right, Section[i].index());			
			const double dr = Section[i].dist(); 
			
			if (dl > dr) {	
				// point belongs to the right corner
				const double diff = dl - dr;
				
				Section[i].dist() = dr;	
				i_belongs_to_left = 0;
				
				g_min_right = std::min(g_min_right, diff);
				Rmax_right = std::max(Rmax_right, dr);
				break;
			} 
			
			// point belongs to the left corner
			const double diff = dr - dl;
			
			Section[i].dist() = dl;	
			
			g_min_left = std::min(g_min_left, diff);
			Rmax_left = std::max(Rmax_left, dl);
		}

		while(j-1 > i) { // either we reached the start of the array or this element was already checked by the i loop
			--j;

			const double dr = Section[j].dist(); // points.distance(center_right, Section[j].index());
			const double dl = points.distance(center_left, Section[j].index());

			if (dr >= dl) {
				// point belongs to the left corner
				const double diff = dr - dl;
				
				Section[j].dist() = dl;	
				j_belongs_to_right = 0;
				
				g_min_left = std::min(g_min_left, diff);
				Rmax_left = std::max(Rmax_left, dl);
				break;
			}
			
			// point belongs to the right corner
			const double diff = dl - dr;
			
			Section[j].dist() = dr;	
			
			g_min_right = std::min(g_min_right, diff);
			Rmax_right = std::max(Rmax_right, dr);
		}			

		if (i == j-1) {
			if ((!i_belongs_to_left) && (!j_belongs_to_right)) {
			        My_Utilities::swap(Section, i,j);		
			} else if (!i_belongs_to_left) {
				i--; j--;
			} else if (!j_belongs_to_right) {	
				i++; j++;
			}
			break;		// finished walking through the array
		} else {
		       My_Utilities::swap(Section, i,j);
		}
	}
	
	childs.first->g_min = g_min_left;
	childs.first->Rmax = Rmax_left;
		
	childs.second->g_min = g_min_right;
	childs.second->Rmax = Rmax_right;
	
	return j;
}				

template<class POINT_SET>
void ATRIA<POINT_SET>::create_tree()
{
	register long k;
	
	if (err) return;
	
	stack<cluster_pointer, cluster_pointer_vector> Stack;	// used for tree construction

	// select random center for root cluster, move this to first position of the indices array
	
	root.center = randindex(Nused);
	permutation_table[0] = neighbor(root.center, 0);

#ifdef VERBOSE
	cout << "Root center : " << root.center << endl;
	cout << "Root starting index  : " << root.start << endl;
	cout << "Root length : " << root.length << endl;
#endif   

	root.Rmax = 0;
	for (k=0; k < root.center; k++) {
		const double d = points.distance(k, root.center);
		permutation_table[k+1] = neighbor(k, d);
		if (d > root.Rmax)
			root.Rmax = d;
	}
	for (k=root.center+1; k < Nused; k++) {
		const double d = points.distance(k, root.center);
		permutation_table[k] = neighbor(k, d);
		if (d > root.Rmax)
			root.Rmax = d;		
	}

	// now create the tree
	
	Stack.push(&root);	// push root cluster on Stack
	
	while(!Stack.empty())	
	{
		cluster* const c = Stack.top(); Stack.pop();
		
		const long c_start = c->start;
		const long c_length = c->length;
		
		neighbor* const Section = permutation_table + c_start;
		
		if (c->length > MINPOINTS) {		// Further divide this cluster ?
			pair<long, long> new_child_centers = find_child_cluster_centers((const cluster*) c, Section, c_length);
			
			if ((new_child_centers.first == -1) || (new_child_centers.second == -1))
				continue;	// cluster could not be divided further and has already been marked as terminal cluster/node
			
			c->left = new cluster(new_child_centers.first);
			c->right = new cluster(new_child_centers.second);
					
			// create two subclusters and set properties
			const long j = assign_points_to_centers(Section, c_length, pair<cluster*, cluster*>(c->left,c->right));			
				
			c->left->start = c_start+1; 	// leave centers out
			c->left->length = j-1;
			
			c->right->start = c_start + j;
			c->right->length = c_length - j - 1;
			
			// process new subclusters (use stacks to avoid recursive call of this function)

 			Stack.push(c->right);
			Stack.push(c->left);

			total_clusters += 2;			
#ifdef VERBOSE		
 			cout << c_start << " " << c_length << " " << c->Rmax << " ";
 			cout << "L " << (c->left)->length << "  ";
 			cout << "R " << (c->right)->length << endl;
#endif				
		} 
		else {		// this is going to be a terminal node
			c->Rmax = - c->Rmax;	// a Rmax value <= 0 marks this cluster as a terminal node of the search tree

#ifdef PROFILE
			terminal_nodes++;
			total_points_in_terminal_node += c_length;
#endif 
#ifdef VERBOSE
 			cout << "Terminal node " << c->length << "  " << c->center << "  " << c->Rmax << endl;
#endif			
		}
	}
}

template<class POINT_SET>
void ATRIA<POINT_SET>::destroy_tree()
{	
	stack<cluster_pointer, cluster_pointer_vector> Stack;	// used for tree construction
			
	if (!root.is_terminal()) {
		Stack.push(root.left);
		Stack.push(root.right);
	}
	
	while (!Stack.empty()) {
		cluster* c = Stack.top(); Stack.pop();
		
		if (c != 0) {
			if (!c->is_terminal()) {
				Stack.push(c->left);
				Stack.push(c->right);
			}
			delete c;
		}
	}
}


template<class POINT_SET> 
template<class ForwardIterator>
long ATRIA<POINT_SET>::search_k_neighbors(vector<neighbor>& v, const long k, ForwardIterator query_point, const long first, const long last, const double epsilon) 
{
#ifdef PROFILE
	number_of_queries++;
#endif

	table.init_search(k);

	search(query_point, first, last, epsilon);
	
	return table.finish_search(v);	// append table items to v, v should be empty, afterwards table is empty
}

template<class POINT_SET>
template<class ForwardIterator>
void ATRIA<POINT_SET>::search(ForwardIterator query_point, const long first, const long last, const double epsilon) 
{
#ifdef PROFILE
	points_searched++;
#endif	
	const double root_dist = points.distance(root.center, query_point);
	
	while(!search_queue.empty()) search_queue.pop();	// clear search queue
				
	// push root cluster as search item into the PR-QUEUE
	search_queue.push(SearchItem(&root, root_dist));
	
	while(!search_queue.empty()) 
	{	
		const SearchItem si = search_queue.top(); search_queue.pop();
		const cluster* const c = si.clusterp();

		if ((table.highdist() > si.dist()) && ((c->center < first) || (c->center > last)))
			table.insert(neighbor(c->center, si.dist()));

		if (table.highdist() >= si.d_min() * (1.0 + epsilon)) {	// approximative (epsilon > 0) queries are supported			
			if (c->is_terminal())  {	
				const neighbor* const Section = permutation_table + c->start;				
#ifdef PROFILE
				terminal_cluster_searched++;
#endif

				if (c->Rmax == 0.0) {
					for (long i=0; i < c->length; i++) {
						const long j = Section[i].index();		

						if (table.highdist() <= si.dist()) 
							break;

						if ((j < first) || (j > last)) 
							table.insert(neighbor(j,si.dist()));
					}	
				} else {
					for (long i=0; i < c->length; i++) { 
						const long j = Section[i].index();		

						if ((j < first) || (j > last)) {
							if (table.highdist() > fabs(si.dist() - Section[i].dist()))
								test(j, query_point, table.highdist());	
						}				
					}
				}
			}
			else {				// this is an internal node
				const double dl = points.distance(c->left->center, query_point);
				const double dr = points.distance(c->right->center, query_point);
#ifdef PROFILE
				points_searched += 2;
#endif		
				// create child cluster search items
				SearchItem si_left = SearchItem(c->left, dl, dr, si);	
				SearchItem si_right = SearchItem(c->right, dr, dl, si);  

				// priority based search				
				search_queue.push(si_right);		
				search_queue.push(si_left);
			}	
		}
	}	
}

template<class POINT_SET> template<class ForwardIterator>
long ATRIA<POINT_SET>::search_range(vector<neighbor>& v, const double radius, ForwardIterator query_point, const long first, const long last) 
{
	long count = 0;

#ifdef PROFILE
	number_of_queries++;
	points_searched++;
#endif

	while (!SearchStack.empty()) SearchStack.pop();	// make shure stack is empty

	SearchStack.push(SearchItem(&root, points.distance(root.center, query_point)));

	while (!SearchStack.empty()) 
	{
		const SearchItem si = SearchStack.top();			
		SearchStack.pop();
		
		if (radius >= si.d_min()) {		
			const cluster* const c = si.clusterp();
		
			if (((c->center < first) || (c->center > last)) && (si.dist() <= radius)) {
				v.push_back(neighbor(c->center, si.dist()));
				count++;
			}
			
			if (c->is_terminal())  {	// this is a terminal node
				const neighbor* const Section = permutation_table + c->start;
				
				if (c->Rmax == 0.0) {		// cluster has zero radius, so all points inside will have the same distance to q
					if (radius >= si.dist()) {					
						for (long i=0; i < c->length; i++) {
							const long j = Section[i].index();		

							if ((j < first) || (j > last)) {
								v.push_back(neighbor(j,si.dist()));
								count++;
							}		
						}
					}	
				} else {
					for (long i=0; i < c->length; i++) { 
						const long j = Section[i].index();		

						if ( ((j < first) || (j > last)) && (radius >= fabs(si.dist() - Section[i].dist())) ) {
#ifdef PARTIAL_SEARCH
							const double d = points.distance(j, query_point, radius);
#else
							const double d = points.distance(j, query_point);
#endif		
							if (d <= radius) {
								v.push_back(neighbor(j,d));
								count++;
							}
#ifdef PROFILE
							points_searched++;
#endif					
						}
					}
				}
#ifdef PROFILE
					terminal_cluster_searched++;
#endif				
			}
			else {				// this is an internal node
				const double dl = points.distance(c->left->center, query_point);
				const double dr = points.distance(c->right->center, query_point);			
#ifdef PROFILE
				points_searched += 2;
#endif
				const SearchItem x = SearchItem(c->left, dl, dr, si);
				const SearchItem y = SearchItem(c->right, dr, dl, si);			
			
				SearchStack.push(x);
				SearchStack.push(y);				
			}
		}
	}
	
	return count;
}

template<class POINT_SET> template<class ForwardIterator>
long ATRIA<POINT_SET>::count_range(const double radius, ForwardIterator query_point, const long first, const long last) 
{
	long count = 0;

#ifdef PROFILE
	number_of_queries++;
	points_searched++;
#endif

	while (!SearchStack.empty()) SearchStack.pop();	// make shure stack is empty

	SearchStack.push(SearchItem(&root, points.distance(root.center, query_point)));

	while (!SearchStack.empty()) 
	{
		const SearchItem si = SearchStack.top();
			
		SearchStack.pop();
		
		if (radius >= si.d_min()) {		
			const cluster* const c = si.clusterp();

			if (((c->center < first) || (c->center > last)) && (si.dist() <= radius)) {
				count++;
			}
				
			if (c->is_terminal())  {	// this is a terminal terminal node
				const neighbor* const Section = permutation_table + c->start;

				if (c->Rmax == 0.0) {		// cluster has zero radius, so all points inside will have the same distance to q
					if (radius >= si.dist()) {					
						for (long i=0; i < c->length; i++) {
							const long j = Section[i].index();		

							if ((j < first) || (j > last)) {
								count++;
							}		
						}
					}	
				} else {
					for (long i=0; i < c->length; i++) {
						const long j = Section[i].index();		/* index of Vergleichspunkt */

						if ( ((j < first) || (j > last)) && (radius >= fabs(si.dist() - Section[i].dist())) ) {
#ifdef PARTIAL_SEARCH
							if (points.distance(j, query_point, radius) <= radius) 
								count++;
#else
							if (points.distance(j, query_point) <= radius) 
								count++;
#endif		

#ifdef PROFILE
							points_searched++;
#endif					
						}
					}
				}
#ifdef PROFILE
					terminal_cluster_searched++;
#endif				
			}
			else {				// this is an internal node
				const double dl = points.distance(c->left->center, query_point);
				const double dr = points.distance(c->right->center, query_point);			
#ifdef PROFILE
				points_searched += 2;
#endif
				const SearchItem x = SearchItem(c->left, dl, dr, si);
				const SearchItem y = SearchItem(c->right, dr, dl, si);			
			
				SearchStack.push(x);
				SearchStack.push(y);				
			}
		}
	}
	
	return count;
}

#ifdef MATLAB_MEX_FILE
#include "nn2matlab.h"
#endif			// ifdef MATLAB_MEX_FILE

#endif			// ifdef NEARNEIGH_SEARCH_H
