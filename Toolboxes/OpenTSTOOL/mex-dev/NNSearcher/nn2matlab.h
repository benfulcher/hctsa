#ifndef NN2MATLAB_H
#define NN2MATLAB_H

// This header file is included from nearneigh_search.h
// It should not be included from any other file.
//
// This code enables to store/retrieve the preproccesing data (as object of type ATRIA)
// to a matlab (version 5.0 - 5.x) structure variable (to avoid multiple preproccesing).

class RetrieveItem	// this class is needed when retrieving the tree from an array-like structure
{
	public:
		long index_of_node_needs_creation; 		// (array)index of the cluster that has to be created
		cluster_pointer parent_needs_update;	// pointer to the parent node that needs an update
		int left_or_right; 						// a flag which indicates if this child	is a left or right node to the parent						
												// == 0 means : left child, == 1 means right child, == -1 means : there's no parent to update

};

typedef vector<RetrieveItem> RetrieveItem_vector; 

template<class POINT_SET>
ATRIA<POINT_SET>::ATRIA(const POINT_SET& p, const mxArray* inStructArr) // create an ATRIA object from data stored in struct array inStructArr
 : nearneigh_searcher<POINT_SET>(p, 0), root(0,Nused), MINPOINTS(0), permutation_table(0), total_clusters(1),
 total_points_in_terminal_node(0), terminal_nodes(0)
{	
	long i;
	 	
	err = 0;
		
	if (!mxIsStruct(inStructArr)) {
		err = 1;
		return;
	}
	
	// params
	
	const mxArray* tmparr = mxGetField(inStructArr, 0, "params");
	
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != 5)) {
		err = 1;
		return;
	}
	const double* tmpptr = mxGetPr(tmparr);
	
	Nused = (long) tmpptr[0]; 
	total_clusters = (long) tmpptr[1]; 
	terminal_nodes = (long) tmpptr[2]; 
	total_points_in_terminal_node = (long) tmpptr[3]; 
	MINPOINTS = (long) tmpptr[4];

	if ((Nused < 2) || ( Nused > p.size())) {
		err = 1;
		return;
	}
	
	permutation_table = new neighbor[Nused];
	if (permutation_table == 0) { 
		cerr << "Out of memory" << endl; 
		err = 1; 
		return;
	}
	// distances
	
	tmparr = mxGetField(inStructArr, 0, "distances");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != Nused)) {
		err = 1;
		return;
	}
	
	tmpptr = mxGetPr(tmparr);
	
	for (i=0; i < Nused; i++) permutation_table[i].dist() = tmpptr[i];	

	// indices
	
	tmparr = mxGetField(inStructArr, 0, "indices");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != Nused)) {
		err = 1;
		return;
	}
	tmpptr = mxGetPr(tmparr);
	
	for (i=0; i < Nused; i++) permutation_table[i].index() = (long) tmpptr[i];	// warning : double to long conversion
	
	// retrieve cluster tree
	
	const mxArray* tree = mxGetField(inStructArr, 0, "tree");
	
	if ((tree == NULL) || (!mxIsStruct(tree))) {
		err = 1;
		return;
	}

	tmparr = mxGetField(tree, 0, "center");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}
    double* center_arr = mxGetPr(tmparr);
	
	tmparr = mxGetField(tree, 0, "Rmax");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}	
    double* Rmax_arr = mxGetPr(tmparr);

	tmparr = mxGetField(tree, 0, "start");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}    
	double* start_arr = mxGetPr(tmparr);
	
	tmparr = mxGetField(tree, 0, "length");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}    
    double* length_arr = mxGetPr(tmparr);
	
	tmparr = mxGetField(tree, 0, "leftchild");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}    
    double* leftchild_arr = mxGetPr(tmparr);
	
	tmparr = mxGetField(tree, 0, "rightchild");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}    
    double* rightchild_arr = mxGetPr(tmparr);
	
	tmparr = mxGetField(tree, 0, "g");
	if ((tmparr == NULL) || (!mxIsDouble(tmparr)) || (mxGetM(tmparr) != total_clusters)) {
		err = 1;
		return;
	}    
    double* g_arr = mxGetPr(tmparr);	
	
	stack<RetrieveItem, RetrieveItem_vector> s;
	RetrieveItem y,z;
	
	root.center = (long) center_arr[0];
	root.Rmax = Rmax_arr[0];
	root.g_min = g_arr[0];
	
	if (root.is_terminal()) {
		root.start = (long) start_arr[0];	
		root.length = (long) length_arr[0];
	} else {
		y.index_of_node_needs_creation = (long) leftchild_arr[0];
		y.parent_needs_update = &root;
		y.left_or_right = 0;
		
		z.index_of_node_needs_creation = (long) rightchild_arr[0];
		z.parent_needs_update = &root;
		z.left_or_right = 1;
		
		s.push(y);
		s.push(z);
	}	
	
	while (!s.empty()) 
	{
		const RetrieveItem x = s.top(); s.pop();
		const cluster_pointer c = new cluster();
		
		if ((x.index_of_node_needs_creation < 0) || (x.index_of_node_needs_creation >= total_clusters)) {
			err = 1;
			return;
		}
		
		c->center = (long) center_arr[x.index_of_node_needs_creation];
		
		if ((c->center < 0) || (c->center >= Nused)) {
			err = 1;
			return;
		}
		
		c->Rmax = Rmax_arr[x.index_of_node_needs_creation];
		c->g_min = g_arr[x.index_of_node_needs_creation];
		
		if (c->is_terminal()) {
			c->start = (long) start_arr[x.index_of_node_needs_creation];
			if ((c->start < 0) || (c->start >= Nused)) {
				err = 1;
				return;
			}
			
			c->length = (long) length_arr[x.index_of_node_needs_creation];
			if ((c->length < 0) || (c->length >= Nused)) {
				err = 1;
				return;
			}
		} else {
			y.index_of_node_needs_creation = (long) leftchild_arr[x.index_of_node_needs_creation];
			y.parent_needs_update = c;
			y.left_or_right = 0;

			z.index_of_node_needs_creation = (long) rightchild_arr[x.index_of_node_needs_creation];
			z.parent_needs_update = c;
			z.left_or_right = 1;

			s.push(y);
			s.push(z);
		}	
		
		switch(x.left_or_right) {
			case 0 :
				x.parent_needs_update->left = c;
				break;
			case 1 :
				x.parent_needs_update->right = c;
				break;
			default :
				; 
				break;
		}
	}
}

class StoreItem	// this class is needed when storing tree into an array-like structure
{
	public:
		cluster_pointer child_to_store;		// a pointer to a cluster that still needs to be stored in the array
		
		unsigned long parent_index; 		// index (0,1,2...) of the parent cluster that is already stored, but needs
											// an update because the array-index at which this child is stored
											// in the array was not known at the time the parent was stored
		
		int left_or_right; 					// a flag which indicates if this child	is a left or right node to the parent						
											// == 0 means : left child, == 1 means right child, == -1 means : there's no parent to update
};
		
typedef vector<StoreItem> StoreItem_vector; 	// needed for SUN CC 4.1

template<class POINT_SET>
mxArray* ATRIA<POINT_SET>::store()	// store an ATRIA object into an MATLAB struct array
{
	long i;
	
	mxArray* output;
	mxArray* tree;
	
	const char *fieldnames[] = {"params",  // (long) vector containing : Nused, total_clusters, terminal_nodes, average_points_in_terminal_node, MINPOINTS
								"distances", //	(double) vector containing distances
								"indices",	// (double) vector containing indices
								"optional", // cell array for optional parameters
								"tree"};	// structure array containing tree of clusters
	
	const char *tree_fieldnames[] = {"center", // (long) vector of center indices
									"Rmax",  // (double) vector of Rmax values
									"g",		// (double) vector
									"start",	// (long) vector of start indices, only valid for terminal nodes	
									"length",		// (long) vector of length values, only valid for terminal nodes	
									"leftchild",	// (long) vector of left child nodes, only valid for non-terminal nodes
									"rightchild"}; // (long) vector of right child nodes, only valid for non-terminal nodes	
	
	double* tmpptr;
	mxArray* tmparr;
			
	output = mxCreateStructMatrix(1, 1, 5, fieldnames);

	// first store scalar paramters of the ATRIA object;
	
	tmparr = mxCreateDoubleMatrix(5, 1, mxREAL);
	tmpptr = (double *) mxGetPr(tmparr);
	
	tmpptr[0] = Nused; tmpptr[1] = total_clusters; tmpptr[2] = terminal_nodes; tmpptr[3] = total_points_in_terminal_node; 
	tmpptr[4] = MINPOINTS;
	
	mxSetField(output, 0, "params", tmparr);

	// store distances
	
	tmparr = mxCreateDoubleMatrix(Nused, 1, mxREAL);
	tmpptr = (double *) mxGetPr(tmparr);
	
	for (i=0; i < Nused; i++) tmpptr[i] = permutation_table[i].dist();	
	
	mxSetField(output, 0, "distances", tmparr);
	
	// store indices
	
	tmparr = mxCreateDoubleMatrix(Nused, 1, mxREAL);
	tmpptr = (double *) mxGetPr(tmparr);
	
	for (i=0; i < Nused; i++) tmpptr[i] = permutation_table[i].index();	// warning : long to double conversion
	
	mxSetField(output, 0, "indices", tmparr);
	
	// store cluster tree
	
	tree = mxCreateStructMatrix(1, 1, 7, tree_fieldnames);
	
	mxSetField(output, 0, "tree", tree);
	
	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* center_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "center", tmparr);
	
	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* Rmax_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "Rmax", tmparr);

	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* start_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "start", tmparr);
	
	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* length_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "length", tmparr);
	
	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* leftchild_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "leftchild", tmparr);
	
	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* rightchild_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "rightchild", tmparr);

	tmparr = mxCreateDoubleMatrix(total_clusters, 1, mxREAL);
	double* g_arr = mxGetPr(tmparr);
	mxSetField(tree, 0, "g", tmparr);

	stack<StoreItem, StoreItem_vector> s;
	StoreItem x;
	long next_free = 0;
	
	x.child_to_store = &root;
	x.parent_index = 0;
	x.left_or_right = -1;
	
	s.push(x);
	
	while (!s.empty()) 
	{
		x = s.top(); s.pop();
		
		center_arr[next_free] = x.child_to_store->center;
		Rmax_arr[next_free] = x.child_to_store->Rmax;
		g_arr[next_free] = x.child_to_store->g_min;
		
		leftchild_arr[next_free] = -1;
		rightchild_arr[next_free] = -1;
		start_arr[next_free] = -1;
		length_arr[next_free] = -1;
		
		if (x.child_to_store->is_terminal()) {
			start_arr[next_free] = x.child_to_store->start;
			length_arr[next_free] = x.child_to_store->length;
		} else {
			StoreItem y,z;
			
			y.child_to_store = x.child_to_store->left;
			y.parent_index = next_free;
			y.left_or_right = 0;
			
			z.child_to_store = x.child_to_store->right;
			z.parent_index = next_free;
			z.left_or_right = 1;
			
			s.push(y);
			s.push(z);
		}
		
		switch(x.left_or_right) {
			case 0 :
				leftchild_arr[x.parent_index] = next_free;
				break;
			case 1 :
				rightchild_arr[x.parent_index] = next_free;
				break;
			default :
				; 
				break;
		}
	
		next_free++;
	}
	
	if (next_free != total_clusters) {
		mexErrMsgTxt("Error storing ATRIA object, number of tree nodes does not match expected value");
		return output;
	}
	
	return output;
}

#endif
