#ifndef TERNARY_SEARCH_TREE
#define TERNARY_SEARCH_TREE

#include <cmath>

#define TST_BUFFERS 128	// number of buffers
#define TST_BUFSIZE 512	// inital buffer size

template<class KEY>
struct Tnode {
    KEY splitkey;
	int level;				/* level of tree */
	unsigned long count;	/* count how often the key leading to this node exists */
    
	Tnode* lokid;
	Tnode* eqkid;
	Tnode* hikid;
};

// class ternary_search_tree stores multikey-data with fixed key length

template<class KEY, class Evaluater>
class ternary_search_tree
{
	protected:
		typedef Tnode<KEY>* Tptr;
		
		// these four variables are used to accelerate allocation of Tnode objects
		
		Tptr buf;			// pointer to current buffer
		long next_buf_size; // size of the buffer that will be allocated next
		long bufn;			// number of next buffer
		long freen;         // number of free nodes in current buffer
		
		Tptr freearr[TST_BUFFERS];	
		
		const long len; 	// key length 
		Tptr root;			// tree root node
		
	public:
		ternary_search_tree(const long keylength);
		~ternary_search_tree();
		int insert(const KEY* const key);		// insert key vector
		long total_nodes() const { return (TST_BUFSIZE * (pow(2, bufn) - 1) - freen); }
		
		void traverse(Evaluater& eval);
};

template<class KEY, class Evaluater>
ternary_search_tree<KEY, Evaluater>::ternary_search_tree(const long keylength) : bufn(0), freen(0), 
	root(0), next_buf_size(TST_BUFSIZE), len(keylength) {}


// return 0 on SUCCESS
template<class KEY, class Evaluater>
int ternary_search_tree<KEY, Evaluater>::insert(const KEY* const key)
{   
	KEY d;
	Tptr pp;
	
    Tptr* p = &root;
	long level = 0;			// level goes up to len-1 
	
    while(pp = *p) {		// as long as we encounter already exisiting nodes, we stay inside this while loop
		if ((d = key[level] - pp->splitkey) == 0) {		// go to next tree level
			pp->count++;
            p = &(pp->eqkid);
			if ((++level) == len) return 0;
        } else if (d < 0) {
            p = &(pp->lokid);	/* move left in the current level */
		}
        else {
            p = &(pp->hikid);	/* move right in the current level */
		}
    }
    for (;;) {	/* once we find a node that is not allocated (==0), we must create every next node */
		if (freen-- == 0) {
			if (bufn == TST_BUFFERS) {
				//mexErrMsgTxt("Ran out of available buffers for tree nodes");
				return -1;		// FAILURE
			}			
			buf = new Tnode<KEY>[next_buf_size]; 
			freearr[bufn++] = buf;
			freen = next_buf_size-1;
			next_buf_size *= 2; 	// double size of the next buffer (this keeps overall number of allocations small)
		}
		*p = buf++;
        pp = *p;
        pp->splitkey = key[level];
		pp->count = 1;					// this node is newly created, so count is set to one 
		pp->level = level;
        pp->lokid = pp->eqkid = pp->hikid = 0;
        if ((++level) == len) return 0;
        p = &(pp->eqkid);
    }
}

// traverse tree in an arbitrary order, execute the given function object on every node that is not empty
template<class KEY, class Evaluater>
void ternary_search_tree<KEY, Evaluater>::traverse(Evaluater& eval)
{   
	long buf_size = TST_BUFSIZE;		// the actual size of the buffer is doubling each iteration
	
	// traverse through all buffers that are completely filled 
	for (long i = 0; i < bufn; i++) {
		long number_nodes = buf_size;
		const Tptr b = (Tptr) freearr[i];
		
		if (i == bufn - 1)		// last buffer may not be completely filled
			number_nodes = buf_size-freen;
		
		buf_size *= 2;
		
		for (long j = 0; j < number_nodes; j++) {
			const Tptr p = b + j;	
			eval(p->count, p->level);
		}
	}
}

template<class KEY, class Evaluater>
ternary_search_tree<KEY, Evaluater>::~ternary_search_tree()
{   
    for (long i = 0; i < bufn; i++)
        delete[] freearr[i];
}

#endif
