#include "nn_aux.h"

long SortedNeighborTable::finish_search(vector<neighbor>& v)
{
	v.reserve(pq.size());

	while (!pq.empty()) { 
		v.push_back(pq.top()); 
		pq.pop(); 
	}

	reverse(v.begin(),v.end());
	return v.size();
}

void SortedNeighborTable::insert(const neighbor& x)
{
	pq.push(x);
	 
	while(pq.size() > NNR) {
		pq.pop();
	} 
	
	if (pq.size() < NNR) {
		hd = INFINITY;
	}	
	else {
		hd = pq.top().dist();
	}
}

#ifdef USE_OWN_CLUSTER_MEMORY_HANDLER
void* cluster::operator new(size_t size)
{
	cluster* p = headOfFreeList;    // p is now a pointer to the head of the free list

	// if p is valid, just move the list head to the
	// next element in the free list
	if (p) 
		headOfFreeList = p->right;
	else {
		// The free list is empty. Allocate a block of memory
		// big enough hold BLOCK_SIZE cluster objects
		
		cluster *newBlock = (cluster*) ::operator new(BLOCK_SIZE * sizeof(cluster));
		
		for (long i = 1; i < BLOCK_SIZE-1; ++i)
			newBlock[i].right = &newBlock[i+1];

		// terminate the linked list with a null pointer
		newBlock[BLOCK_SIZE-1].right = 0;

		// set p to front of list, headOfFreeList to
		// chunk immediately following
		p = newBlock;
		headOfFreeList = &(newBlock[1]);

		// cmerk April 99 : Fibonacci like scheme to increase BLOCK_SIZE each time a new block is demanded
		const long NEW_BLOCK_SIZE = BLOCK_SIZE + OLD_BLOCK_SIZE;
		OLD_BLOCK_SIZE = BLOCK_SIZE;
		BLOCK_SIZE = NEW_BLOCK_SIZE; 
	}

	return p;
}

// operator delete is passed a memory chunk, which,
// if it's the right size, is just added to the
// front of the list of free chunks
void cluster::operator delete(void *deadObject, size_t size)
{
	if (deadObject == 0) return;           

	cluster *carcass = (cluster*) deadObject;

	carcass->right = headOfFreeList;
	headOfFreeList = carcass;
}

cluster* cluster::headOfFreeList = 0;
    
long cluster::OLD_BLOCK_SIZE = 0;
long cluster::BLOCK_SIZE = 512;  
#endif
