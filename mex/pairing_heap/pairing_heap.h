/* the pairing heap file.
 *
 * used for double ended priority queue.
 *
 */
#ifndef __PAIRINGHEAP_H__ 
#define __PAIRINGHEAP_H__
#include <vector>
#include <queue>
#include <iostream>

// This node stores the relationship in heap and data information.
template<typename T> class DataNode;
template<typename T>
class DataNode{
public:

	// member functions, constructor and deconstructor.
	DataNode();
	DataNode(T key);
	~DataNode();
	
	bool						type(); //true if it is a leaf node.

	DataNode<T>*		getPrevious();
	DataNode<T>*		getSibling();
	DataNode<T>*  	getChild();
	DataNode<T>*  	getParent();
	DataNode<T>*		getCorrespond(); // get the corresponding node in another heap.
	T								getKey();

	void 						setCorrespond(DataNode<T>* node);

	DataNode<T>*		correspond;
	DataNode<T>*		previous;
	DataNode<T>*		parent;
	DataNode<T>*		sibling;
	DataNode<T>*		child;
	T								key;
	
	// define new operators
	bool						operator= ( DataNode<T>& rhs);
	bool						operator> ( DataNode<T>& rhs);
	bool						operator< ( DataNode<T>& rhs);
	bool						operator>= ( DataNode<T>& rhs);
	bool						operator<= ( DataNode<T>& rhs);
	bool						operator== ( DataNode<T>& rhs);
};

template<typename T>
class PairingHeap{

public:
	//constructor and deconstructor
	PairingHeap();
	~PairingHeap();

	//public member functions
	bool						empty();
	size_t 					size();
	bool    				type(); 				// type of current heap. true: min heap. false: max heap.
	DataNode<T>* 		Insert(DataNode<T>* new_node);
	void          	Merge(PairingHeap<T>* rhs);
	void 						DeleteRoot();
	DataNode<T>*		Remove(DataNode<T>* node);
	bool						getMerged(); // true if thie heap is merged into another one. Used while deleting
	void						ClearSize();

	T								getRootKey();

	PairingHeap<T>* getCorrespondHeap();
	void 						setCorrespondHeap(PairingHeap<T>* correspond_heap);

	void						setTypeOfHeap(bool type);
	DataNode<T>*		getRoot();
	void						setMerged(bool merge_value = true);

	
	//testing functions
	void 					PrintHeap();
	
private:
	// free memory.
	void  					ReclaimMemory(DataNode<T>* current);
	// merge related functions.
	DataNode<T>*		TwoPassMerge(DataNode<T>* current);
	void  					CompareAndLink(DataNode<T>* &first, DataNode<T>* second);

	// decrease/increase number of nodes in this heap.
	void						DecreaseNodeNumber(size_t number = 1);
	void						IncreaseNodeNumber(size_t number = 1);

	size_t 					number_of_nodes_;
	bool 						type_of_heap_;
	bool						merged_;
	DataNode<T>* 		root_;
	PairingHeap<T>* correspond_heap_;
};

#include "data_node.tpp"
#include "pairing_heap.tpp"
#endif
