/* the double ended priority queue head file
 *
 * define all the methods and variables for a double ended priority queue.
 * here, we use leaf correspondence.
 *
 */

#ifndef __DOUBLEENDEDPRIORITYQUEUE_H__
#define __DOUBLEENDEDPRIORITYQUEUE_H__
#include "pairing_heap.h"

template<typename T>
class DoubleEndedPriorityQueue{

public:
	//constructor and deconstructor
	DoubleEndedPriorityQueue();
	~DoubleEndedPriorityQueue();

	//public member functions.
	bool						empty();
	size_t					size();

	// clean the heap. free memories used.
	void						ClearHeap();
	void						ClearSize();
	T								getMin();
	T								getMax();
	void						Insert(T new_element);

	// merge two heaps
	void						Merge(DoubleEndedPriorityQueue<T>* depq);
	void						RemoveMin();
	void						RemoveMax();
	PairingHeap<T>* getMinHeap();
	PairingHeap<T>* getMaxHeap();

	// get the node in the buffer.
	DataNode<T>*		getBuffer();

	void						PrintQueue();

private:
	// insert node in the buffer.
	void						BufferInsert(DataNode<T>* new_node);
	bool						heap_created_;
	PairingHeap<T>* min_heap_;
	PairingHeap<T>* max_heap_;
	DataNode<T>*		buffer_;
};

#include "double_ended_priority_queue.tpp"
#endif
