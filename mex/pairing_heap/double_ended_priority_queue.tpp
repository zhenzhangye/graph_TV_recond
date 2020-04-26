/* The double ended priority queue definition file
 *
 * here we use leaf correspondence.
 *
 */

// constructor: create the min and max heaps.
template <typename T>
DoubleEndedPriorityQueue<T>::DoubleEndedPriorityQueue(){
	min_heap_ = new PairingHeap<T>();
	min_heap_ ->setTypeOfHeap(true);
	max_heap_ = new PairingHeap<T>();
	max_heap_ ->setTypeOfHeap(false);
	buffer_ = NULL;
}

// deconstructor: delete two heaps.
template <typename T>
DoubleEndedPriorityQueue<T>::~DoubleEndedPriorityQueue(){
	delete min_heap_;
	delete max_heap_;
}

template <typename T>
bool DoubleEndedPriorityQueue<T>::empty(){
	if(buffer_==NULL)
		return min_heap_->size()+max_heap_->size() == 0;
	else
		return false;
}

template <typename T>
size_t DoubleEndedPriorityQueue<T>::size(){
	if(buffer_ == NULL)
		return min_heap_->size() + max_heap_->size();
	else
		return min_heap_->size() + max_heap_->size() + 1;
}

// clear the max and min heap.
template <typename T>
void DoubleEndedPriorityQueue<T>::ClearHeap(){
	while(!empty()){
		RemoveMin();
	}
	if(buffer_ != NULL){
		delete buffer_;
		buffer_ = NULL;
	}
}

// clean the size. don't earse the memory.
template <typename T>
void DoubleEndedPriorityQueue<T>::ClearSize(){
	buffer_ = NULL;
	min_heap_->ClearSize();
	max_heap_->ClearSize();
}

template <typename T>
T DoubleEndedPriorityQueue<T>::getMin(){
	// if min heap is empty, return the node in buffer.
	if (min_heap_->empty()){
		return buffer_->getKey();
	}else{
		// if buffer is empty, return the root node in min heap
		if(buffer_ == NULL)
			return min_heap_->getRootKey();
		else
			// otherwise, we have to compare the root node and buffer.
			if (*(buffer_->getKey())<=*(min_heap_->getRootKey()))
				return buffer_->getKey();
			else
				return min_heap_->getRootKey();
	}
}

template <typename T>
T DoubleEndedPriorityQueue<T>::getMax(){
	// if max heap is empty, return the node in buffer.
	if (max_heap_->empty()){
		return buffer_->getKey();
	}else{
		// if buffer is empty, return the root node in max heap.
		if(buffer_ == NULL)
			return max_heap_->getRootKey();
		else
			// otherwise, we have to compare the root node and buffer.
			if (*(buffer_->getKey())>=*(max_heap_->getRootKey()))
				return buffer_->getKey();
			else
				return max_heap_->getRootKey();
	}
}

template <typename T>
void DoubleEndedPriorityQueue<T>::Insert(T new_element){
	DataNode<T>* new_node = new DataNode<T>(new_element);
	// if buffer is empty, insert it into buffer.
	if (buffer_ == NULL){
		buffer_ = new_node;
	}else{
		DataNode<T> *smaller,*larger;
		// otherwise, compare the node with buffer.
		if(*(new_node->getKey()) < *(buffer_->getKey())){
			smaller = new_node;
			larger = buffer_;
		}else{
			smaller = buffer_;
			larger = new_node;
		}
		// the smaller node is inserted to min heap.
		min_heap_->Insert(smaller);
	
		// if the type of smaller node is a leaf node
		// insert the larger to max heap and create correspondency.
		if (smaller->type()){
			max_heap_->Insert(larger);
			smaller->setCorrespond(larger);
			larger->setCorrespond(smaller);
			buffer_ = NULL;
		}else{
			// otherwise, insert it into buffer.
			buffer_ = larger;
		}
	}
}

// merge two queues.
template <typename T>
void DoubleEndedPriorityQueue<T>::Merge(DoubleEndedPriorityQueue<T>* rhs){
	// min heap merges with min one.
	min_heap_->Merge(rhs->getMinHeap());
	rhs->getMinHeap()->setMerged(true);

	// max heap merges with max one.
	max_heap_->Merge(rhs->getMaxHeap());
	rhs->getMaxHeap()->setMerged(true);

	// special operation for buffer.
	if(buffer_==NULL && rhs->getBuffer() != NULL){
		buffer_ = rhs->getBuffer();
	}else if(buffer_!=NULL && rhs->getBuffer() != NULL){
		BufferInsert(rhs->getBuffer());
	}
}

// Insert a value into the buffer. Treat it as inserting a new node.
template <typename T>
void DoubleEndedPriorityQueue<T>::BufferInsert(DataNode<T>* new_node){
	DataNode<T> *smaller, *larger;
	if(*(new_node->getKey()) < *(buffer_->getKey())){
		smaller = new_node;
		larger = buffer_;
	}else{
		smaller = buffer_;
		larger = new_node;
	}
	min_heap_->Insert(smaller);
	if (smaller->type()){
		max_heap_->Insert(larger);
		smaller->setCorrespond(larger);
		larger->setCorrespond(smaller);
		buffer_ = NULL;
	}else{
		buffer_ = larger;
	}
}

template <typename T>
void DoubleEndedPriorityQueue<T>::RemoveMax(){
//case 1: buffer is empty, root corresponds to a non-leaf node OR corresponds to NULL.
//				#1 delete root directly and set correspond to NULL.
//case 2: buffer is empty, root corresponds to a leaf node.
//				2.1 the parent of correspond will become a leaf.
//						2.1.1 the parent has a correspond or is NULL. #2.1.1 move the correspond to buffer.
//						2.1.2 the parent doesn't have a correspond.   #2.1.2 set correspond and insert correspond to min heap.
//				2.2 the parent of correspond will not become a leaf.
//						#2.2 move the correspond to buffer.
//case 3: buffer is not empty.
//				3.1 compare buffer with root. If buffer is lager, delete buffer.
//				3.2 root corresponds to a non-leaf node OR corresponds to NULL.
//				#3.2 set correspond to NULL. delete root.
//				3.3 root corresponds to a leaf node.
//					3.3.1 if buffer is larger than root corresponce, insert buffer into max heap, set correspond. buffer = null.
//					3.3.2 if not, remove corresponce from min and insert it into max.
//						#3.3.2.1 if parent or the inserted becomes a leaf, and parent doesn't have corresponce, set correspond.
//						#3.3.2.2 if inserted becomes a leaf and parent has corresponce, insert buffer to min and set correspond.
//						#3.3.2.2 Otherwise, set inserted corresponce = null.

	if(buffer_ == NULL){
		DataNode<T>* y = max_heap_->getRoot();
		DataNode<T>* correspond = y->getCorrespond();
		max_heap_->DeleteRoot();

		if(correspond ==NULL){
			return;
		}

		if(!correspond->type()){
			correspond->setCorrespond(NULL);
		}else{
			DataNode<T>* p = correspond->getParent();
			y = min_heap_ -> Remove(correspond);
			y -> setCorrespond(NULL);

			if(p==NULL){
				buffer_ = y;
				y->setCorrespond(NULL);
				return;
			}


			if(p->getCorrespond()==NULL && p->type()){
				y -> setCorrespond(p);
				p -> setCorrespond(y);
				max_heap_->Insert(y);
				return;
			}else{
				buffer_ = y;
				y-> setCorrespond(NULL);
				return;
			}
		}
	}else{
		DataNode<T>* y = max_heap_->getRoot();

		if(y==NULL){
			delete buffer_;
			buffer_ = NULL;
			return;
		}

		if(*(y->getKey()) <= *(buffer_->getKey())){
			delete buffer_;
			buffer_ = NULL;
			return;
		}

		if (y->getCorrespond() == NULL){
			max_heap_->DeleteRoot();
			return;
		}

		DataNode<T>* correspond = y->getCorrespond();
		max_heap_->DeleteRoot();

		if (!correspond->type()){
			correspond->setCorrespond(NULL);
			return;
		}
		
		if (correspond->type()){
			if(*(buffer_->getKey())>=*(correspond->getKey())){
				max_heap_->Insert(buffer_);
				buffer_->setCorrespond(correspond);
				correspond->setCorrespond(buffer_);
				buffer_ = NULL;
				return;
			}else{
				DataNode<T>* p = correspond->getParent();
				if(p == NULL){
					min_heap_->Remove(correspond);
					min_heap_->Insert(buffer_);
					max_heap_->Insert(correspond);
					correspond->setCorrespond(buffer_);
					buffer_->setCorrespond(correspond);
					buffer_ = NULL;
					return;
				}

				min_heap_->Remove(correspond);
				max_heap_->Insert(correspond);
				if((p->type() || correspond->type()) && p->getCorrespond() == NULL){
					p->setCorrespond(correspond);
					correspond->setCorrespond(p);
					return;
				}else{
					if(correspond->type()){
						min_heap_->Insert(buffer_);
						buffer_->setCorrespond(correspond);
						correspond->setCorrespond(buffer_);
						buffer_ = NULL;
						return;
					}else{
						correspond->setCorrespond(NULL);
						return;
					}
				}
			}
		}
	}
}

// The same considerations as remove max. see above.
template <typename T>
void DoubleEndedPriorityQueue<T>::RemoveMin(){
	if(buffer_ == NULL){
		DataNode<T>* y = min_heap_->getRoot();
		DataNode<T>* correspond = y->getCorrespond();
		min_heap_->DeleteRoot();

		if(correspond ==NULL){
			return;
		}

		if(!correspond->type()){
			correspond->setCorrespond(NULL);
		}else{
			DataNode<T>* p = correspond->getParent();
			y = max_heap_ -> Remove(correspond);
			y->setCorrespond(NULL);

			if(p==NULL){
				buffer_ = y;
				y->setCorrespond(NULL);
				return;
			}

			if(p->getCorrespond()==NULL && p->type()){
				y -> setCorrespond(p);
				p -> setCorrespond(y);
				min_heap_->Insert(y);
				return;
			}else{
				buffer_ = y;
				y->setCorrespond(NULL);
				return;
			}
		}
	}else{

		DataNode<T>* y = min_heap_->getRoot();
		if(y==NULL){
			delete buffer_;
			buffer_ = NULL;
			return;
		}

		if(*(y->getKey()) >= *(buffer_->getKey())){
			delete buffer_;
			buffer_ = NULL;
			return;
		}

		if (y->getCorrespond() == NULL){
			min_heap_->DeleteRoot();
			return;
		}

		DataNode<T>* correspond = y->getCorrespond();
		min_heap_->DeleteRoot();

		if (!correspond->type()){
			correspond->setCorrespond(NULL);
			return;
		}
		
		if (correspond->type()){
			if(*(buffer_->getKey())<=*(correspond->getKey())){
				min_heap_->Insert(buffer_);
				buffer_->setCorrespond(correspond);
				correspond->setCorrespond(buffer_);
				buffer_ = NULL;
				return;
			}else{
				DataNode<T>* p = correspond->getParent();
				if(p == NULL){
					max_heap_->Remove(correspond);
					max_heap_->Insert(buffer_);
					min_heap_->Insert(correspond);
					correspond->setCorrespond(buffer_);
					buffer_->setCorrespond(correspond);
					buffer_ = NULL;
					return;
				}

				max_heap_->Remove(correspond);
				min_heap_->Insert(correspond);
				if((p->type() || correspond->type()) && p->getCorrespond() == NULL){
					p->setCorrespond(correspond);
					correspond->setCorrespond(p);
					return;
				}else{
					if(correspond->type()){
						max_heap_->Insert(buffer_);
						buffer_->setCorrespond(correspond);
						correspond->setCorrespond(buffer_);
						buffer_ = NULL;
						return;
					}else{
						correspond->setCorrespond(NULL);
						return;
					}
				}
			}
		}
	}
}

template <typename T>
PairingHeap<T>* DoubleEndedPriorityQueue<T>::getMinHeap(){
	return min_heap_;
}

template <typename T>
PairingHeap<T>* DoubleEndedPriorityQueue<T>::getMaxHeap(){
	return max_heap_;
}

template <typename T>
DataNode<T>* DoubleEndedPriorityQueue<T>::getBuffer(){
	return buffer_;
}

// test function: print the information of this queue.
template <typename T>
void DoubleEndedPriorityQueue<T>::PrintQueue(){
	std::cout<<"Buffer: ";
	if(buffer_!=NULL){
		std::cout<<buffer_->getKey()<<std::endl;
		std::cout<<"lambda, slope, shift: "<<buffer_->getKey()->lambda_<<" "<<buffer_->getKey()->slope_<<" "<<buffer_->getKey()->shift_<<std::endl;
	}else{
		std::cout<<"NULL"<<std::endl;
	}
	min_heap_->PrintHeap();
	max_heap_->PrintHeap();
}
