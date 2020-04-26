//default constructor
template <typename T>
PairingHeap<T>::PairingHeap(){
	root_ = NULL;
	type_of_heap_ = true;
	number_of_nodes_ = 0;
	correspond_heap_ = NULL;
	merged_ = false;
}

//deconstructor, release all memory space
template <typename T>
PairingHeap<T>::~PairingHeap(){
	if(!merged_ || root_!=NULL){
		ReclaimMemory(root_);
		root_ = NULL;
		type_of_heap_ = true;
		number_of_nodes_ = 0;
		correspond_heap_ = NULL;
	}
}

//return if the heap is empty
template <typename T>
bool PairingHeap<T>::empty(){
	return number_of_nodes_ == 0;
}

//return the number of nodes in the heap
template <typename T>
size_t PairingHeap<T>::size(){
	return number_of_nodes_;
}

//return the type of this heap. true: min heap, false: max heap;
template <typename T>
bool PairingHeap<T>::type(){
	return type_of_heap_;
}

//decrease number of nodes by one.
template <typename T>
void PairingHeap<T>::DecreaseNodeNumber(size_t number){
	if (number_of_nodes_ < number){
		std::cout<<"DecreaseNodeNumber: nodes in heap is less than decreased number"<<std::endl;
	}else{
		number_of_nodes_ -= number;
	}
}

//Increase number of nodes by one.
template <typename T>
void PairingHeap<T>::IncreaseNodeNumber(size_t number){
	number_of_nodes_ += number;
}

//set if this heap is merged into another one
template <typename T>
void PairingHeap<T>::setMerged(bool merge_value){
	merged_ = merge_value;
}

//get the root node of this heap.
template <typename T>
DataNode<T>* PairingHeap<T>::getRoot(){
	if (empty()){
		return NULL;
	}
	return root_;
}

//get the root key of this heap.
template <typename T>
T PairingHeap<T>::getRootKey(){
	if (empty()){
		std::cout<<"getRootKey: This heap is empty"<<std::endl;
		T t;
		return t;
	}
	return root_ -> getKey();
}

//insert new element into the heap
template <typename T>
DataNode<T>* PairingHeap<T>::Insert(DataNode<T>* new_node){
	// update the size of the heap
	IncreaseNodeNumber();

	// if current heap is empty, new element is root.
	if (root_ == NULL){
		root_ = new_node;
	}else{
		CompareAndLink(root_, new_node);
	}

	return root_;
}

//Merge a new heap with current heap
template <typename T>
void PairingHeap<T>::Merge(PairingHeap<T>* rhs){
	if (rhs->type() != type_of_heap_){
		std::cout<<"Merge: type does not match"<<std::endl;
		return;
	}
	if (root_ == NULL && rhs -> getRoot() == NULL)
		return;
	if (root_ == NULL){
		root_ = rhs -> getRoot();
		number_of_nodes_ = rhs -> size();
		return;
	}
	if (rhs -> getRoot() == NULL){
		return;
	}

	rhs->setMerged();
	bool current_case = false;
	// two heaps are min ones and this heap has smaller key  OR
	// two heaps are max ones and this heap has larger key
	if ((type_of_heap_ && (*root_->getKey()) <= (*rhs->getRoot()->getKey())) ||
		 (!type_of_heap_ && (*root_->getKey()) >= (*rhs->getRoot()->getKey())))
		current_case = true;

	// this heap's root becomes new root
	if (current_case){
		rhs-> getRoot() -> sibling = root_ -> child;
		rhs-> getRoot() -> previous = root_;

		//update the child node with smaller key.
		if (root_ -> child != NULL){
			root_ -> child -> previous = rhs-> getRoot();
			root_ -> child -> parent = NULL;
		}
		
		//set new child node.
		root_ -> child = rhs-> getRoot();
		rhs-> getRoot() -> parent = root_;
	}else{
		root_ -> sibling = rhs-> getRoot() -> child;
		root_ -> previous = rhs-> getRoot();
	
		//update the child node with smaller key.
		if(rhs -> getRoot() -> child != NULL){
			rhs -> getRoot() -> child -> previous = root_;
			rhs -> getRoot() -> child -> parent = NULL;
		}

		//set new child node.
		rhs -> getRoot() -> child = root_;
		root_ -> parent = rhs ->getRoot();

		//defaultly, merge rhs to current heap.
		root_ = rhs -> getRoot();
	}

	IncreaseNodeNumber(rhs->size());
}

//delete root node.
//if it has children, perform two pass merging on sub trees.
template <typename T>
void PairingHeap<T>::DeleteRoot(){
	if (empty()){
		delete root_;
		DecreaseNodeNumber();
		return;
	}

	DataNode<T> *old_root = root_;
	if (root_-> child == NULL){
		root_ = NULL;
		delete old_root;
	}else{
		root_ = TwoPassMerge(root_->child);
		root_->previous = NULL;
		root_->parent = NULL;
		delete old_root;
	}
	DecreaseNodeNumber();
}

//reclaim all the memory in this heap
template <typename T>
void PairingHeap<T>::ReclaimMemory(DataNode<T> *current){
	if (current != NULL){
		ReclaimMemory(current->sibling);
		ReclaimMemory(current->child);
		delete current;
	}
}

//two-pass merging, used in delete root.
template <typename T>
DataNode<T>* PairingHeap<T>::TwoPassMerge(DataNode<T> *current){
	if (current == NULL)
		return current;
	if (current->sibling == NULL){
		current->parent = NULL;
		return current;
	}

	std::vector<DataNode<T>*> tree_array(50);
	size_t num_siblings = 0;
	for (; current != NULL; num_siblings++){
		current -> parent = NULL;
		if (num_siblings == tree_array.size())
			tree_array.resize(num_siblings * 2);
		tree_array[num_siblings] = current;
		if (current -> previous != NULL)
			current -> previous -> sibling = NULL;
		current = current->sibling;
	}

	if (num_siblings == tree_array.size())
		tree_array.resize(num_siblings + 1);
	tree_array[num_siblings] = NULL;

	//forward passing (left to right)
	size_t i = 0;
	for (; i+1 < num_siblings; i += 2)
		CompareAndLink(tree_array[i], tree_array[i+1]);

	//backward passing (right to left)
	size_t j = i-2;
	if (j == num_siblings-3)
		CompareAndLink(tree_array[j], tree_array[j+2]);
	for (; j>=2; j-=2)
		CompareAndLink(tree_array[j-2], tree_array[j]);

	return tree_array[0];
}

//compare two nodes and make the smaller as the root.
template <typename T>
void PairingHeap<T>::CompareAndLink(DataNode<T>* &first, DataNode<T>* second){
	if (second == NULL)
		return;

	bool current_case = false;
	// two heaps are min ones and this heap has smaller key  OR
	// two heaps are max ones and this heap has larger key
	if ((type_of_heap_ && (*first->getKey()) <= (*second->getKey())) ||
		 (!type_of_heap_ && (*first->getKey()) >= (*second->getKey())))
		current_case = true;

	if (!current_case){
		second -> previous = first -> previous;
		second -> parent = first -> parent;
		first -> sibling = second -> child;
		if (second->child != NULL)
			second -> child -> parent = NULL;
		if (first -> sibling != NULL)
			first -> sibling -> previous = first;
		first -> previous = second;
		second -> child = first;
		first -> parent = second;
		first = second;
	}else{
		first -> previous = second -> previous;
		first -> parent = second -> parent;
		second -> sibling = first -> child;
		if (first->child != NULL)
			first -> child -> parent = NULL;
		if (second -> sibling != NULL)
			second -> sibling -> previous = second;
		second -> previous = first;
		first -> child = second;
		second -> parent = first;
	}
}

//remove corresponding node in corresponding heap.
template <typename T>
DataNode<T>* PairingHeap<T>::Remove(DataNode<T>* node){
	if (node == NULL)
		return NULL;

	//if this node has previous node, it cannot be a root.
	if (node->previous != NULL){
		if (node -> previous -> child == node){
			if (node->sibling == NULL && node->child == NULL){
				node-> previous -> child = NULL;
				node -> parent = NULL;
			}else{
				if (node->sibling != NULL && node->child == NULL){
					node->previous -> child = node->sibling;
					node->sibling -> previous = node->previous;
					node->sibling ->parent = node->parent;
				}else{
					node -> previous -> child = TwoPassMerge(node->child);
					node -> previous -> child -> parent = node -> previous;
					node -> previous -> child -> sibling = node -> sibling;
					node -> previous -> child -> previous = node -> previous;
					node -> previous -> child -> sibling -> previous = node -> previous -> child;
				}
			}
		}else{
			if (node->sibling == NULL && node->child == NULL){
				node -> previous -> sibling = NULL;
			}else{
				if (node->sibling != NULL && node->child == NULL){
					node -> previous -> sibling = node->sibling;
					node -> sibling -> previous = node -> previous;
				}else{
					node -> previous -> sibling = TwoPassMerge(node->child);
					node -> previous -> sibling -> sibling = node -> sibling;
					node -> sibling -> previous = node -> previous -> sibling;
					node -> previous -> sibling -> previous = node -> previous;
				}
			}
		}
	
		DecreaseNodeNumber(1);
		node -> previous = NULL;
		node -> child = NULL;
		node -> sibling = NULL;
		node -> parent = NULL;
		return node;
	//in this case, just delete the root.
	}else{
		DecreaseNodeNumber(1);
		DataNode<T>* return_node = root_;
		root_ = NULL;
		return return_node;
	}
}	

//set the pointer to corresponding heap.
template <typename T>
void PairingHeap<T>::setCorrespondHeap(PairingHeap<T>* correspond_heap){
	if(type_of_heap_ != correspond_heap -> type() ){
		correspond_heap_ = correspond_heap;
	}else{
		std::cout<<"setCorrespondHeap: the type of two heaps are same"<<std::endl;
	}
}

//get the corresponding heap.
template <typename T>
PairingHeap<T>* PairingHeap<T>::getCorrespondHeap(){
	return correspond_heap_;
}

//testing function: print the structure of current heap
template <typename T>
void PairingHeap<T>::PrintHeap(){
	std::queue<DataNode<T>*> keys;
	if(!empty())
		keys.push(root_);
	std::cout<<"---------------------------------------------"<<std::endl;
	std::cout<<"type: ";
	if(type_of_heap_){
		std::cout<<"min heap";
	}else{
		std::cout<<"max heap";
	}
	std::cout<<std::endl;

	std::cout<<"size: "<<number_of_nodes_<<std::endl;
	if(root_!=NULL)
		std::cout<<root_->getKey()->lambda_<<std::endl;

	while(!keys.empty()){
		DataNode<T>* current = keys.front();
		keys.pop();

		if (current->child != NULL)
			keys.push(current->child);

		if (current->sibling != NULL)
			keys.push(current->sibling);

		std::cout<<current->getKey()<<": "<<std::endl;
		std::cout<<"lambda, slope, shift: "<<current->getKey()->lambda_<<" "<<current->getKey()->slope_<<" "<<current->getKey()->shift_<<std::endl;

		std::cout<<"child: ";
		if (current->child != NULL)
			std::cout<<current->child->getKey();
		else
			std::cout<<"null";

		std::cout<<", linked node: ";
		if(current->correspond!=NULL)
			std::cout<<current->correspond->getKey();
		else
			std::cout<<"null";

		std::cout<<std::endl;

		std::cout<<"previous: ";
		if (current->previous != NULL)
			std::cout<<current->previous->getKey();
		else
			std::cout<<"null";
		std::cout<<", sibling: ";
		if (current->sibling != NULL)
			std::cout<<current->sibling->getKey();
		else
			std::cout<<"null";
		std::cout<<", parent: ";
		if (current->parent != NULL)
			std::cout<<current->parent->getKey();
		else
			std::cout<<"null";

		std::cout<<std::endl<<std::endl;
	}
}

template <typename T>
void PairingHeap<T>::setTypeOfHeap(bool type){
	type_of_heap_ = type;
}

template <typename T>
bool PairingHeap<T>::getMerged(){
	return merged_;
}

template <typename T>
void PairingHeap<T>::ClearSize(){
	root_ = NULL;
	number_of_nodes_ = 0;
}
