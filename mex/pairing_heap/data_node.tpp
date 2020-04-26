//default constructor
template<typename T>
DataNode<T>::DataNode(){
	this->correspond 	= NULL;
	this->previous 		= NULL;
	this->parent 			= NULL;
	this->sibling 		= NULL;
	this->child 			= NULL;
	this->key 				= NULL;
}

//constructor with certain key
template<typename T>
DataNode<T>::DataNode(T key){
	this->key					= key;
	this->correspond 	= NULL;
	this->previous 		= NULL;
	this->parent 			= NULL;
	this->sibling 		= NULL;
	this->child 			= NULL;
}

//default destructor
template<typename T>
DataNode<T>::~DataNode(){
	delete this->key;
	this->correspond 	= NULL;
	this->previous 		= NULL;
	this->parent 			= NULL;
	this->sibling 		= NULL;
	this->child 			= NULL;
}

//get previous node
template<typename T>
DataNode<T>* DataNode<T>::getPrevious(){
	return this->previous;
}

//get sibling node
template<typename T>
DataNode<T>* DataNode<T>::getSibling(){
	return this->sibling;
}

//get child node
template<typename T>
DataNode<T>* DataNode<T>::getChild(){
	return this->child;
}

//get parent node
template<typename T>
DataNode<T>* DataNode<T>::getParent(){
	DataNode<T>* current = this;
	while(current->parent==NULL){
		current = current->getPrevious();
		if (current == NULL){
			return NULL;
		}
	}
	return current->parent;
}

//get correspond node in another heap
template<typename T>
DataNode<T>* DataNode<T>::getCorrespond(){
	return this->correspond;
}

//get the key of this node
template<typename T>
T DataNode<T>::getKey(){
	return this->key;
}

//get the type of this node, true: data node
template<typename T>
bool DataNode<T>::type(){
	return this->child==NULL;
}

//set the correspond node in another heap
template<typename T>
void DataNode<T>::setCorrespond(DataNode<T>* node){
	//std::cout<<(node==NULL)<<std::endl;
	if (node == NULL)
		this->correspond = NULL;
	else
		this->correspond = node;
}

template<typename T>
bool DataNode<T>::operator<(DataNode<T>& rhs){
	return this->getKey() < rhs.getKey();
}

template<typename T>
bool DataNode<T>::operator>(DataNode<T>& rhs){
	return this->getKey() > rhs.getKey();
}

template<typename T>
bool DataNode<T>::operator<=(DataNode<T>& rhs){
	return this->getKey() <= rhs.getKey();
}

template<typename T>
bool DataNode<T>::operator>=(DataNode<T>& rhs){
	return this->getKey() >= rhs.getKey();
}

template<typename T>
bool DataNode<T>::operator==(DataNode<T>& rhs){
	return this->getKey() == rhs.getKey();
}
