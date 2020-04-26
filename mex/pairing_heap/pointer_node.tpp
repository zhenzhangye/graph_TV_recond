//default constructor
template<typename T>
PointerNode<T>::PointerNode(){
	this->correspond 	= NULL;
	this->previous 		= NULL;
	this->sibling 		= NULL;
	this->child 			= NULL;
}

//constructor with certain key
template<typename T>
PointerNode<T>::PointerNode(DataNode<T>* &node){
	this->correspond 	= node;
	node->setCorrespond(this);
	this->previous 		= NULL;
	this->sibling 		= NULL;
	this->child 			= NULL;
}

//default deconstructor
template<typename T>
PointerNode<T>::~PointerNode(){
}

//get previous node
template<typename T>
Node<T>* PointerNode<T>::getPrevious(){
	return this->previous;
}

//get sibling node
template<typename T>
Node<T>* PointerNode<T>::getSibling(){
	return this->sibling;
}

//get child node
template<typename T>
Node<T>* PointerNode<T>::getChild(){
	return this->child;
}

//get correspond node in another heap
template<typename T>
Node<T>* PointerNode<T>::getCorrespond(){
	return this->correspond;
}

//get the key of this node
template<typename T>
T PointerNode<T>::getKey(){
	return this->correspond->getKey();
}

//get the type of this node, false: pointer node
template<typename T>
bool PointerNode<T>::getType(){
	return false;
}
