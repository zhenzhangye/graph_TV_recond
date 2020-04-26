/* tree node, implementation file.
 *
 * implement the functions of a tree node.
 *
 */
#include "tree_node.h"
#include "../pairing_heap/break_point.h"

// contsructor: set the parent,
// 							the index in original graph,
// 							the edge index in original graph (to parent),
// 							the index of f.
TreeNode::TreeNode(TreeNode* parent, size_t node_index, size_t edge_index, size_t f_index){
	parent_ 											= parent;
	node_index_ 									= node_index;
	edge_index_										= edge_index;
	f_index_											= f_index;
	message_ 											= NULL;
	message_node_									= NULL;
	computed_number_of_children_ 	= 0;
	number_of_children_						= 0;
}

// deconstructor: clean the children and message array.
TreeNode::~TreeNode(){
	if (number_of_children_ > 0){
		children_.clear();
	}
	if (number_of_children_ != 1){
		delete message_;
	}
}

// add a child node to current node.
// input: child_node (tree node pointer): the pointer to the node to be added as child.
void TreeNode::AddChildNode(TreeNode* child_node){
	children_.push_back(child_node);
	number_of_children_++;
}

// add a parent node to current node.
// input: parent node (tree node pointer): the pointer to the node to be added as parent.
void TreeNode::AddParentNode(TreeNode* parent_node){
	parent_ = parent_node;
}

// set the index of the edge in original graph to its parent.
void TreeNode::setEdgeIndex(size_t edge_index){
	edge_index_ = edge_index;
}

size_t TreeNode::getEdgeIndex(){
	return edge_index_;
}

size_t TreeNode::getNodeIndex(){
	return node_index_;
}

size_t TreeNode::getNumberChildren(){
	return number_of_children_;
}

// decrease the number of children by 1.
void TreeNode::DecreaseNumberChildren(){
	if (number_of_children_ == 0){
		std::cout<<"TreeNode: DecreaseNumberChildren: this node has no child node"<<std::endl;
		return;
	}
	number_of_children_ --;
}

// if current node stores message, initialize the double ended priority queue for it.
void TreeNode::InitializeDEPQ(){
		message_ = new DoubleEndedPriorityQueue<BreakPoint*>();
}

// if current node doesn't store message, set the node which stores the message.
void TreeNode::setMessageNode(TreeNode* node){
	message_node_ = node;
	message_			= node->message_;
}

TreeNode* TreeNode::getMessageNode(){
	return message_node_;
}

TreeNode* TreeNode::getFirstChild(){
	if (number_of_children_==0){
		std::cout<<"TreeNode: getFirstChild: this node has no child node"<<std::endl;
		return NULL;
	}
	return children_[0];
}

// sum over the unary and the messages from its child. Used in forward passing.
void TreeNode::SumMessageWithChildren(){
	//if current node is a leaf node or only has one child
	//do nothing else.
	
	//if current node has more than one child
	//merge all the message of children nodes
	if (number_of_children_ > 1){
		for (size_t it = 0; it < children_.size(); ++it){
			message_->Merge((children_[it])->getMessageNode()->message_);
		}
	}
}

// after sum, clip the message from both sides. Used in forward passing.
// input: a (|V|*1, array): the weight of unary.
// 				b (|V|*1, array): the unary.
// 				w (|E|*1, array): the weight of each edge.
void TreeNode::ClipMessage(double a, double b, double* w){
	//if current node is a leaf node
	//create two break points and insert them
	//update the lambda_neg lambda_pos for backward passing
	if (number_of_children_ == 0){
		lambda_neg_ = -b - w[edge_index_] / a;
		lambda_pos_ = -b + w[edge_index_] / a;
		BreakPoint*  bp1 = new BreakPoint(lambda_neg_, a, a*b+w[edge_index_]);
		message_->Insert(bp1);
		BreakPoint* bp2 = new BreakPoint(lambda_pos_, -a, w[edge_index_]-a*b);
		message_->Insert(bp2);
		//std::cout<<"After Clip: "<<std::endl;
		//PrintMessage();
	}

	// if current node has at least one child
	// perform clip from the left side
	// y = slope x + shift 
	if (number_of_children_ > 0){
		double sum_children_weight = 0;
		for (size_t it = 0; it < children_.size(); ++it){
			sum_children_weight += w[children_[it]->getEdgeIndex()];
		}

		ClipLeftSide(a, b, w, -w[edge_index_], sum_children_weight);

		//std::cout<<"After ClipLeft: "<<std::endl;
		//PrintMessage();
		ClipRightSide(a, b, w, w[edge_index_], sum_children_weight);
		//std::cout<<"After ClipRight: "<<std::endl;
		//PrintMessage();
	}
}

// clip message on the left side. used in clip message. i.e. find a break point with m(\lambda)>-weight.
// input: a (|V|*1, array): the weight in front of unary.
// 				b (|V|*1, array): the unary.
// 				w (|E|*1, array): the weight of each edge.
// 				target : - weight. where we get \lambda_neg.
// 				sum_children_weight: sum over the weights of all child nodes.
void TreeNode::ClipLeftSide(double a, double b, double* w, double target, double sum_children_weight){
	// compute the left most shift and slope
	double current_slope = a;
	double current_shift = a*b - sum_children_weight;
	
	// compute the value of the first break point.
	BreakPoint* bp = message_->getMin();
	double current_value = current_slope * bp->lambda_ + current_shift;

	double current_lambda = bp->lambda_;
	while(true){
		// if find a break point with value larger than -w, insert it back and break.
		if (current_value > target){
			break;
		}

		while(bp->lambda_ == current_lambda){
			current_slope += bp->slope_;
			current_shift += bp->shift_;
			message_->RemoveMin();
			if (message_->empty()){
				break;
			}
			bp = message_->getMin();
		}

		if (message_->empty()){
			break;
		}
		current_lambda = bp->lambda_;
		current_value = current_slope * current_lambda + current_shift;
	}

	//compute the lambda_negative.
	lambda_neg_ = (target - current_shift) / current_slope;

	// if target is not 0, insert the break point into message.
	// otherwise, we only need lambda_neg, don't update the message.
	if(target != 0)
		message_->Insert(new BreakPoint(lambda_neg_, current_slope, current_shift - target));
}

// compute the opimtal value for root node. Used after forward passing.
// input: a: the weight in front of unary for root node.
// 				b: the unary of root node.
// 				w (|E|*1, array): the weight of each edge.
void TreeNode::ComputeOptimalRootValue(double a, double b, double* w){
		double sum_children_weight = 0;
		for (size_t it = 0; it < children_.size(); ++it){
			sum_children_weight += w[children_[it]->getEdgeIndex()];
		}
		ClipLeftSide(a, b, w, 0, sum_children_weight);
}

// clip the message from right side. Find the break point with message(lambda_pos)=w
// input: a: the weight in front of unary of current node
// 				b: the unary of the root node.
// 				w (|E|*1, array): the weight of each edge.
// 				target: weight w, wehere we compute lambda_pos
// 				sum_children_weight: the sum over weights of all child nodes.
void TreeNode::ClipRightSide(double a, double b, double* w, double target, double sum_children_weight){
	//clip right side
	double current_slope = a;
	double current_shift = a*b + sum_children_weight;

	// compute the value of  the largest break point.
	BreakPoint* bp = message_->getMax();
	double current_value = current_slope * bp->lambda_ + current_shift;
	double current_lambda = bp->lambda_;
	while(true){
		//std::cout<<"right: "<<current_value<<" "<<target<<std::endl;
		//std::cout<<"here: "<<current_value<<" "<<current_slope<<" "<<current_shift<<" "<<bp->lambda_<<" "<<bp->slope_<<" "<<bp->shift_<<std::endl;
		if(current_value < target){
			break;
		}

		while(bp->lambda_ == current_lambda){
			current_slope -= bp->slope_;
			current_shift -= bp->shift_;
			message_->RemoveMax();
			if (message_->empty()){
				break;
			}
			bp = message_->getMax();
		//	std::cout<<"new message:"<<std::endl;
	//		PrintMessage();
		}

		if (message_->empty()){
			break;
		}
		current_lambda = bp->lambda_;
		current_value = current_slope * current_lambda + current_shift;
	}

	lambda_pos_ = (target - current_shift) / current_slope;
	message_->Insert(new BreakPoint(lambda_pos_, -current_slope, target - current_shift));
}

// refresh the memory of message and computed number of children for next iteration.
void TreeNode::RefreshMemory(){
	if (number_of_children_ != 1){
		message_->ClearSize();
	}
	computed_number_of_children_ = 0;
}

// clear the message heap. used in destructor.
void TreeNode::ClearHeap(){
	message_->ClearHeap();
}

// test function: print the information of current node.
void TreeNode::PrintTreeNodeInformation(){
	std::cout<<"Node index: "<<node_index_<<", ";
	std::cout<<"Edge index: "<<edge_index_<<", ";
	std::cout<<"Number of Children: "<<number_of_children_<<std::endl;
	std::cout<<"Children's node index: ";

	for(size_t it = 0; it < children_.size(); it++){
		std::cout<<children_[it]->getNodeIndex()<<", ";
	}
	std::cout<<std::endl;

	std::cout<<"Parent: ";
	if(parent_!=NULL){
		std::cout<<parent_->getNodeIndex()<<std::endl;
	}else{
		std::cout<<"Null"<<std::endl;
	}

	std::cout<<"Message node index: ";
	if(message_node_ != NULL){
		std::cout<<message_node_->getNodeIndex()<<std::endl;
	}else{
		std::cout<<"NULL"<<std::endl;
	}

	std::cout<<"lambda pos: " <<lambda_pos_<<" lambda neg: "<<lambda_neg_<<std::endl;
	std::cout<<"#######################################"<<std::endl;
}

// test function: print the passed message of current node.
void TreeNode::PrintMessage(){
	std::cout<<"Node #"<<node_index_<<std::endl;
	message_->PrintQueue();
	std::cout<<"---------------------------------"<<std::endl;
}
