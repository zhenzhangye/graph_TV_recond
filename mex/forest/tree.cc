/* tree cc file
 *
 * the definition of functions in one tree.
 *
 */

#include <vector>
#include <list>
#include <deque>

#include "tree.h"

// default constructor
Tree::Tree(){
	root_ = NULL;
}

// return the size of current tree
size_t Tree::size(){
	return size_;
}

// default deconstructor
Tree::~Tree(){
	// if it is a chain, free the memory of arrays.
	if(is_chain_){
		delete [] chain_message_;
		delete [] chain_tm_;
		delete [] chain_tp_;
		delete [] f_idx_;
		delete [] p_idx_;
	}

	// perform a bfs.
	size_t front = 0;
	size_t back = 0;
	bfs_[back] = root_;
	back = (++back)%bfs_size_;

	// clear all the nodes in this tree by bfs.
	while(front!=back){
		TreeNode* current = bfs_[front];
		front = (++front)%bfs_size_;

		if(current -> getNumberChildren() > 0){
			for(size_t it = 0; it < current->children_.size(); ++it) {
				bfs_[back] = current->children_[it];
				back = (++back)%bfs_size_;
			}
		}

		delete current;
	}
	delete bfs_;
	delete visit_leaves_queue_;
}

// create the tree structure.
// input: graph(pointer): the original graph.
// 				source_node(array): the indexes of the source nodes for edges.
// 				target_node(array): the indexes of the target nodes for edges.
// 				weights: the weight of each edge.
// 				number_forests: the forests created until now.
// 				node_index: the index of root node.
// 				number_vertice_tree: the number of vertices in this tree.
// 				use_chain: ture: if treat chain structure with efficient backward solver.
Tree::Tree(Graph* graph, size_t* source_node, size_t* target_node, 
		WeightMap& weights, size_t number_forests, size_t node_index, size_t number_vertice_tree, bool use_chain){

	// the number of vertices in current tree.
	size_ = number_vertice_tree + 1;

	// initialize the root node in the tree.
	root_ = new TreeNode(NULL, node_index, 0, number_forests);
	
	// assume current tree is a chain.
	is_chain_ = use_chain;

	//create the vector which record the node in the tree and the vertice in the graph.
	size_t number_vertice_graph = boost::num_vertices(*graph);
	std::vector<TreeNode*> link_graph_to_tree(number_vertice_graph, NULL);
	link_graph_to_tree[node_index] = root_;

	TreeNode* new_parent_node;
	TreeNode* new_child_node;

	//transverse all the nodes in current minimal spanning tree.
	//convert it into custom structure.
	for(size_t prim_index = 0; prim_index<number_vertice_tree; ++prim_index){
		double current_weight = boost::get(weights, 
																				boost::edge(source_node[prim_index], 
																										target_node[prim_index], 
																										*graph).first
																			 ).inactive_weight_;
		
		//record the edge_index of these two nodes.
		size_t edge_index = boost::get(weights,
																	 boost::edge(source_node[prim_index],
																							 target_node[prim_index],
																							 *graph).first
																	).original_index_;

		//check if the parent is in the tree.
		if (link_graph_to_tree[source_node[prim_index]] == NULL){
			// in this case, create a new tree node.
			new_parent_node = new TreeNode(NULL, source_node[prim_index], 0, number_forests);

			// link the tree node to the graph vertice
			link_graph_to_tree[source_node[prim_index]] = new_parent_node;

			// check if the child is in the tree
			if (link_graph_to_tree[target_node[prim_index]] == NULL){
				// in this case, create a new tree node.
				new_child_node = new TreeNode(new_parent_node, target_node[prim_index], edge_index, number_forests);

				// link the tree node to the graph vertice
				link_graph_to_tree[target_node[prim_index]] = new_child_node;

				// add the child node to the parent node
				new_parent_node->AddChildNode(new_child_node);
			}else{
				// in this case, no need to create a new tree node.
				// add the child node to the parent node.
				new_parent_node->AddChildNode(link_graph_to_tree[target_node[prim_index]]);
				
				//add parent to this child node
				link_graph_to_tree[target_node[prim_index]] -> AddParentNode(new_parent_node);

				// set the edge index in the child node.
				link_graph_to_tree[target_node[prim_index]]->setEdgeIndex(edge_index);
			}	
		}else{
			// in this case parent node is in the tree.
			new_parent_node = link_graph_to_tree[source_node[prim_index]];

			// check if current tree is a chain.
			// if the root node has more than 3 children, this tree cannot be a chain.
			if (new_parent_node == root_ && new_parent_node->getNumberChildren() >= 2){
				is_chain_ = false;
			}
			// if other node has more than 2 children, this tree cannot be a chain.
			if (new_parent_node != root_ && new_parent_node->getNumberChildren() >= 1){
				is_chain_ = false;
			}

			//check if the child is in the tree
			if (link_graph_to_tree[target_node[prim_index]] == NULL){
				// in this case, create a new tree node.
				new_child_node = new TreeNode(new_parent_node, target_node[prim_index], edge_index, number_forests);

				// link the tree node to the graph vertice
				link_graph_to_tree[target_node[prim_index]] = new_child_node;

				// add the child node to the parent node
				new_parent_node->AddChildNode(new_child_node);
			}else{
				// in this case, add the child node to the parent node.
				new_parent_node->AddChildNode(link_graph_to_tree[target_node[prim_index]]);
				
				//add parent to this child node
				link_graph_to_tree[target_node[prim_index]] -> AddParentNode(new_parent_node);

				//set the edge index in the child node.
				link_graph_to_tree[target_node[prim_index]]->setEdgeIndex(edge_index);
			}
		}
	}

	// store all the leaf nodes of this tree
	bfs_size_ = size_*4;
	bfs_ = new TreeNode*[bfs_size_];
	SearchLeavesInTree();
	queue_length_ = leaves_.size()*4;
	visit_leaves_queue_ = new TreeNode*[queue_length_];

	// if the tree is a chain and use 1D solver,
	// reconstruct this tree (root might not be the end of the chain).
	if (is_chain_ && use_chain){
		ReconstructTree();
	}else{
		//if it is not a chain, allocate the memory space for each 
		AllocateQueueMemory();
	}
}

// create the tree structure.
// input: graph(pointer): the original graph.
// 				remaining_graph(pointer): remove the edges from this graph after insertion.
// 				source_node(array): the indexes of the source nodes for edges.
// 				target_node(array): the indexes of the target nodes for edges.
// 				weights: the weight of each edge.
// 				number_forests: the forests created until now.
// 				node_index: the index of root node.
// 				number_vertice_tree: the number of vertices in this tree.
// 				use_chain: ture: if treat chain structure with efficient backward solver.
Tree::Tree(Graph* graph, Graph* remaining_graph, size_t* source_node, size_t* target_node, 
		WeightMap& weights, size_t number_forests, size_t node_index, size_t number_vertice_tree, bool use_chain){

	// the number of vertices in current tree.
	size_ = number_vertice_tree + 1;
	//initialize the root node in the tree.
	root_ = new TreeNode(NULL, node_index, 0, number_forests);
	is_chain_ = use_chain;

	//create the vector which record the node in the tree and the vertice in the graph.
	size_t number_vertice_graph = boost::num_vertices(*graph);
	std::vector<TreeNode*> link_graph_to_tree(number_vertice_graph, NULL);
	link_graph_to_tree[node_index] = root_;

	TreeNode* new_parent_node;
	TreeNode* new_child_node;

	//transverse all the nodes in current minimal spanning tree.
	//convert it into custom structure.
	for(size_t prim_index = 0; prim_index<number_vertice_tree; ++prim_index){
		double current_weight = boost::get(weights, 
																				boost::edge(source_node[prim_index], 
																										target_node[prim_index], 
																										*graph).first
																			 ).inactive_weight_;
		
		//remove this edge from graph. Insert it into custom tree.
		boost::remove_edge(source_node[prim_index], target_node[prim_index], *remaining_graph);

		//record the edge_index of these two nodes.
		size_t edge_index = boost::get(weights,
																	 boost::edge(source_node[prim_index],
																							 target_node[prim_index],
																							 *graph).first
																	).original_index_;

		//check if the parent is in the tree.
		if (link_graph_to_tree[source_node[prim_index]] == NULL){
			// in this case, create a new tree node.
			new_parent_node = new TreeNode(NULL, source_node[prim_index], 0, number_forests);

			// link the tree node to the graph vertice
			link_graph_to_tree[source_node[prim_index]] = new_parent_node;

			// check if the child is in the tree
			if (link_graph_to_tree[target_node[prim_index]] == NULL){
				// in this case, create a new tree node.
				new_child_node = new TreeNode(new_parent_node, target_node[prim_index], edge_index, number_forests);

				// link the tree node to the graph vertice
				link_graph_to_tree[target_node[prim_index]] = new_child_node;

				// add the child node to the parent node
				new_parent_node->AddChildNode(new_child_node);
			}else{
				// in this case, no need to create a new tree node.
				// add the child node to the parent node.
				new_parent_node->AddChildNode(link_graph_to_tree[target_node[prim_index]]);
				
				//add parent to this child node
				link_graph_to_tree[target_node[prim_index]] -> AddParentNode(new_parent_node);

				// set the edge index in the child node.
				link_graph_to_tree[target_node[prim_index]]->setEdgeIndex(edge_index);
			}	
		}else{
			// in this case parent node is in the tree.
			new_parent_node = link_graph_to_tree[source_node[prim_index]];

			// check if current tree is a chain.
			// if the root node has more than 3 children, this tree cannot be a chain.
			if (new_parent_node == root_ && new_parent_node->getNumberChildren() >= 2){
				is_chain_ = false;
			}
			// if other node has more than 2 children, this tree cannot be a chain.
			if (new_parent_node != root_ && new_parent_node->getNumberChildren() >= 1){
				is_chain_ = false;
			}

			//check if the child is in the tree
			if (link_graph_to_tree[target_node[prim_index]] == NULL){
				// in this case, create a new tree node.
				new_child_node = new TreeNode(new_parent_node, target_node[prim_index], edge_index, number_forests);

				// link the tree node to the graph vertice
				link_graph_to_tree[target_node[prim_index]] = new_child_node;

				// add the child node to the parent node
				new_parent_node->AddChildNode(new_child_node);
			}else{
				// in this case, add the child node to the parent node.
				new_parent_node->AddChildNode(link_graph_to_tree[target_node[prim_index]]);
				
				//add parent to this child node
				link_graph_to_tree[target_node[prim_index]] -> AddParentNode(new_parent_node);

				//set the edge index in the child node.
				link_graph_to_tree[target_node[prim_index]]->setEdgeIndex(edge_index);
			}
		}
	}

	// store all the leaf nodes of this tree
	bfs_size_ = size_*4;
	bfs_ = new TreeNode*[bfs_size_];
	SearchLeavesInTree();
	queue_length_ = leaves_.size()*4;
	visit_leaves_queue_ = new TreeNode*[queue_length_];

	// if the tree is a chain and use 1D solver,
	// reconstruct this tree (root might not be the end of the chain).
	if (is_chain_ && use_chain){
		ReconstructTree();
	}else{
		//if it is not a chain, allocate the memory space for each 
		AllocateQueueMemory();
	}
}

// search all the leaf nodes in current tree. i.e. the nodes with 0 children.
void Tree::SearchLeavesInTree(){
	//start from root of the tree, perform BFS
	size_t bfs_front = 0;
	size_t bfs_back = 0;
	bfs_[bfs_back] = root_;
	bfs_back = (++bfs_back) % bfs_size_;
	TreeNode* current;

	while(bfs_front!=bfs_back){
		current = bfs_[bfs_front];
		bfs_front = (++bfs_front) % bfs_size_;

		// if current node has child node, push them into bfs queue
		// otherwise, push this node into leaves_ queue.
		if(current->getNumberChildren()>0){
			for(size_t it =0; it < current->children_.size(); ++it){
				bfs_[bfs_back] = current->children_[it];
				bfs_back = (++bfs_back) % bfs_size_;
			}
		}else{
			leaves_.push_back(current);
		}
	}
}

// if current tree is a chain, find the new root (end point). (because current root might make it as a tree).
void Tree::ReconstructTree(){
	// create the array for chain solver.
	chain_message_ = new double[16*size_];
	chain_tm_      = new double[4*size_];
	chain_tp_			 = new double[4*size_];
	TreeNode* current;

	// if the root_ has one child, it must be the end of this chain.
	// Nothing needs to do.
	if(root_->getNumberChildren()>1){

		// if not, we start from a leaf node
		current = leaves_[0];

		// add the parent as a child of the leaf node.
		TreeNode* parent = current->parent_;
		current->AddChildNode(parent);
		current->parent_ = NULL;

		//modify the edge index
		size_t old_edge_index = current->getEdgeIndex();
		size_t temp_edge_index;
		current->setEdgeIndex(0);

		//modify all intermediate nodes
		while(parent != root_){
			current = parent;
			parent = current->parent_;

			current->parent_ = current->getFirstChild();
			current->children_.erase(current->children_.begin());
			current->children_.push_back(parent);

			temp_edge_index = current->getEdgeIndex();
			current->setEdgeIndex(old_edge_index);
			old_edge_index = temp_edge_index;
		}

		//modify the root node.
		if (root_->getFirstChild()==current){
			root_->children_.erase(root_->children_.begin());
		}else{
			root_->children_.pop_back();
		}

		root_->DecreaseNumberChildren();

		root_->parent_ = current;
		root_->setEdgeIndex(old_edge_index);

		root_ = leaves_[0];
		leaves_.pop_back();
	}

	// save the arrays used for chain solver.
	f_idx_ = new size_t[size_];
	p_idx_ = new size_t[size_-1];
	current = root_;
	f_idx_[0] = current->getNodeIndex();
	for(size_t i = 1; i<size_; ++i){
		current = current->getFirstChild();
		f_idx_[i] = current->getNodeIndex();
		p_idx_[i-1] = current->getEdgeIndex();
	}

}

// allocate the queue memory for leaf/branch nodes. Used for storing messages in backward solver.
void Tree::AllocateQueueMemory(){
	// copy the leaf nodes.
	std::copy(leaves_.begin(), leaves_.end(), visit_leaves_queue_);
	TreeNode* current;

	// perform a bfs.
	size_t front = 0;
	size_t back = leaves_.size();
	while(front!=back){
		current = visit_leaves_queue_[front];
		front = (++front) % (queue_length_);

		// if current node has child nodes more than 1, initialize the space to store message
		if (current->getNumberChildren()!=1){
			current->InitializeDEPQ();
			current->setMessageNode(current);
		}else{
			// otherwise, the message will be stored in the first child node.
			current->setMessageNode(current->getFirstChild()->getMessageNode());
		}

		if (current == root_){
			break;
		}

		// push parent nodes into the queue.
		if (++current->parent_->computed_number_of_children_ == current->parent_->getNumberChildren()){
			current->parent_->computed_number_of_children_ = 0;
			visit_leaves_queue_[back] = current->parent_;
			back = (++back) % (queue_length_);
		}
	}
}

// the backward solver on this tree.
// solving: v = argmin a/2 ||v + b||^2 + ||Kv||_1.
double Tree::BackwardSolver(double* a, double* b, double* w, double* v){
	ForwardPassing(a, b, w);
	BackwardPassing(a, b, w, v);
	return 0;
}

// the forward passing phrase in backward solver. From leaves to root.
// input: a (|V|*1, array): the weight in front of unary.
// 				b (|V|*1, array): the unary.
// 				w (|E|*1, array): the weight of each edge.
double Tree::ForwardPassing(double* a, double* b, double* w){
	// copy the leaf nodes in current tree to current queue.
	std::copy(leaves_.begin(), leaves_.end(), visit_leaves_queue_);
	TreeNode* current;
	size_t front = 0;
	size_t back = leaves_.size();

	// perform a bfs from leaf to root
	while(1){
		current = visit_leaves_queue_[front];
		front = (++front)% queue_length_;

		//insert the parent node into the list,
		//when all children node is computed.
		if(current!=root_){
			++current->parent_->computed_number_of_children_;
			if(current->parent_->computed_number_of_children_ == current->parent_->getNumberChildren()){
				current->parent_->computed_number_of_children_ = 0;
				visit_leaves_queue_[back] = current->parent_;
				back = (++back)% queue_length_;
			}
		}

		// summing over the messages from child nodes.
		current->SumMessageWithChildren();
		//std::cout<<"After Sum:"<<std::endl;
		//current->PrintMessage();

		if(current==root_){
			break;
		}

		//clip message
		current->ClipMessage(a[current->getNodeIndex()], b[current->getNodeIndex()], w);
		//current->PrintTreeNodeInformation();

	}
	return 0;
}

// the backward passing phrase in backward solver. from root to leaves.
// input: a (|V|*1, array): the weight in front of unary.
// 				b (|V|*1, array): the unary.
// 				w (|E|*1, array): the weight of each edge.
// 				v (|V|*1, array): the array to store result.
double Tree::BackwardPassing(double* a, double* b, double *w, double *v){
	// compute the optimal value of root node.
	root_->ComputeOptimalRootValue(a[root_->getNodeIndex()], b[root_->getNodeIndex()], w);
	v[root_->getNodeIndex()] = root_->lambda_neg_;

	// perform a bfs from root to leaves. Insert child nodes of root.
	size_t front = 0;
	size_t back = 0;
	for (size_t it = 0; it < root_->children_.size(); ++it){
		bfs_[back] = root_->children_[it];
		back = (++back) % bfs_size_;
	}

	TreeNode* current;
	size_t parent_node_index;
	size_t current_node_index;
	while(front!=back){
		current = bfs_[front];
		front = (++front) % bfs_size_;

		// add child nodes into bfs queue.
		for (size_t it = 0; it < current->children_.size(); ++it){
			bfs_[back] = current->children_[it];
			back = (++back) % bfs_size_;
		}

		// v_{i-1} >= lambda+ : v_i = lambda+
		// v_{i-1} <= lambda- : v_i = lambda-
		// otherwise: v_i = v_{i-1}
		parent_node_index = current->parent_->getNodeIndex();
		current_node_index = current->getNodeIndex();
		if (v[parent_node_index] >= current->lambda_pos_){
			v[current_node_index] = current->lambda_pos_;
		}else if (v[parent_node_index] <= current->lambda_neg_){
			v[current_node_index] = current->lambda_neg_;
		}else{
			v[current_node_index] = v[parent_node_index];
		}
	}
	return 0;
}

// solve the linear system to recover p^{k+1} from v (dual solution of dual update). nabla^T p = u+f
// input: b (|V|*1, array): the unary.
// 				w (|E|*1, array): the weight of each edge
// 				v (|V|*1, array): the result of dual problem of dual update.
// 				r (|E|*1, array): p^{k+1}. the real dual solution.
void Tree::SolveLinearSystem(double* b, double* w, double* v, double* r){
	// perform a bfs on current tree.
	size_t front = 0;
	size_t back = 0;
	bfs_[back] = root_;
	back = (++back) % bfs_size_;

	// start from root to leaves, compute the v.
	TreeNode* current;
	size_t it;
	while(front!=back){
		current = bfs_[front];
		front = (++front) % bfs_size_;

		if(current->getNumberChildren() > 0){
			for (it = 0; it < current->children_.size(); ++it){
				bfs_[back] = current->children_[it];
				back = (++back) % bfs_size_;
			}
		}

		v[current->getNodeIndex()] += b[current->getNodeIndex()];
	}

	// start from the leaves to root, compute the result r.
	double child_value;
	std::copy(leaves_.begin(), leaves_.end(), visit_leaves_queue_);
	front = 0;
	back = leaves_.size();
	while(front!=back){
		current = visit_leaves_queue_[front];
		front = (++front) % queue_length_;

		if(current == root_){
			break;
		}

		child_value = 0;
		if(current->getNumberChildren() > 0){
			for (it = 0; it < current->children_.size(); ++it){
				if (current->getNodeIndex() > (current->children_[it])->getNodeIndex())
					child_value -= r[(current->children_[it])->getEdgeIndex()];
				else
					child_value += r[(current->children_[it])->getEdgeIndex()];

			}
		}

		if(current->getNodeIndex() > current->parent_->getNodeIndex()){
			r[current->getEdgeIndex()] = child_value + v[current->getNodeIndex()];
		}else{
			r[current->getEdgeIndex()] = -(child_value + v[current->getNodeIndex()]);
		}

		if(++current->parent_->computed_number_of_children_ == current->parent_->getNumberChildren()){
			visit_leaves_queue_[back] = current->parent_;
			back = (++back) % queue_length_;
			current->parent_->computed_number_of_children_ = 0;
		}
	}

	// get rid of the weight to compute the final result r.
	front = 0;
	back = 0;
	bfs_[back] = root_;
	back = (++back) % bfs_size_;
	while(front!=back){
		current = bfs_[front];
		front = (++front) % bfs_size_;

		if(current->getNumberChildren() > 0){
			for (it = 0; it < current->children_.size(); ++it){
				bfs_[back] = current->children_[it];
				back = (++back) % bfs_size_;
			}
		}
		if(current==root_)
			continue;

		r[current->getEdgeIndex()] /= w[current->getEdgeIndex()];
	}
}

// refresh the memory for next iteration.
void Tree::RefreshMemory(){
	size_t front = 0;
	size_t back = 0;
	bfs_[back] = root_;
	back = (++back) % bfs_size_;
	// clear the message heap for root node.
	root_ -> ClearHeap();

	// perform a bfs.
	// if a node has place to store message, clean it.
	TreeNode* current;
	while(front!=back){
		current = bfs_[front];
		front = (++front)%bfs_size_;

		if(current->getNumberChildren() > 0){
			for (size_t it = 0; it < current->children_.size(); ++it){
				bfs_[back] = current->children_[it];
				back = (++back) % bfs_size_;
			}
		}

		current->RefreshMemory();
	}
}

// compute the f for backward solver. f = K^T*p_l^k + u_bar/t
// input: u_bar (|V|*1, array): the over-relaxation of primal variable.
// 				pk (|E|*1, array): the dual variable at iteration k.
// 				w (|E|*1, array): the weight of each edge.
// 				t: the weight in front of regularization.
// 				b: the result f.
void Tree::ComputeF(double* u_bar, double* pk, double* w, double t, double* b){
	// starting from the root.
	size_t node_index = root_->getNodeIndex();
	b[node_index] =  u_bar[node_index] / t;

	// perform a bfs on current tree. compute the b for each node.
	size_t it;
	size_t front = 0;
	size_t back = 0;
	for (it = 0; it < root_->children_.size(); ++it){
		if(node_index < (root_->children_[it])->getNodeIndex())
			b[node_index] -= w[(root_->children_[it])->getEdgeIndex()] * pk[(root_->children_[it])->getEdgeIndex()];
		else
			b[node_index] += w[(root_->children_[it])->getEdgeIndex()] * pk[(root_->children_[it])->getEdgeIndex()];
		bfs_[back] = root_->children_[it];
		back = (++back) % bfs_size_;
	}

	TreeNode* current;
	while(front!=back){
		current = bfs_[front];
		front = (++front)%bfs_size_;

		node_index = current->getNodeIndex();
		b[node_index] = u_bar[node_index]/t;

		if(node_index > current->parent_->getNodeIndex())
			b[node_index] += w[current->getEdgeIndex()] * pk[current->getEdgeIndex()];
		else
			b[node_index] -= w[current->getEdgeIndex()] * pk[current->getEdgeIndex()];

		if(current->getNumberChildren()>0){
			for (it = 0; it < current->children_.size(); ++it){
				if(node_index < (current->children_[it])->getNodeIndex())
					b[node_index] -= w[(current->children_[it])->getEdgeIndex()] * pk[(current->children_[it])->getEdgeIndex()];
				else
					b[node_index] += w[(current->children_[it])->getEdgeIndex()] * pk[(current->children_[it])->getEdgeIndex()];
				bfs_[back] = current->children_[it];
				back = (++back) % bfs_size_;
			}
		}
	}
}

// test function: print the structure of current tree. root, leaves and so on.
void Tree::PrintTreeStructure(){
	size_t front = 0;
	size_t back = 0;
	bfs_[back] = root_;
	back = (++back) % bfs_size_;

	TreeNode* current;
	while(front!=back){
		current = bfs_[front];
		front = (++front) % bfs_size_;

		if(current->getNumberChildren()>0){
			for(size_t it = 0; it < current->children_.size(); ++it){
				bfs_[back] = current->children_[it];
				back = (++back) % bfs_size_;
			}
		}

		current->PrintTreeNodeInformation();
	}

	std::cout<<"leaves in this tree: ";
	std::copy(leaves_.begin(), leaves_.end(), visit_leaves_queue_);
	front = 0;
	back = leaves_.size();
	while (front!=back){
		std::cout<<visit_leaves_queue_[front]->getNodeIndex()<<", ";
		front = (++front) % queue_length_;
	}
	std::cout<<std::endl;

	if(is_chain_){
		std::cout<<"f indexes: ";
		for(size_t i = 0; i<size_; ++i){
			std::cout<<f_idx_[i]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"p indexes: ";
		for(size_t i = 0; i<size_-1; ++i){
			std::cout<<p_idx_[i]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"-----------------------------"<<std::endl;

}

// compute the unary in dual problem for chain. f = K^T * p_l^k + u_bar/t
// input: u_bar (|V|*1, array): the over-relaxation of primal variable.
// 			  pk (|E|*1, array): the dual variable at iteration k.
// 			  w (|E|*1, array): the weight of each edge.
// 			  t: the weight in front of regularization.
// 			  b (|V|*1, array): the unary, final result.
// 			  chain_weight (|E|*1, array): the weight for chains.
void Tree::ComputeFChain(double* u_bar, double* pk, double* w, double t, double* b, double* chain_weight){
	if(!is_chain_){
		std::cout<<"Tree: ComputeFChain: this tree is not a chain"<<std::endl;
	}
	// iterate over the chain, compute the unary based on previous node.
	if (f_idx_[0]<f_idx_[1])
		b[0] = -w[p_idx_[0]] * pk[p_idx_[0]] + u_bar[f_idx_[0]] / t;
	else
		b[0] = w[p_idx_[0]] * pk[p_idx_[0]] + u_bar[f_idx_[0]] / t;

	size_t prev_p_ind;
	size_t p_ind;
	size_t f_ind;
	for (size_t i = 1; i<size_-1; ++i){
		prev_p_ind = p_idx_[i-1];
		p_ind = p_idx_[i];
		f_ind = f_idx_[i];

		chain_weight[i-1] = w[prev_p_ind];

		if(f_ind < f_idx_[i+1])
			b[i] = -w[p_ind]*pk[p_ind] + u_bar[f_ind] /t;
		else
			b[i] = w[p_ind]*pk[p_ind] + u_bar[f_ind] /t;

		if(f_ind > f_idx_[i-1])
			b[i] += w[prev_p_ind] * pk[prev_p_ind];
		else
			b[i] -= w[prev_p_ind] * pk[prev_p_ind];
	}

	// compute the unary of the end node.
	p_ind = p_idx_[size_-2];
	f_ind = f_idx_[size_-1];
	chain_weight[size_-2] = w[p_ind];
	
	if(f_ind > f_idx_[size_-2])
		b[size_-1] = w[p_ind]*pk[p_ind] + u_bar[f_ind]/t;
	else
		b[size_-1] = -w[p_ind]*pk[p_ind] + u_bar[f_ind]/t;
}

// The backward solver for chain, much efficient.
// Code is from Vladimir Kolmogorov, see https://pub.ist.ac.at/~vnk/papers/TV_TREE.html.
void Tree::TV1DWeightedKPR(double* ya, double* yb, double* lambda, double* beta, size_t n){
	double* message = chain_message_;
	double* tm = chain_tm_;
	double* tp = chain_tp_;
	double mean_slope, mean_offset;
	int l, r;

	l = 3 * n + 1;
	r = 3 * n + 4;

	mean_slope = ya[0];
	mean_offset = yb[0];

	// message after first iteration
	message[l-2] = -mean_slope; // slope 1
	message[l-1] = -lambda[0]-yb[0]; // intercept 1
	message[l] = (-lambda[0] - yb[0]) / ya[0]; // BP 1
	message[l+1] = 0; // slope 2
	message[l+2] = 0; // intercept 2
	message[l+3] = ( lambda[0] - yb[0]) / ya[0]; // BP 2 
	message[l+4] = -mean_slope; // slope 3
	message[l+5] = lambda[0]-yb[0]; // intercept 3

	// lambda^- and lambda^+ in the paper
	tm[0] = (-lambda[0] - yb[0]) / ya[0];
	tp[0] = ( lambda[0] - yb[0]) / ya[0];

	int lo, hi;
	bool found;
	for(int k=1; k<n-1; ++k) {

		// add unary term
		mean_slope += ya[k];
		mean_offset += yb[k];

		// do clipping

		// find breakpoint where segment is larger -lambda[k]
		found = false;
		for(lo = l; lo <= r; lo += 3) {
			// compute breakpoint
			const double pt = (-lambda[k] - (message[lo - 1] + mean_offset)) / (message[lo - 2] + mean_slope);

			if(pt < message[lo]) {
				message[lo - 3] = pt;
				message[lo - 4] = -mean_offset - lambda[k];
				message[lo - 5] = -mean_slope;
				found = true;
				break;
			}
		}

		if(!found) {
			l = r;
			const double pt = (-lambda[k] - (message[l+2] + mean_offset)) / (message[l+1] + mean_slope);

			message[l] = pt;
			message[l-1] = -mean_offset - lambda[k];
			message[l-2] = -mean_slope;

			tm[k] = message[l];
		}
		else {
			l = lo - 3;
			tm[k] = message[l];
		}

		// find breakpoint where segment is smaller than lambda[k]
		found = false;
		for(hi = r; hi >= l; hi -= 3) {
			// compute breakpoint
			const double pt = (lambda[k] - (message[hi + 2] + mean_offset)) / (message[hi + 1] + mean_slope);

			if(pt > message[hi]) {
				message[hi+3] = pt;
				message[hi+4] = -mean_slope;
				message[hi+5] = -mean_offset + lambda[k];
				found = true;
				break;
			}
		}

		// something went wrong...
		/*
			 if(!found) {
			 return false;
			 }
			 */

		r = hi + 3;
		tp[k] = message[r];
	}

	mean_slope += ya[n-1];
	mean_offset += yb[n-1];

	// compute beta[n-1]
	double a, b;
	for(lo = l-3; lo <= r; lo += 3) {
		a = message[lo+1]+mean_slope;
		b = message[lo+2]+mean_offset;

		if(a * message[lo+3] + b > 0) {
			break;
		}
	}
	beta[n-1] = -b / a;

	// Compute the rest of the coefficients, by the
	// back-pointers
	for (int k=n-2; k>=0; k--) {
		if (beta[k+1]>tp[k]) beta[k] = tp[k];
		else if (beta[k+1]<tm[k]) beta[k] = tm[k];
		else beta[k] = beta[k+1];
	}

	/*
		 delete[] tm;
		 delete[] tp;
		 delete[] message;
		 */

}

// if current tree is a chain, solve the linear system on it to recover p^{k+1}. nalba_w^T p = u+f.
// input: b (|V|*1, array): unary (f in above equation).
// 			  w (|E|*1, array): the weight of each edge.
// 			  v (|V|*1, array): the result of dual problem for dual update.
// 			  r (|E|*1, array): the final result p^{k+1}.
void Tree::SolveLinearSystemChain(double* b, double* w, double* v, double* r){
	// iterate over this chain, compute the final r.
	size_t p_ind = p_idx_[size_-2];
	size_t f_ind = f_idx_[size_-1];
	r[p_ind] = (v[size_-1] + b[size_-1])/w[p_ind];

	double prev_p = r[p_ind] * w[p_ind];
	
	if (f_ind < f_idx_[size_-2])
		r[p_ind] = -r[p_ind];

	for(size_t i = size_-2; i>=1; i--){
		p_ind = p_idx_[i-1];
		r[p_ind] = (v[i] + b[i] + prev_p) / w[p_ind];
		prev_p = r[p_ind] * w[p_ind];
		if (f_idx_[i] < f_idx_[i-1])
			r[p_ind] = -r[p_ind];
	}
}
