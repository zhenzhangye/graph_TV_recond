/* tree header file
 *
 * define the functions and variables for one tree
 *
 *
 */
#ifndef __TREE_H__
#define __TREE_H__

#include <vector>

#include "tree_node.h"
#include "../typedefs.h"

class Tree{
public:
	//contructor and deconstructor
	Tree();

	// the tree constructor, given a partitioned tree, create the corresponding tree data structure.
	// each edge is removed from original graph after insertion into the tree.
	// here we consider the remaining edge in the graph
	Tree(Graph* graph, Graph* remaining_graph, size_t* source_node, size_t* target_node, 
		WeightMap& weights, size_t number_forests, size_t node_index, size_t number_vertice_tree, bool use_chain);

	// the tree constructor, given a partitioned tree, create the corresponding tree data structure.
	// each edge is removed from original graph after insertion into the tree.
	Tree(Graph* graph, size_t* source_node, size_t* target_node, WeightMap& weights, 
		size_t number_forests, size_t node_index, size_t number_vertice_tree, bool use_chain);
	
	// delete all the tree nodes,
	// the array message if this tree is a chain,
	// leaves queue and otheres
	~Tree();

	//public member function
	// the size of this tree (number of vertices)
	size_t				size();

	// used to compute the f in dual problem for dual update. f = K^T*p_l^k + u_bar /t
	void					ComputeF(double* u_bar, double* pk, double* w, double t, double* b);

	// used to compute the f in dual problem for dual update. This f is for chain solver. f = K^T*p_l^k + u_bar / t.
	void					ComputeFChain(double* u_bar, double* pk, double* w, double t, double* b, double* chain_weight);

	// the backward solver for current tree. u = argmin a/2 *(u+f)^2 + ||Ku||_1.
	double				BackwardSolver(double* a, double* b, double* w, double* v);

	// the fast backward solver for chain. u = argmin a/2 * (u+f)^2 + ||Ku||_1
	void					TV1DWeightedKPR(double* ya, double* yb, double* lambda, double* beta, size_t n);

	// using the dual result to recover the solution of original p^{k+1}. nabla_w^T P = u+f
	void 					SolveLinearSystem(double* b, double* w, double* v, double* r);

	// using the dual result to recover the solution of original p^{k+1} in chain. nabla_w^T P = u+f
	void					SolveLinearSystemChain(double* b, double* w, double* v, double* r);
	
	// after one outer iteration, the messages should be cleaned.
	void					RefreshMemory();

	// test function
	void					PrintTreeStructure();

	//public member variable
	bool					is_chain_;
	
	// chain solver related. If this tree is a chain,
	// f_idx, p_idx_ stores the indexes of the nodes and edges in this tree.
	size_t*				f_idx_;
	size_t*				p_idx_;
	
private:
	// search the leaf nodes in current tree.
	void										SearchLeavesInTree();

	// If current tree is a chain, find the new root (should be an end point).
	// Also prepare some arrays for chain solver.
	void										ReconstructTree();

	// for each node in the tree, set  the message node or create new queue to store the message.
	void										AllocateQueueMemory();

	// pass message for backward solver.
	double 									ForwardPassing(double* a, double* b, double* w);
	double 									BackwardPassing(double* a, double* b, double *w, double *v);

	// tree structure related
	TreeNode*								root_;
	size_t									size_;
	std::vector<TreeNode*>	leaves_;

	// arrays used for some tree operations
	size_t									number_nodes_in_graph_;
	TreeNode**							bfs_;
	TreeNode**							visit_leaves_queue_;
	size_t 									queue_length_;
	size_t									bfs_size_;

	// the arrays used for chain backward solver.
	double*									chain_tm_;
	double*									chain_tp_;
	double*									chain_message_;
};

#endif
