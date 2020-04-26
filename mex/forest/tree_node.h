/* node in the tree, head file
 *
 * the method and varaibles of a node in the tree structure.
 *
 */

#ifndef __TREENODE_H__
#define __TREENODE_H__

#include <vector>

#include "../pairing_heap/double_ended_priority_queue.h"
#include "../pairing_heap/break_point.h"

class TreeNode{
public:
	// constructor and deconstructor
	TreeNode();
	// constructor: set parent, node index, the index of edge to parent and unary index.
	TreeNode(TreeNode* parent, size_t node_index, size_t edge_index, size_t f_index);
	
	~TreeNode();

	
	//public functions
	//tree structure related functions
	void				AddChildNode(TreeNode* child_node);
	void				AddParentNode(TreeNode* parent_node);
	void				setEdgeIndex(size_t edge_index);
	size_t			getEdgeIndex();
	size_t			getNodeIndex();
	size_t			getNumberChildren();

	// decrease the number of children by 1.
	void				DecreaseNumberChildren();

	// initialize the double ended priority queue for storing message.
	void				InitializeDEPQ();
	
	// set the node where current node can store/get message.
	void				setMessageNode(TreeNode* node);

	TreeNode*		getMessageNode();
	TreeNode*		getFirstChild();
	
	// refresh the memory of current tree node.
	void				RefreshMemory();
	
	// clear the heap of message.
	void				ClearHeap();

	//backward passing related functions
	void				SumMessageWithChildren();
	void				ClipMessage(double a, double b, double* w);
	void				ComputeOptimalRootValue(double a, double b, double* w);	

	// test functions. used for printing information.
	void				PrintTreeNodeInformation();
	void				PrintMessage();

	//public member variables
	TreeNode* 															parent_;
	std::vector<TreeNode*>									children_;
	DoubleEndedPriorityQueue<BreakPoint*>*	message_;

	// variables used in backward solver.
	// lambdas are the varaibles for clipping
	// computed children used to record when the parent should be inserted to bfs queue.
	double																	lambda_pos_;
	double																	lambda_neg_;
	size_t																	computed_number_of_children_;

private:
	// used during forward passing. clip the message from left side and right side.
	void				ClipLeftSide(double a, double b, double* w, double target, double sum_children_weight);
	void				ClipRightSide(double a, double b, double* w, double target, double sum_children_weight);

	// node index in original graph
	size_t			node_index_;
	// the index of edge to parent in original graph
	size_t 			edge_index_;
	// the index in f array
	size_t			f_index_;
	size_t			number_of_children_;

	// the node storing message
	TreeNode*		message_node_;
};

#endif
