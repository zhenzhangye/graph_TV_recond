#ifndef __FOREST_PRECONDITION_H__
#define __FOREST_PRECONDITION_H__

#include <iostream>
#include <map>
#include <queue>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>

#include "tree.h"
#include "../typedefs.h"


// the algorithms used for partition a graph into forests
enum NestedForestGreedyAlgorithm {
  NF_KRUSKAL,
  NF_BREADTH_FIRST,
  NF_DEPTH_FIRST, 
};

// the forest precondition class
class ForestPrecondition{
public:
	// default constructor.
	ForestPrecondition();
	// default destructor.
	virtual ~ForestPrecondition();	

	// if given graph is a grid, create chains along each dimension.
	bool CreateChains(Graph& graph);

	// set inactive weight for each edge in graph and generate the MST.
	void CreateInitialForest(Graph& graph, bool use_chain, double* pk);

	// partition the graph into forests, each tree is forced to be a chain.
	void CreateInitialForestForceChain(Graph& graph);

	// parallel call the backward solver on each tree. to get p^{k+1}
	void BackwardStepPDHG(double* pk, double* u_bar, double* w, double t, double* result_p, double* a);

	// allocate memory space for arrays used in backward solver.
	void AllocateMemory(Graph& graph);

	// if new partition is generated, release all the memory for previous trees.
	void ClearMemory();

	// (test function) generate the two subgraphs, one with inactive edges, another with active.
	void GenerateSubgraphs(double* pk);

	// Given a subgraph, run DFS to partition it into forests (VERY SLOW)
	void InitialForestForceChain(Graph* forest);

	// print the tree structure, for debuging.
	void PrintTreeInformation();

	// return the number of forests for this graph.
	size_t GetNumberOfForests();

	size_t 	number_edges() 			const { return number_edges_; }
	size_t 	number_vertices() 	const { return number_vertices_; }
	int 	 	number_forests()		const { return number_trees_; }

protected:
	// store the corresponding graph
	std::shared_ptr<Graph> 											graph_; 				
	size_t																			number_edges_;		
	size_t																			number_vertices_;
	size_t																			number_trees_;

	//store the subgraph of each forest partitioned from this graph
	std::vector<std::shared_ptr<Graph> >				forests_;

private:
	// compute the weight of each edge for patition. w=1-|1-|p^k||
	void									SetInactiveWeight(double* pk);

	// correspond the edge to the forest. (i.e. give edge, return the index of that forest).
	std::map<Edge, int> 	edge_to_forest_map_;

	// sotre the pointers of each tree.
	std::vector<Tree*>		trees_;

	// partition the graph by using Prim. (VERY SLOW)
	size_t 								PrimMST(size_t start_vertice, Graph& g);
	size_t								number_forests_;

	// some variables used for improving the partition speed.
	bool*									vertices_visited_;
	bool*									edges_visited_;
	size_t*								org_to_tree_ind_;
	size_t*								prim_source_;
	size_t*								prim_target_;
	size_t*								degree_in_active_graph_;
	size_t*								degree_in_inactive_graph_;
	Graph*								active_graph_;
	Graph*								inactive_graph_;
	static const size_t 	INF_ = 0x3f3f3f3f3f;

	// arrays used for backward solver. 2-dim due to parallel.
	double**							unary_b_;
	double**							chain_weight_;
	double**							dual_result_;
};
#endif
