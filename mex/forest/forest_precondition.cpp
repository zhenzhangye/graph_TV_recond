/* forest precondition cpp file
 *
 * define the functions in forest precondition class.
 *
 */

#include <boost/graph/connected_components.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include <mex.h>
#include <time.h>
#include <queue>
#include <omp.h>
#include <cmath>
#include <stack>

#include "forest_precondition.h"
#include "../utils.h"

// default constructor
// set all as 0 or null
ForestPrecondition::ForestPrecondition(){
	graph_ = NULL;
	number_trees_ = 0;
	number_forests_ = 0;
	vertices_visited_ = NULL;
}

// default deconstructor
// delete all allocated arrays.
ForestPrecondition::~ForestPrecondition(){
	for(size_t i = 0; i<trees_.size(); ++i){
		delete trees_[i];
	}

	for(size_t i =0; i<MAX_CORE; ++i){
		delete [] unary_b_[i];
		delete [] chain_weight_[i];
		delete [] dual_result_[i];
	}
	delete [] unary_b_;
	delete [] chain_weight_;
	delete [] dual_result_;
	delete [] vertices_visited_;
	delete [] org_to_tree_ind_;
	delete [] prim_source_;
	delete [] prim_target_;
	forests_.clear();
}

// if the graph is a grid, partition it into chains along the dimension.
// input: the graph
// output: partitioned chains, stored in forests.
// 				 Return false if this graph is not a grid.
bool ForestPrecondition::CreateChains(Graph& graph){
	// allocate some memory first
	AllocateMemory(graph);

	WeightMap weights = boost::get(boost::edge_weight_t(), *graph_);

	boost::graph_traits<Graph>::adjacency_iterator ai, a_end;
	boost::graph_traits<Graph>::adjacency_iterator ai2, a_end2;
	boost::graph_traits<Graph>::out_edge_iterator ei, e_end;

	// find vertex with minimum degree
	size_t min_degree = boost::degree(0, *graph_);
	size_t vertex_idx = 0;

	for(int i = 1; i < boost::num_vertices(*graph_); ++i) 
		if(min_degree > boost::degree(i, *graph_)) {
			vertex_idx = i;
			min_degree = boost::degree(i, *graph_);
		}

	//mexPrintf("min_degree: %d\n", min_degree);

	int max_distance = 0;
	int *distances = new int[boost::num_vertices(*graph_)];
	std::map<Edge, int> edge_labels;
	std::map<int, std::vector<int> > bfs_levels;
	std::fill_n(distances, boost::num_vertices(*graph_), -1);

	std::queue<size_t> bfs;

	// run BFS starting from that vertex, record distance to start vertex
	bfs.push(vertex_idx);
	distances[vertex_idx] = 0;
	while(!bfs.empty()){
		size_t current = bfs.front();
		bfs.pop();

		bfs_levels[distances[current]].push_back(current);
		max_distance = std::max(max_distance, distances[current]);

		boost::tie(ai,a_end) = boost::adjacent_vertices(current, *graph_);
		for(;ai != a_end; ++ai){
			if(distances[*ai] != -1)
				continue;

			distances[*ai] = distances[current] + 1;
			bfs.push(*ai);
		}
	}

	// give incident edges of root node different labels
	int color = 0;
	for (boost::tie(ei, e_end) = boost::out_edges(vertex_idx, *graph_); ei != e_end; ++ei) {
		edge_labels[*ei] = color++;
	} 
	// something went very wrong if this is not true ...
	if(color != min_degree)
		return false;

	for(int i = 2; i <= max_distance; ++i) {

		//mexPrintf("distance %d: vertices=%d\n", i, bfs_levels[i].size());
		for(int j = 0; j < bfs_levels[i].size(); ++j) {

			int u = bfs_levels[i][j];
			std::vector<int> down_vertices;

			// let D = { (u,v) | dist(u) == i, dist(v) == i - 1 } be the down edges
			for (boost::tie(ei, e_end) = boost::out_edges(u, *graph_); ei != e_end; ++ei) {
				int v1 = boost::source(*ei, *graph_);
				int v2 = boost::target(*ei, *graph_);

				if(v1 != u)
					std::swap(v1, v2);

				if(distances[v2] == distances[u] - 1)
					down_vertices.push_back(v2);
			}

			//mexPrintf("  vertex %d: down_vertices=%d\n", j, down_vertices.size());

			// if |D| == 1 then
			if(down_vertices.size() == 1) {
				int v = down_vertices[0];

				for (boost::tie(ei, e_end) = boost::out_edges(v, *graph_); ei != e_end; ++ei) {
					int v1 = boost::source(*ei, *graph_);
					int v2 = boost::target(*ei, *graph_);

					if(v1 != v)
						std::swap(v1, v2);

					if(distances[v2] == distances[v] - 1) {
						if(edge_labels.count(*ei) == 0) // down edge of v is unlabeled -- shouldn't happen
							return false;

						// label edge (u, v) in the same color as (v, v2)
						edge_labels[boost::edge(u, v, *graph_).first] = edge_labels[*ei];
						//mexPrintf("  -> labeled edge (%d,%d)\n", u, v);
					}
				}
			}
			else {
				// for each unlabeled (u,v) in D
				for(int k = 0; k < down_vertices.size(); ++k) {
					int v = down_vertices[k];

					// continue if the edge (u, v) is already labeled
					if(edge_labels.count(boost::edge(u, v, *graph_).first) != 0)
						continue;

					// try to find edges (u, v) and (u, w) which have common neighbour x 
					int x = -1, w;
					for(int l = 0; l < down_vertices.size(); ++l) {
						if(l == k)
							continue; 

						w = down_vertices[l];

						// check if v,w2 have a common neighbour x != u
						boost::graph_traits<Graph>::adjacency_iterator ai_start, ai_start2;
						boost::tie(ai_start, a_end) = boost::adjacent_vertices(v, *graph_);
						boost::tie(ai_start2, a_end2) = boost::adjacent_vertices(w, *graph_);

						for(ai = ai_start; ai != a_end; ++ai) 
							for(ai2 = ai_start2; ai2 != a_end2; ++ai2) 
								if(*ai == *ai2 && *ai != u) {
									x = *ai;
									goto found_x;
								}
					}

					//mexPrintf("    -> Couldn't find common neighbour x.\n");
					// couldn't find a common neighbour
					return false;

found_x:
					//mexPrintf("    u=%d, v=%d, w=%d, x=%d\n", u, v, w, x);

					if(edge_labels.count(boost::edge(v, x, *graph_).first) == 0)
						return false;

					if(edge_labels.count(boost::edge(w, x, *graph_).first) == 0)
						return false;

					int wx = edge_labels[boost::edge(w, x, *graph_).first];
					int vx = edge_labels[boost::edge(v, x, *graph_).first];

					// check if there is inconsistent labeling
					if(edge_labels.count(boost::edge(u, v, *graph_).first) != 0)
						if(wx != edge_labels[boost::edge(u, v, *graph_).first])
							return false; 

					if(edge_labels.count(boost::edge(u, w, *graph_).first) != 0)
						if(vx != edge_labels[boost::edge(u, w, *graph_).first])
							return false; 

					// label the edge (u, v) into the same as (w, x) and (u,w) same as (v,x)
					edge_labels[boost::edge(u, v, *graph_).first] = wx;
					edge_labels[boost::edge(u, w, *graph_).first] = vx;
					//mexPrintf("  -> labeled edges (%d,%d) and (%d,%d)\n", u, v, u, w);
				}
			}
		}
	}

	// create min_degree forests
	for(int i = 0; i < min_degree; ++i) {
		Graph *forest = new Graph(boost::num_vertices(*graph_));
		forests_.push_back( std::shared_ptr<Graph>(forest) );
	}

	// put edge with color i into forest i
	auto es = boost::edges(*graph_);
	for (auto eit = es.first; eit != es.second; ++eit) {
		const int v1 = boost::source(*eit, *graph_);
		const int v2 = boost::target(*eit, *graph_);

		const Edge& edge_in_graph = boost::edge(v1, v2, *graph_).first;

		if(edge_labels.count(edge_in_graph) == 0)
			return false;

		const int forest_idx = edge_labels[edge_in_graph];

		if(forest_idx < 0 || forest_idx >= min_degree)
			return false;

		if(v1 < v2)
			boost::add_edge(v2, v1, EdgeWeightProperty( boost::get(weights, edge_in_graph) ), *(forests_[forest_idx]));
		else
			boost::add_edge(v1, v2, EdgeWeightProperty( boost::get(weights, edge_in_graph) ), *(forests_[forest_idx]));

		edge_to_forest_map_[edge_in_graph] = forest_idx;
	}

	// verify if result is a grid
	int N = 1;
	int *counts = new int[boost::num_vertices(*graph_)];
	for(int i = 0; i < min_degree; ++i) {
		memset(counts, 0, sizeof(int) * boost::num_vertices(*graph_));

		std::vector<int> component(boost::num_vertices(*forests_[i]));
		int num = boost::connected_components(*forests_[i], &component[0]);

		for(int i = 0; i < component.size(); ++i)
			++counts[component[i]];

		for(int i = 0; i < num - 1; ++i)
			if(counts[i] != counts[i+1])
				mexErrMsgTxt("create_chains: Failure! Chains in same forest have different length.");

		N *= num;
	}

	if(N != boost::num_vertices(*graph_))
		mexErrMsgTxt("create_chains: Failure! Most likely input was not a grid.");

	delete [] distances;
	delete [] counts;

	number_forests_ = forests_.size();

	// now we create the trees for each forest
	size_t visited_edges, total_edges;
	Graph* g;
	size_t tree_number_vertices;
	size_t current;
	bool use_chain = true;
	while(!bfs.empty()){
		bfs.pop();
	}

	for(size_t i =0; i<forests_.size(); ++i){
		g = forests_[i].get();
		weights = boost::get(boost::edge_weight_t(), *g);
		visited_edges = boost::num_edges(*g);
		total_edges = visited_edges;
		std::fill_n(vertices_visited_, number_vertices_, false);

		//For current graph, stop if all vertices are visited or all edges are visited.
		for (size_t j =0; j<number_vertices_ && visited_edges>0; j++){
			tree_number_vertices = 0;
			// if current vertex is not visited
			if (!vertices_visited_[j]){
				boost::tie(ai, a_end) = boost::adjacent_vertices(j, *g);

				// if this vertex is not isolated in this graph.
				if (ai!=a_end){
					// perform a bfs starting from current node.
					bfs.push(*ai);
					while(!bfs.empty()){
						current = bfs.front();
						bfs.pop();
						vertices_visited_[current] = true;
						
						boost::tie(ai, a_end) = boost::adjacent_vertices(current, *g);
						for(; ai!=a_end; ++ai){
							if (vertices_visited_[*ai]){
								continue;
							}
							bfs.push(*ai);
							if (*ai<current){
								prim_source_[tree_number_vertices] = *ai;
								prim_target_[tree_number_vertices] = current;
							}else{
								prim_source_[tree_number_vertices] = current;
								prim_target_[tree_number_vertices] = *ai;
							}
							tree_number_vertices++;
						}
					}
					trees_.push_back(new Tree(graph_.get(), prim_source_, prim_target_, weights, trees_.size(), j, tree_number_vertices, use_chain));
				}
			}
		}
	}
	number_trees_ = trees_.size();
	return true;
}

// partition the graph into chains.
// input: graph (pointer), the boost library graph pointer.
// output: partitioned graph, each tree is a chain, stored in trees_.
void ForestPrecondition::CreateInitialForestForceChain(Graph& graph){
	// allocate memory for the arrays
	AllocateMemory(graph);

	WeightMap weights = boost::get(boost::edge_weight_t(), *graph_);
	bool *visited = new bool[number_vertices_];

	Graph g(*graph_);

	while(boost::num_edges(g)>0){
		Graph *forest = new Graph(number_vertices_);
		std::vector<Edge> spanning_tree;
		boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

		for (int i = 0; i<spanning_tree.size(); ++i){
			const Edge& e = spanning_tree[i];

			int v1 = boost::source(e, g);
			int v2 = boost::target(e, g);

			const Edge& edge_in_graph = boost::edge(v1, v2, *graph_).first;
			if (v1 < v2)
				boost::add_edge(v2, v1, EdgeWeightProperty( boost::get(weights, edge_in_graph) ), *forest);
			else
				boost::add_edge(v1, v2, EdgeWeightProperty( boost::get(weights, edge_in_graph) ), *forest);

			edge_to_forest_map_[edge_in_graph] = forests_.size();
		}

		// turn spanning tree into linear forest by removing edges
		std::map<int, int>  degree_map;
		std::map<Edge, bool> delete_map;

		auto vs = boost::vertices(*forest);
		for (auto vit = vs.first; vit != vs.second; ++vit)
			degree_map[*vit] = boost::degree(*vit, *forest);

		// flag edges for deletion
		auto eds = boost::edges(*forest);
		for (auto eit = eds.first; eit != eds.second; ++ eit){
			const int v1 = boost::source(*eit, *forest);
			const int v2 = boost::target(*eit, *forest);

			if (degree_map[v1] > 2 || degree_map[v2] > 2){
				delete_map[*eit] = true;
			}else{
				delete_map[*eit] = false;
			}
		}

		// remove edges
		for (auto it = delete_map.begin(); it != delete_map.end(); ++it)
			boost::remove_edge(it->first, *forest);

		// check if edges from g can be added to the linear forest without destroying it
		auto es = boost::edges(g);
		for (auto eit = es.first; eit != es.second; ++eit){
			const int v1 = boost::source(*eit, g);
			const int v2 = boost::target(*eit, g);

			// check if edge is in forest already, if true, continue.
			if (boost::edge(v1, v2, *forest).second)
				continue;

			if (v1 == v2)
				continue;

			const int d1 = boost::degree(v1, *forest);
			const int d2 = boost::degree(v2, *forest);

			if (d1 + d2 <= 1){
				const Edge& edge_in_graph = boost::edge(v1, v2, *graph_).first;
				if (v1 < v2)
					boost::add_edge(v2, v1, EdgeWeightProperty(boost::get(weights, edge_in_graph)), *forest);
				else
					boost::add_edge(v1, v2, EdgeWeightProperty(boost::get(weights, edge_in_graph)), *forest);
				edge_to_forest_map_[edge_in_graph] = forests_.size();
			}else if (d1 == 1 & d2 == 1){
				// Check if there is a path from v1 to v2 in forest
				// if not -> add edge
				memset(visited, 0, number_vertices_ * sizeof(bool));
				boost::graph_traits<Graph>::adjacency_iterator ai, a_end;
				std::queue<int> bfs;
				bfs.push(v1);
				int cur_v;

				while(!bfs.empty()){
					cur_v = bfs.front();
					bfs.pop();

					boost::tie(ai,a_end) = boost::adjacent_vertices(cur_v, *forest);
					for(; ai != a_end; ai++){
						if (visited[*ai])
							continue;

						visited[*ai] = true;
						bfs.push(*ai);
					}
				}

				if (cur_v != v2){
					const Edge& edge_in_graph = boost::edge(v1, v2, *graph_).first;
					if (v1<v2)
						boost::add_edge(v2, v1, EdgeWeightProperty(boost::get(weights, edge_in_graph)), *forest);
					else
						boost::add_edge(v1, v2, EdgeWeightProperty(boost::get(weights, edge_in_graph)), *forest);
					edge_to_forest_map_[edge_in_graph] = forests_.size();
				}
			}
		}

		// subtract forest from graph g
		es = boost::edges(*forest);
		for (auto eit = es.first; eit != es.second; ++eit){
			const int v1 = boost::source(*eit, *forest);
			const int v2 = boost::target(*eit, *forest);
			boost::remove_edge(v1, v2, g);
		}

		forests_.push_back( std::shared_ptr<Graph>(forest) );
		number_forests_++;
		InitialForestForceChain(forest);

	}

	delete [] visited;
}
	
// partition the graph and initialize the trees.
// input: graph (pointer), the boost library graph pointer
// 			  use_chain (bool), if it is true, when a tree is a chain, call more efficient backward solver.
// 			  pk (|E|*1 vector), here is all zero. The partition will be bfs.
// output: partitioned graph, each tree is stored in trees_.
void ForestPrecondition::CreateInitialForest(Graph& graph, bool use_chain, double* pk){
	// allocate memmory for the arrays
	AllocateMemory(graph);

	// set the weight of each edge for partitioning.
	SetInactiveWeight(pk);

	WeightMap weights = boost::get(boost::edge_weight_t(), *graph_);

	Graph g(*graph_);

	// iterate all edges in the graph.
	while(boost::num_edges(g)>0){
		Graph *forest = new Graph(number_vertices_);
		std::vector<Edge> spanning_tree;
		boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
		for (int i =0; i<spanning_tree.size(); ++i){
			const Edge& e = spanning_tree[i];

			int v1 = boost::source(e,g);
			int v2 = boost::target(e,g);

			const Edge& edge_in_graph = boost::edge(v1, v2, *graph_).first;
			if (v1 < v2)
				boost::add_edge(v2, v1, EdgeWeightProperty( boost::get(weights, edge_in_graph) ), *forest);
			else
				boost::add_edge(v1, v2, EdgeWeightProperty( boost::get(weights, edge_in_graph) ), *forest);

			edge_to_forest_map_[edge_in_graph] = forests_.size();
		}

		auto es = boost::edges(*forest);
		for (auto eit = es.first; eit!= es.second; ++eit){
			const int v1 = boost::source(*eit, *forest);
			const int v2 = boost::target(*eit, *forest);
			boost::remove_edge(v1, v2, g);
		}
		forests_.push_back( std::shared_ptr<Graph> (forest));
		number_forests_ ++;

		InitialForestForceChain(forest);

		
		/*
		std::fill_n(vertices_visited_, number_vertices_, false);

		for (size_t i =0; i< number_vertices_; ++i){
			if(vertices_visited_[i])
				continue;

			if(boost::out_degree(i,g) == 0)
				continue;

			// perform prim to get MST on current graph.
			size_t tree_number_vertices = PrimMST(i,g);
		
			// initialize the tree and add it into the trees vector.
			trees_.push_back(new Tree(graph_.get(), &g, prim_source_, prim_target_, weights, trees_.size(), i, tree_number_vertices, use_chain)); 	
			++number_trees_;
		}
		*/
	}
	//number_forests_ = forests_.size();
}

// call the backward solver for the dual update.
// input: pk (|E|*1, array): the dual variable p at iteration k.
// 			  u_bar (|V|*1, array): the over-relaxation primal variabl u.
// 			  w (|E|*1, array): the weight of each edge.
// 			  t (scalar): the stepsize for dual updating.
// 			  result_p (|E|*1, array): the result after duall update. p^{k+1}.
// 			  a (|V|*1, array): the weight infront of data term.
void ForestPrecondition::BackwardStepPDHG(double* pk, double* u_bar, double* w, double t, double* result_p, double* a){
// parallel compute on each tree.
#pragma omp parallel for num_threads(MAX_CORE)
	for(size_t i = 0; i<trees_.size(); ++i){
		double* this_unary_b;
	 	double*	this_chain_weight; 
		double* this_dual_result;
		size_t core_ind;

		// avoid the memory conflict in parallel.
		core_ind = omp_get_thread_num()%MAX_CORE;
		this_unary_b = unary_b_[core_ind];
		this_chain_weight = chain_weight_[core_ind];
		this_dual_result = dual_result_[core_ind];

		// if current tree is a chain
		if (trees_[i]->is_chain_){
			// compute the dual form of the dual subproblem. f = K^T * p_l^k + u_bar/t
			trees_[i] -> ComputeFChain(u_bar, pk, w, t, this_unary_b, this_chain_weight);

			// apply efficient backward solver on chain. u = argmin a/2 *(u+f)^2 + ||Ku||_1
			trees_[i] -> TV1DWeightedKPR(a, this_unary_b, this_chain_weight, this_dual_result, trees_[i]->size());

			// recover the dual variabl p^{k+1} from its dual solution. nabla_w^T p = u+f
			trees_[i] -> SolveLinearSystemChain(this_unary_b, w, this_dual_result, result_p);
		}else{
			// compute the dual form of the dual subproblem.  f = K^T * p_l^k + u_bar/t
			trees_[i] -> ComputeF(u_bar, pk, w, t, this_unary_b);

			//std::cout<<"#############################################"<<std::endl;
			// apply backward solver on a tree. u = argmin a/2 * (u+f)^2 + ||Ku||_1
		  trees_[i] -> BackwardSolver(a, this_unary_b, w, this_dual_result);

			// recover the dual variabl p^{k+1} from its dual solution. nabla_w^T p = u+f
			trees_[i] -> SolveLinearSystem(this_unary_b, w, this_dual_result, result_p);

			// clear some memory for next iteration.
			trees_[i] -> RefreshMemory();
		}
	}
}

// allocate memory for current preconditioner.
// input: graph.
// output: the arrays might be used.
void ForestPrecondition::AllocateMemory(Graph& graph){
	// if the graph is not created
	if(graph_ == NULL){
		// graph memory
		graph_ = std::shared_ptr<Graph>(new Graph(graph));

		// get the number of vertices, edges in the graph.
		number_vertices_ = boost::num_vertices(*graph_);
		number_edges_    = boost::num_edges(*graph_);
		number_trees_  = 0;

/*		// two subgraphs. One stores the active edges and another with inactive.
		inactive_graph_ = new Graph(graph);
		active_graph_ = new Graph(number_vertices_);
		degree_in_active_graph_ = new size_t[number_vertices_];
		degree_in_inactive_graph_ = new size_t[number_vertices_];
		for(size_t i = 0; i<number_vertices_; ++i){
			degree_in_inactive_graph_[i] = boost::out_degree(i, *inactive_graph_);
			degree_in_active_graph_[i] = 0;
		}
		*/

		// the arrays used for efficient backward solver on a chain.
		unary_b_ 			= new double*[MAX_CORE];
		chain_weight_ = new double*[MAX_CORE];
		dual_result_  = new double*[MAX_CORE];
		for(int i = 0; i<MAX_CORE; ++i){
			unary_b_[i] 			= new double[number_vertices_];
			chain_weight_[i] 	= new double[number_vertices_];
			dual_result_[i]		= new double[number_vertices_];
		}

		// these arrays used for prim to generate MST.
		vertices_visited_ = new bool[number_vertices_];
		edges_visited_		= new bool[number_edges_];
		org_to_tree_ind_	= new size_t[number_vertices_];
		prim_source_			= new size_t[number_vertices_];
		prim_target_			= new size_t[number_vertices_];
	}
}

// given a p^k, set the partition weight of each edge.
// input: p^k (|E|*1, array): the dual variable at iteration k.
// output: compute the partition weight, w_partition = 1- ||pk|-1|
void ForestPrecondition::SetInactiveWeight(double* pk){
	WeightMap weights = boost::get(boost::edge_weight_t(), *graph_);

	boost::graph_traits<Graph>::edge_iterator ei, e_end;

	// iterate over all edges and compute the inactive_weight.
	for(boost::tie(ei, e_end) = boost::edges(*graph_); ei != e_end; ++ei){
		size_t index = boost::get(weights, *ei).original_index_;
		double inactive_weight = 1 - fabs(fabs(pk[index]) - 1);
		boost::get(weights, *ei).inactive_weight_ = inactive_weight;
	}
}

// Given a start vertice and the graph. Generate a minimum spanning tree by prim.
// The start vertices is included.
// input: start_vertice(int): the index of starting vertice
// 				graph
// output: number of vertices in this tree.
// 				 results are stored in prim_source_ prim_target_.
size_t ForestPrecondition::PrimMST(size_t start_vertice, Graph& g){
	// get the weights of each edge
	WeightMap weights = boost::get(boost::edge_weight_t(), g);

	// the priority queue used in prim
	std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair> > pq;

	std::vector<double> key(number_vertices_, INF_);

	// used to store the indexes of vertices in result tree.
	std::fill_n(org_to_tree_ind_, number_vertices_, number_vertices_+5);
	
	// push the start vertice to the pq.
	pq.push(std::make_pair(0, start_vertice));
	key[start_vertice] = 0;
	size_t vertice_ind = 0;
	size_t src_tag_ind = 0;

	while(!pq.empty()){
		size_t u = pq.top().second;
		pq.pop();

		vertices_visited_[u] = true;

		boost::graph_traits<Graph>::adjacency_iterator ai, a_end;
		boost::tie(ai, a_end) = boost::adjacent_vertices(u, g);

		for(;ai != a_end; ++ai){
			vertice_ind = *ai;
			double inactive_weight = boost::get(weights, boost::edge(vertice_ind, u, g).first).inactive_weight_;

			// if current vertice is not visited and the cost is smaller.
			if(vertices_visited_[vertice_ind] == false && key[vertice_ind] > inactive_weight){
				// update the weight and priority queue
				key[vertice_ind] = inactive_weight;
				pq.push(std::make_pair((double)key[vertice_ind], vertice_ind));

				// save the indexes  into source and target. we need to ensure one edge always starts from vertice with smaller index.
				if(org_to_tree_ind_[vertice_ind] > number_vertices_){
					org_to_tree_ind_[vertice_ind] = src_tag_ind;
					prim_source_[src_tag_ind] = u;
					prim_target_[src_tag_ind] = vertice_ind;
					++src_tag_ind;
				}else{
					prim_source_[org_to_tree_ind_[vertice_ind]] = u;
				}
			}
		}
	}	

	return src_tag_ind;
}

// Test function: print the structure of trees.
void ForestPrecondition::PrintTreeInformation(){
	std::cout<<"There are "<<number_forests_<<" forests"<<std::endl;
	std::cout<<"There are "<<number_trees_<<"("<<trees_.size()<<") trees"<<std::endl;
	for(size_t i =0; i<trees_.size(); ++i){
		std::cout<<"Tree#"<<i<<": "<<std::endl;
		trees_[i]->PrintTreeStructure();
	}
}

// Clear the memory of trees.
void ForestPrecondition::ClearMemory(){
	for(size_t i = 0; i<trees_.size(); ++i){
		// call tree deconstructor
		delete trees_[i];
	}
	trees_.clear();
	forests_.clear();
	number_trees_ = 0;
	number_forests_ = 0;
}

// generate the two subgraphs from original one.
// one graph with active edges and another with inactive.
// input: p^k (|E|*1, array): the dual variable at iteration k.
// output: two subgraphs stored in inactive_graph_ and active_graph.
void ForestPrecondition::GenerateSubgraphs(double* pk){
	WeightMap weights = boost::get(boost::edge_weight_t(), *graph_);

	boost::graph_traits<Graph>::edge_iterator ei, e_end;

	for(boost::tie(ei, e_end) = boost::edges(*graph_); ei != e_end; ++ei){
		size_t index = boost::get(weights, *ei).original_index_;
		boost::get(weights, *ei).inactive_weight_ = 1 - fabs(fabs(pk[index])-1) ;
		size_t source = boost::source(*ei, *graph_);
		size_t target = boost::target(*ei, *graph_);
		//if this edge is an inactive one. Check if it is in the inactive graph.
		//If yes do nothing. Otherwise, delete it from active graph and insert into inactive one.
		if(fabs(1-pk[index])>THRESHOLD){
			if(!boost::edge(source, target, *inactive_graph_).second){
				boost::add_edge(source, target, EdgeWeightProperty(boost::get(weights, *ei)), *inactive_graph_);
				boost::remove_edge(source, target, *active_graph_);
			}	
		}else{
			if(!boost::edge(source, target, *active_graph_).second){
				boost::add_edge(source, target, EdgeWeightProperty(boost::get(weights, *ei)), *active_graph_);
				boost::remove_edge(source, target, *inactive_graph_);
			}
		}

		degree_in_inactive_graph_[source] = boost::out_degree(source, *inactive_graph_);
		degree_in_inactive_graph_[target] = boost::out_degree(target, *inactive_graph_);
		degree_in_active_graph_[source] = boost::out_degree(source, *active_graph_);
		degree_in_active_graph_[target] = boost::out_degree(target, *active_graph_);

	}

}

// used for partitioning subgraphs into forests. allocate memory as well.
// input: active_inactive(bool): true if we want to initialize active graph. false if inactive
// 				use_chain(bool): true if we call efficient backward solver for chains.
// output: partitioned subgraphs.
void ForestPrecondition::InitialForestForceChain(Graph* forest){
	WeightMap weights = boost::get(boost::edge_weight_t(), *graph_);

	std::stack<size_t> dfs;

	size_t number_edges_partitioned = 0;
	size_t number_edges_graph = boost::num_edges(*forest);

	boost::graph_traits<Graph>::adjacency_iterator ai, a_end;

	std::vector<size_t> degree_array(number_vertices_, 0);

	for (size_t i = 0; i<number_vertices_; ++i){
		degree_array[i] = boost::out_degree(i, *forest);
	}

	std::fill_n(edges_visited_, number_edges_, false);
	
	// iterate all the edges in the forest
	while(number_edges_partitioned < number_edges_graph){
		std::fill_n(vertices_visited_, number_vertices_, false);
		// check the degree of each vertice in remaining graph
		for(size_t i = 0; i < number_vertices_; ++i){
			if(boost::out_degree(i, *forest)==0)
				continue;

			// find one edge which is not visited.
			for(boost::tie(ai, a_end) = boost::adjacent_vertices(i, *forest); ai!= a_end; ++ai){
				size_t edge_index = boost::get(weights, boost::edge(i, *ai, *forest).first).original_index_;
				// if this edge is not visited, add current vertice
				if (!edges_visited_[edge_index]){
					dfs.push(i);
					dfs.push(*ai);
					vertices_visited_[i] = true;
					vertices_visited_[*ai] = true;
					edges_visited_[edge_index] = true;
					break;
				}
			}

			// perform dfs staring from that edge
			size_t tree_number_vertices = 0;
			while (!dfs.empty()){
				size_t current_node = dfs.top();
				
				boost::tie(ai, a_end) = boost::adjacent_vertices(current_node, *forest);
				bool found_one_edge = false;
				if(degree_array[current_node]>0){
					for(; ai!= a_end; ++ai){
						size_t edge_index = boost::get(weights, boost::edge(current_node, *ai, *graph_).first).original_index_;
						if (!edges_visited_[edge_index] && !vertices_visited_[*ai]){
							dfs.push(*ai);
							edges_visited_[edge_index] = true;
							vertices_visited_[*ai] = true;
							found_one_edge = true;
							break;
						}
					}
				}

				if(dfs.size() == 1){
					dfs.pop();
					break;
				}

				// if we cannot find at least one edge, we get 
				if(!found_one_edge){
					prim_target_[tree_number_vertices] = dfs.top();
					--degree_array[dfs.top()];
					dfs.pop();
					prim_source_[tree_number_vertices]= dfs.top();
					++tree_number_vertices;
					--degree_array[dfs.top()];
					number_edges_partitioned ++;
				}
			}

			// if this tree contains at least one nodes, we create a new tree for it.
			if(tree_number_vertices>0){
				trees_.push_back(new Tree(graph_.get(), prim_source_, prim_target_, weights, number_forests_, prim_source_[tree_number_vertices-1], tree_number_vertices, true)); 	
				number_trees_++;
			}

			if(tree_number_vertices+1 == number_vertices_)
				break;
		}
	}
}

// return the number of forests
size_t ForestPrecondition::GetNumberOfForests(){
	return number_forests_;
}
