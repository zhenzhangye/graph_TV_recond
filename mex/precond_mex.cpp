/* Mex function file
 * 
 * The gate way function for matlab.
 * Contains c++ functions are called in matlab files.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <queue>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/connected_components.hpp>

#include <mex.h>

#include "utils.h"
#include "forest/forest_precondition.h"

using std::shared_ptr;
using std::string;
using std::map;
using std::vector;
using std::function;
using std::stringstream;

static vector<shared_ptr<Graph> > graphs_;
static vector<shared_ptr<ForestPrecondition> > forest_preconds_;

#define MEX_ARGS int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs

// Given a sparse adjacency matrix, construct the graph
// input: matrix (|V|*|E|): a sparse adjacency matrix of the graph.
// output: graph pointer: a graph structure defined by boost library.
static void CreateGraphFromAdjacency(MEX_ARGS){
	if (!mxIsSparse(prhs[0]))
		mexErrMsgTxt("CreateGraphFromAdjacency: adjacency matrix must be sparse!");

	// read sparse matrix from matlab.
	// size: |V| * |E|
	const mxArray *pm = prhs[0];
	mwIndex *ind = mxGetIr(pm);
	mwIndex *ptr = mxGetJc(pm);
	double *pr = mxGetPr(pm);
	const mwSize *dims = mxGetDimensions(pm);
	const size_t numVertices = dims[0];

	// construct graph
	size_t count = 0;
	Graph *graph = new Graph(numVertices);
	for(size_t idx = 0; idx < numVertices; idx++){
		for (size_t j = ptr[idx];j<ptr[idx+1];j++){
			size_t v1 = idx;
			size_t v2 = ind[j];
			if(v1>v2){
				boost::add_edge(v1, v2, EdgeWeightProperty(EdgeData(pr[j],count++, 0)), *graph);
				//std::cout<<v2<<"->"<<v1<<std::endl;
			}
		}
	}

	graphs_.push_back(shared_ptr<Graph>(graph));
	plhs[0] = ptr_to_handle(graph);
}

// Given a .bk file, construct the graph.
// input: file name (string): directory to the bk file.
// output: graph pointer: a grph structure defined by boost library.
static void CreateGraphFrombkFile(MEX_ARGS){
  std::string binary_name(mxArrayToString(prhs[0]));
  binary_name += ".graph";

  std::ifstream ifs;
  ifs.open(binary_name);

  if(ifs.good()) {

    // read graph and unaries from binary file
    Graph *graph = new Graph;
    graphs_.push_back(shared_ptr<Graph>(graph));
    plhs[0] = ptr_to_handle(graph);
    
    boost::archive::binary_iarchive ia(ifs);
    ia >> *graph;

    //std::cout << "read graph (" << boost::num_vertices(*graph) << " vertices)" << std::endl;

    plhs[1] = mxCreateDoubleMatrix(boost::num_vertices(*graph), 1, mxREAL);
    double *unaries = (double *)mxGetPr(plhs[1]);
    for(int i = 0; i < boost::num_vertices(*graph); i++)
      unaries[i] = (*graph)[i].unary_;

    double energy_offset;
    ia >> energy_offset;    
    plhs[2] = mxCreateDoubleScalar(energy_offset);
    ifs.close();
  }
  else {
    ifs.open(mxArrayToString(prhs[0]));

    if(!ifs.good())
      mexErrMsgTxt("create_graph_from_bk_file: could not open file!");
  
    int flow_nodes, flow_arcs;
    int node1, node2, capacity, capacity2;
    std::string line;
    char c;
    Graph *graph = 0;
    double *unaries = 0;
    double energy_offset = 0;
    size_t count = 0;
    while( std::getline(ifs, line) ) {
      if(line.empty())
	continue;

      switch(line[0]) {
      case 'c':
      case '\n':
      case '\0':
      default:
	break;

      case 'p':
	sscanf(line.c_str(), "%c %d %d", &c, &flow_nodes, &flow_arcs);
	graph = new Graph(flow_nodes);
	graphs_.push_back(shared_ptr<Graph>(graph));
	plhs[0] = ptr_to_handle(graph);
	break;

      case 'n':
	sscanf(line.c_str(), "%c %d %d %d", &c, &node1, &capacity, &capacity2);
	if(capacity == 0 && capacity2 == 0)
	  break;

	if(graph == 0)
	  mexErrMsgTxt("create_graph_from_bk_file: graph not initialized!");

	if(node1 < 0 || node1 >= flow_nodes)
	  mexErrMsgTxt("create_graph_from_bk_file: bad node ID");

	(*graph)[node1].unary_ += (-capacity + capacity2);
	energy_offset += capacity;
	break;

      case 'a':
	sscanf(line.c_str(), "%c %d %d %d %d", &c, &node1, &node2, &capacity, &capacity2);

	if(graph == 0)
	  mexErrMsgTxt("create_graph_from_bk_file: graph not initialized!");

	if(node1 < 0 || node1 >= flow_nodes || node2 < 0 || node2 >= flow_nodes)
	  mexErrMsgTxt("create_graph_from_bk_file: bad node ID");

	bool edge_in_graph =
	  boost::edge(node1, node2, *graph).second ||
	  boost::edge(node2, node1, *graph).second;

	if(!edge_in_graph) {
	  double weight = ((double)capacity + (double)capacity2) / 2.0;

	  if(capacity < 0 || capacity2 < 0)
	    mexErrMsgTxt("create_graph_from_bk_file: negative weight encountered!");

	  if(weight > 0) {
	    if(node1 > node2)
	      boost::add_edge(node1, node2, EdgeWeightProperty( EdgeData(weight, count++, weight) ), *graph);
	    else 
	      boost::add_edge(node2, node1, EdgeWeightProperty( EdgeData(weight, count++, weight) ), *graph);
          
	    (*graph)[node1].unary_ += ((double)capacity - (double)capacity2) / 2.0;
	    (*graph)[node2].unary_ -= ((double)capacity - (double)capacity2) / 2.0;
	  }
	}
	else {
	  mexErrMsgTxt("create_graph_from_bk_file: duplicate edge!");
	}
      
	break;
      }
    }

    ifs.close();
    
    // remove nodes with degree 0 
    bool done = false;
    while(!done) {
      done = true;
      for(int i = 0; i < boost::num_vertices(*graph); i++)
	if(boost::degree(i, *graph) == 0) {
	  if((*graph)[i].unary_ < 0) {
	    energy_offset += (*graph)[i].unary_;
	  }

	  boost::remove_vertex(i, *graph);
	  done = false;
	  break;
	}
    }
  
    plhs[1] = mxCreateDoubleMatrix(boost::num_vertices(*graph), 1, mxREAL);
    unaries = (double *)mxGetPr(plhs[1]);

    for(int i = 0; i < boost::num_vertices(*graph); i++)
      unaries[i] = (*graph)[i].unary_;
  
    plhs[2] = mxCreateDoubleScalar(energy_offset);

    // store computed graph in binary format
    std::ofstream ofs;
    ofs.open(binary_name);
    if(!ofs.good()) mexErrMsgTxt("create_graph_from_bk_file: failed to write binary.");
    
    boost::archive::binary_oarchive oa(ofs);
    oa << *graph;
    oa << energy_offset;

    ofs.close();
  }  

}

// Partition the graph and intialize the forest data structure.
// input: graph pointer: the boost graph.
//				use_chain (bool): true: a fatser backward solver will be called later.
//				into_chain (bool): true: partition the graph into chains only for grid graphs.
// output: preconditioner pointer: the forest preconditioner handler.
static void PartitionInitialGraph(MEX_ARGS){
	if(nrhs < 4){
		mexErrMsgTxt("PartitionInitialGraph: expected four arguments \n");
	}

	Graph& graph = *handle_to_ptr<Graph>(prhs[0]);

	// when use_chain is true, an efficient back ward solver will be applied for the chains in the forests.
	bool use_chain = (bool)mxGetScalar(prhs[1]);

	// when into_chain is true, grid will be partitioned into chains.
	bool into_chain = (bool)mxGetScalar(prhs[2]);

	// during the initialization, pk is all zero.
	double *pk = (double*)mxGetPr(prhs[3]);

	ForestPrecondition *precond;
	precond = new ForestPrecondition;

	// initialize the forest by (p^k).
	if (into_chain){
		if (precond-> CreateChains(graph)){
			std::cout<< "The graph is a grid, chains are created!"<<std::endl;
		}else{
			mexErrMsgTxt("This is not a grid, failed to create chains");
		}
	}else{
		precond-> CreateInitialForest(graph, use_chain, pk);
		std::cout<< "The graph is partitioned into trees"<<std::endl;
	}
	
	//precond->PrintTreeInformation();

	plhs[0] = ptr_to_handle(precond);

}

// Partition the graph and initialize the forest, each tree is forced to be a chain.
// input: graph pointer: the boost graph.
// output: preconditioner pointer: the forest preconditioner handler.
static void PartitionInitialGraphForceChain(MEX_ARGS){
	if(nrhs < 1){
		mexErrMsgTxt("PartitionInitialGraphForceChain: expected one argument \n");
	}

	Graph& graph = *handle_to_ptr<Graph>(prhs[0]);

	ForestPrecondition *precond;
	precond = new ForestPrecondition;

	precond -> CreateInitialForestForceChain(graph);

	plhs[0] = ptr_to_handle(precond);
}

// Partition the graph by active edges. Used for reconditioning.
// input: graph pointer: the graph contsrtucted before.
// 			  p^k(|E|*1, array), the dual variable at iteration k.
//				use_chain (bool): true: a fatser backward solver will be called later.
// output(two subgraphs pointers): one subgraph has active edges, another has inactive.
static void PartitionGraphActiveSets(MEX_ARGS){
	if(nrhs < 4){
		mexErrMsgTxt("PartitionGraphActiveSets: expected four arguments \n");
	}

	Graph& graph = *handle_to_ptr<Graph>(prhs[0]);

	ForestPrecondition *precond = handle_to_ptr<ForestPrecondition>(prhs[1]);
	
	// the dual variable at k-th interation.
	double *pk = (double*)mxGetPr(prhs[2]);

	// when use_chain is true, an efficient back ward solver will be applied for the chains in the forests.
	bool use_chain = (bool)mxGetScalar(prhs[3]);

	// release the memory of previous trees
	precond-> ClearMemory();

	// create the new trees
	precond-> CreateInitialForest(graph, use_chain, pk);

	//precond->PrintTreeInformation();

	plhs[0] = ptr_to_handle(precond);
}

// return the number of partitioned forests.
static void GetNumberOfForests(MEX_ARGS){
	ForestPrecondition *precond = handle_to_ptr<ForestPrecondition>(prhs[0]);
	plhs[0] = mxCreateDoubleScalar(precond->GetNumberOfForests());
}

// apply backward solver to dual update in PDHG.
// solve p^{k+1} = argmin F*(p)- <Ku_bar, p> + t/2 ||p-p^k||^2
// input: preconditioner pointer, pointer of forest preconditioner
// 			  p^k(|E|*1, array), the dual variable at iteration k.
// 			  u_bar(|V|*1, array), the over-relaxation of primal variable.
// 			  w(|E|*1, array), weight of each edge.
// 			  t(double), step size.
//				a(|V|*1, array), weight in front of data term.
// output: p^{k+1}(|E|*1, array), the dual variable at iteration k+1.
static void ForestBackwardPDHG(MEX_ARGS){
	ForestPrecondition *precond = handle_to_ptr<ForestPrecondition>(prhs[0]);

	size_t number_edges = precond->number_edges();
	size_t number_vertices = precond->number_vertices();

	double *pk = (double*)mxGetPr(prhs[1]);
	double *u_bar = (double*) mxGetPr(prhs[2]);
	double *w = (double*)mxGetPr(prhs[3]);
	double t = (double)mxGetScalar(prhs[4]);
	double *a = (double*)mxGetPr(prhs[5]);

	const mwSize *dims = mxGetDimensions(prhs[1]);
	if (dims[0] != number_edges)
		mexErrMsgTxt("ForestBackwardPDHG: pk should be a matrix of size number_edges * 1");

	//create the result vector (|E|*1) to store p^{k+1}
	plhs[0] = mxCreateDoubleMatrix(number_edges, 1, mxREAL);
	double *result = mxGetPr(plhs[0]);

	precond->BackwardStepPDHG(pk, u_bar, w, t, result, a);
}

// given the graph, return the nabla matrix.
// input: graph pointer, the graph structure in boost library
// output: sparse matrix, a sparse adjacent matrix
static void GraphToNabla(MEX_ARGS){
	Graph *graph = handle_to_ptr<Graph>(prhs[0]);
	const size_t num_edges = boost::num_edges(*graph);
	const size_t num_vertices = boost::num_vertices(*graph);

	plhs[0] = mxCreateSparse(
			num_edges,
			num_vertices,
			num_edges*2,
			mxREAL);

	plhs[1] = mxCreateDoubleMatrix(boost::num_edges(*graph), 1, mxREAL);

	std::vector<double> csr_data;
	std::vector<mwIndex> csr_rowstart;
	std::vector<mwIndex> csr_colidx;

	csr_rowstart.push_back(0);
	
	double *weight_vec = mxGetPr(plhs[1]);
	WeightMap weights = boost::get(boost::edge_weight_t(), *graph);
	boost::graph_traits<Graph>::edge_iterator ei, ei_end;

	for(boost::tie(ei, ei_end) = boost::edges(*graph); ei != ei_end; ei++){
		const int src = boost::source(*ei, *graph);
		const int dst = boost::target(*ei, *graph);

		*weight_vec++ = boost::get(weights, *ei).weight_;
		csr_rowstart.push_back( csr_rowstart.back() + 2);

		csr_data.push_back(1);
		csr_data.push_back(-1);

		csr_colidx.push_back(src);
		csr_colidx.push_back(dst);
	}

	std::vector<double> csc_data; csc_data.resize(num_edges*2);
	std::vector<mwIndex> csc_colstart; csc_colstart.resize(num_vertices + 1);
	std::vector<mwIndex> csc_rowidx; csc_rowidx.resize(num_edges * 2);

	csr2csc<double, mwIndex>(num_edges, num_vertices, num_edges*2,
			(double*)&csr_data[0], (mwIndex*)&csr_colidx[0], (mwIndex*)&csr_rowstart[0],
			(double*)&csc_data[0], (mwIndex*)&csc_rowidx[0], (mwIndex*)&csc_colstart[0]);

	double *mat_pr = (double*)mxGetPr(plhs[0]);
	mwIndex *mat_ir = (mwIndex*)mxGetIr(plhs[0]);
	mwIndex *mat_jc = (mwIndex*)mxGetJc(plhs[0]);

	std::copy(csc_data.begin(), csc_data.end(), mat_pr);
	std::copy(csc_rowidx.begin(), csc_rowidx.end(), mat_ir);
	std::copy(csc_colstart.begin(), csc_colstart.end(), mat_jc);
}

// list of all functions
const static map<string, function<void(MEX_ARGS)>> cmd_reg = {
	{ "CreateGraphFromAdjacency", 			 CreateGraphFromAdjacency },
	{ "CreateGraphFrombkFile",					 CreateGraphFrombkFile},
	{ "GraphToNabla",										 GraphToNabla },
	{ "GetNumberOfForests",							 GetNumberOfForests },

	{ "PartitionInitialGraph", 					 PartitionInitialGraph },
	{ "PartitionGraphActiveSets",				 PartitionGraphActiveSets },
	{ "PartitionInitialGraphForceChain", PartitionInitialGraphForceChain},
	{ "ForestBackwardPDHG",							 ForestBackwardPDHG },	

};

// the gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if (nrhs == 0)
		mexErrMsgTxt("Usage: precond_mex(command, arg1, arg2, ...);");

	char *cmd = mxArrayToString(prhs[0]);
	bool executed = false;

	for(auto&  c :cmd_reg){
		if(c.first.compare(cmd)==0){
			c.second(nlhs, plhs, nrhs-1, prhs+1);
			executed = true;
			break;
		}
	}

	if(!executed){
		stringstream msg;
		msg<<"Unkown command' " <<cmd <<"'. List of supported commands:";
		for(auto&c : cmd_reg)
			msg<<"\n - " << c.first.c_str();

		mexErrMsgTxt(msg.str().c_str());
	}

	mexEvalString("drawnow;");
}
