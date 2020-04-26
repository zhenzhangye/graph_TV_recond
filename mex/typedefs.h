/* Classes definition file 
 *
 * define the edge, vertice classes used in graph
 *
 */

#ifndef __TYPEDEFS_H__
#define __TYPEDEFS_H__
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>

#define MAX_CORE 32
#define THRESHOLD 1e-10

// vertice class used in graph
class VertexData{
private:
	friend class boost::serialization::access;

	template<class Archive>
		void serialize(Archive & ar, const unsigned int version){
			ar & unary_;
		}

public: 
	VertexData(): unary_(0) {}
	// unary variable, storing the data term.
	double unary_; 
};

// edge class used in graph
class EdgeData {
private:
	friend class boost::serialization::access;
	
	template<class Archive>
		void serialize(Archive & ar, const unsigned int version){
			ar & weight_;
			ar & original_index_;
			ar & inactive_weight_;
		}

public:
	EdgeData() : weight_(0), original_index_(0), inactive_weight_(0){ }
	EdgeData(double w, size_t i, double i_w) : weight_(w), original_index_(i), inactive_weight_(i_w) { }

	// weight of the edge in the graph.
	double weight_;

	// the index of this edge in the graph.
	size_t original_index_;

	// this weight is computed by p^k, used for partitioning only.
	double inactive_weight_;

	// define some comparison operators for edges (used for graph partition).
	bool operator>(const EdgeData& rhs) const{
		return inactive_weight_ > rhs.inactive_weight_;
	}

	bool operator>=(const EdgeData& rhs) const{
		return inactive_weight_ >= rhs.inactive_weight_;
	}

	bool operator==(const EdgeData& rhs) const{
		return inactive_weight_ == rhs.inactive_weight_;
	}

	bool operator<=(const EdgeData& rhs) const{
		return inactive_weight_ <= rhs.inactive_weight_;
	}

	bool operator<(const EdgeData& rhs) const{
		return inactive_weight_ < rhs.inactive_weight_;
	}
};

// define some alias for boost graph related types.
typedef boost::property<boost::edge_weight_t, EdgeData> EdgeWeightProperty;
typedef boost::property<boost::edge_weight_t, double> EdgeWeightDouble;
typedef std::pair<double, size_t> iPair;
 
typedef boost::adjacency_list<
		    boost::hash_setS,
				boost::vecS,
				boost::undirectedS,
				VertexData,
				EdgeWeightProperty> Graph;

typedef boost::adjacency_list<
			  boost::vecS,
			  boost::vecS,
			  boost::undirectedS,
			  boost::property<boost::vertex_distance_t, int>,
			  EdgeWeightDouble> GraphPrim;

typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;
typedef boost::property_map<GraphPrim, boost::edge_weight_t>::type WeightMapPrim;


#endif

