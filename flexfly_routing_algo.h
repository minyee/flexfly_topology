#include <vector>
#include <stack>
#include <queue>
#include <cassert>
#include <utility>
#include <climits>
#include <cstdlib>
#include <unordered_map>
#include "flexfly_topology.h"

namespace sstmac {
namespace hw{
// Use stack and returns the distance from src_switch to dst_switch
int depth_first_search(std::vector< std::vector<int> >& connectivity_matrix, 
														int src_switch, 
														int dst_switch); 

// use queue and returns the distance from src_switch to dst_switch
int breadth_first_search(std::vector< std::vector<int> >& connectivity_matrix,
														int src_switch, 
														int dst_switch);

struct distance_optical_id_pair {
	int distance;
	int optical_switch; // through which switch
};

bool is_optical_switch(int i, int num_groups, int switches_per_group);

/**
 * Bellman ford routing scheme to routing entire topology 
 * physical_connection is the matrix with entry i,j meaning the number of links connecting switch i
 * to switch j.
 */
void bf_routing_minimal(int src, int dst, 
							std::vector<std::vector<int>>& topology_physical_connection,
							std::vector<std::vector<int>>& optical_inout_configuration,
							int num_groups, 
							int switches_per_group);

/**
 * This is a subroutine called by the dijkstra's routing algorithm and not called
 * by anything else.
 * TODO: Finish up add path tonight
 **/
void add_to_path(std::vector<flexfly_path *> &path_collection, 
				int dst_switch, 
				flexfly_path* basis,
				std::vector<int> &outports,
				std::vector<int> &switch_ids,
				bool useRand);

/**
 * Deletes tha path information from source switch to dest switch.
 * This is a subroutine called by the Dijkstra_minimal_routing algorithm
 **/
void clear_path_collections(std::vector<flexfly_path *> &path_collection, 
								int dst_switch);

/**
 * Create 
 **/
void dijkstra_minimal_route(int src, 
				const std::unordered_map<switch_id, std::vector<switch_link *>> &topology_physical_connection_outport,
				const std::unordered_map<switch_id, std::vector<switch_link *>> &topology_physical_connection_inport,
				std::vector<int> &distance_vector, // should this be distance vector instead?
				std::vector< flexfly_path * > &path_collection,
				std::unordered_map<int, std::vector<int> > &optical_inout_configuration, 
				int num_groups,
				int switches_per_group);




}
}
