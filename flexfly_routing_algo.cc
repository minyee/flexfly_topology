#include "flexfly_routing_algo.h"
#include <sstmac/hardware/topology/structured_topology.h>
#include <queue>

namespace sstmac {
namespace hw{

// Use stack and returns the distance from src_switch to dst_switch
int depth_first_search(std::vector< std::vector<int> >& connectivity_matrix, 
														int src_switch, 
														int dst_switch) {
	int curr_switch = src_switch;
	std::stack<int> switches_stack;
	switches_stack.push(src_switch);
	std::unordered_map<int, int> distance_switches;
	distance_switches.insert({curr_switch,0});
	bool found = false;
	while(!found && !switches_stack.empty()) {
		curr_switch = switches_stack.top();
		switches_stack.pop();
		auto dist_iter_curr_switch = distance_switches.find(curr_switch);
		int distance_to_curr_switch = dist_iter_curr_switch->second;
		std::vector<int>& connectivity_vector = connectivity_matrix[curr_switch];
		for (int i = 0; i < connectivity_vector.size(); i++) {
			int target_switch = connectivity_vector[i];
			auto searched_node = distance_switches.find(target_switch);
			// if node has never been visited before
			if (searched_node != distance_switches.end()) {
				int distance_to_target_switch = searched_node->second;
				searched_node->second = std::min(distance_to_curr_switch + 1, distance_to_target_switch);
			} else {
				// push node into stack if 
				switches_stack.push(target_switch);
				distance_switches.insert({target_switch, distance_to_curr_switch + 1});
			}
		}
	}
	auto dist_iter_dst_switch = distance_switches.find(dst_switch);
	assert(dist_iter_dst_switch != distance_switches.end());
	return dist_iter_dst_switch->second;
}; 

// use queue and returns the distance from src_switch to dst_switch
int breadth_first_search(std::vector< std::vector<int> >& connectivity_matrix,
														int src_switch, 
														int dst_switch) {
	int curr_switch = src_switch;
	std::queue<int> switches_queue;
	switches_queue.push(src_switch);
	std::unordered_map<int, int> distance_switches;
	distance_switches.insert({curr_switch,0}) ;
	bool found = false;
	while(!found && !switches_queue.empty()) {
		curr_switch = switches_queue.front();
		switches_queue.pop();
		auto dist_iter_curr_switch = distance_switches.find(curr_switch);
		int distance_to_curr_switch = dist_iter_curr_switch->second;
		std::vector<int>& connectivity_vector = connectivity_matrix[curr_switch];
		for (int i = 0; i < connectivity_vector.size(); i++) {
			int target_switch = connectivity_vector[i];
			auto searched_node = distance_switches.find(target_switch);
			// if node has never been visited before
			if (searched_node != distance_switches.end()) {
				int distance_to_target_switch = searched_node->second;
				searched_node->second = std::min(distance_to_curr_switch + 1, distance_to_target_switch);
			} else {
				// push node into stack if 
				switches_queue.push(target_switch);
				distance_switches.insert({target_switch, distance_to_curr_switch + 1});
			}
		}
	}
	auto dist_iter_dst_switch = distance_switches.find(dst_switch);
	assert(dist_iter_dst_switch != distance_switches.end());
	return dist_iter_dst_switch->second;
};

bool is_optical_switch(int i, int num_groups, int switches_per_group) {
	return i >= (num_groups * switches_per_group);
};

/**
 * Bellman ford routing scheme to routing entire topology 
 * physical_connection is the matrix with entry i,j meaning the number of links connecting switch i
 * to switch j.
 */
void bf_routing_minimal(int src, int dst, 
							std::vector<std::vector<int>>& topology_physical_connection,
							std::vector<std::vector<int>>& optical_inout_configuration,
							int num_groups, 
							int switches_per_group) {
	int total_switches = topology_physical_connection.size();
	// START WITH SOME PREPROCESSING, CHANGING THE 
	// inout port connection to a topology adjacency matrix
	std::vector<std::vector<int>> adjacency_matrix;
	adjacency_matrix.resize(total_switches); 
	for (int j = 0; j < total_switches; j++) {
		adjacency_matrix[j].resize(total_switches);
	}
	for (int i = 0; i < total_switches; i++) {
		std::vector<int>& curr_switch_outgoing_links = topology_physical_connection[i];
		for (int j = 0; i < total_switches; j++) {
			adjacency_matrix[i][j]++;
		}
	}
	// DONE WITH PREPROCESSING
	bool converged = false;
	std::vector<std::vector<distance_optical_id_pair>> distance_matrix;

	// initialize the distance matrix
	for (int i = 0; i < total_switches; i++) {
		for (int j = 0; i < total_switches; j++) {
			if (adjacency_matrix[i][j] > 0) {
				distance_matrix[i][j].distance = 1;
				//distance_vector[i][j].optical_switch = -1;	
			} else {
				distance_matrix[i][j].distance = INT_MAX;
			}
			distance_matrix[i][j].optical_switch = -1;
		}
	}

	// start
	while (!converged) {
		// NOTE: i - src_switch, j - dst_switch
		// first iterate through all the switches
		bool changed = false;
		for (int i = 0; i < total_switches; i++) { // first toggle through all the switches
			bool is_curr_switch_optical = is_optical_switch(i, num_groups, switches_per_group);
			for (auto neighbor : topology_physical_connection[i]){
				for (int j = 0; j < total_switches; j++) {
					bool is_target_switch_optical = is_optical_switch(j, num_groups, switches_per_group);
					//distance_matrix[i][j] = min(distance_matrix[i][j], );
				}
			}
		}

		converged = !changed;
	}
}

/**
 * This is a subroutine called by the dijkstra's routing algorithm and not called
 * by anything else.
 **/
void add_to_path(std::vector<flexfly_path *> &path_collection, 
				int dst_switch, 
				flexfly_path* basis,
				std::vector<int> &outports,
				std::vector<int> &switch_ids,
				bool useRand) {
	uint32_t seed = 100;
	std::srand(seed);
	int rand_num = rand() % 2;
	assert(outports.size() == switch_ids.size());
	flexfly_path* new_path = new flexfly_path;
	if (!path_collection[dst_switch]) {
		
		new_path->path_length = outports.size();
		if (basis) {
			new_path->path_length += basis->path_length;
			new_path->path.resize(new_path->path_length);
			for (int i = 0; i < new_path->path_length; i++) {
				switch_port_pair *spp = new switch_port_pair();
				if (i < basis->path.size()) {
					spp->switch_id = basis->path[i]->switch_id;
					spp->outport = basis->path[i]->outport;
				} else {
					spp->switch_id = switch_ids[i - basis->path_length];
					spp->outport = outports[i - basis->path_length];
				}
				new_path->path[i] = spp;
			}
		} else {
			new_path->path.resize(new_path->path_length);
			for (int i = 0; i < new_path->path_length; i++) {
				switch_port_pair *spp = new switch_port_pair();
				spp->switch_id = switch_ids[i];
				spp->outport = outports[i];
				new_path->path[i] = spp;
			}
		}
	} else {
		assert(basis);
		if (!useRand && rand_num == 0) {
			// don't do anything in this case
			return;
		} 
		new_path->path_length = basis->path_length + outports.size();
		new_path->path.resize(new_path->path_length);
		for (int i = 0; i < new_path->path_length; i++) {
			switch_port_pair *spp = new switch_port_pair();
			if (i < basis->path.size()) {
				spp->switch_id = basis->path[i]->switch_id;
				spp->outport = basis->path[i]->outport;
				//switch_port_pair *spp = new switch_port_pair();
			} else {
				//switch_port_pair *spp = new switch_port_pair();
				spp->switch_id = switch_ids[i - basis->path_length];
				spp->outport = outports[i - basis->path_length];
				
			}
			new_path->path[i] = spp;
		}
	}
	assert(new_path->path.size() == new_path->path_length);
	path_collection[dst_switch] = new_path;
	return;
};

/**
 * Deletes tha path information from source switch to dest switch.
 * This is a subroutine called by the Dijkstra_minimal_routing algorithm
 **/
void clear_path_collections(std::vector<flexfly_path *> &path_collection, 
								int dst_switch) {
	std::cout << "clear_path_colletions in flexfly_routing_algo" << std::endl;
	
	if (path_collection[dst_switch] != nullptr) {
		//assert
		std::cout << "clear_path_colletions in flexfly_routing_algo1" << std::endl;
		//assert(path_collection[dst_switch]->path_length == path_collection[dst_switch]->path.size());
		int path_len = path_collection[dst_switch]->path_length;
		assert(path_len == path_collection[dst_switch]->path.size());
		for (int i = 0; i < path_collection[dst_switch]->path.size(); i++) {
			std::cout << "clear_path_colletions in flexfly_routing_algo2" << std::endl;
			if (path_collection[dst_switch]->path[i] != nullptr)
				delete path_collection[dst_switch]->path[i];
		}
		
		path_collection[dst_switch]->path.resize(0);
		path_collection[dst_switch]->path_length = 0;

		delete path_collection[dst_switch];
		path_collection[dst_switch] = nullptr;
	}

};

void print_switch_link(switch_link* sl) {
	if (sl) {
		std::string link_type = "";
		if (sl->type == flexfly_topology::Link_Type::electrical)
			link_type += "ELECTRICAL";
		else 
			link_type += "OPTICAL";

		std::string printout = "";
		printout += ("src_id: " + std::to_string(sl->src_sid) + "\n");
		printout += ("dst_id: " + std::to_string(sl->dest_sid) + "\n");
		printout += ("src_outport " + std::to_string(sl->src_outport) + "\n");
		printout += ("dst_inport " + std::to_string(sl->dest_inport) + "\n");
		printout += ("Link type: " + link_type);
		std::cout << printout << std::endl;
		
	}
}

/**
 * Create 
 **/
void dijkstra_minimal_route(int src, 
				const std::unordered_map<switch_id, std::vector<switch_link *>>& topology_physical_connection_outport,
				const std::unordered_map<switch_id, std::vector<switch_link *>>& topology_physical_connection_inport,
				std::vector<int> &distance_vector, // should this be distance vector instead?
				std::vector< flexfly_path * > &path_collection,
				std::unordered_map<int, std::vector<int> > &optical_inout_configuration, 
				int num_groups,
				int switches_per_group) {
	int total_switches = num_groups * switches_per_group;//topology_physical_connection_inport.size();
	std::unordered_map<int,bool> visited_nodes;
	//std::cout << "The size of inport vector is: " + std::to_string(topology_physical_connection_inport.size()) << std::endl;
	for (int i = 0; i < total_switches; i++) {
		visited_nodes[i] = false;
	}
	distance_vector.resize(total_switches);
	for (int i = 0; i < total_switches; i++) {
		distance_vector[i] = 10000;
		path_collection[i] = nullptr;
	}
	distance_vector[src] = 0; 
	visited_nodes[src] = true;
	std::queue<int> switch_queue;
	switch_queue.push(src);
	int curr_switch;
	std::cout << "In Dijkstra's minimal routing for switch: " + std::to_string(src) << std::endl;
	while (!switch_queue.empty()) {
		curr_switch = switch_queue.front();
		switch_queue.pop();
		visited_nodes[curr_switch] = true;



		auto outport_connections = topology_physical_connection_outport.find(curr_switch);
		
		if (outport_connections != topology_physical_connection_outport.end()) {
			const std::vector<switch_link *>& curr_switch_outport_connection = outport_connections->second;
			int outport_count = curr_switch_outport_connection.size();
			//std::cout << "The outport is " + std::to_string(outport_count) << std::endl; 
			for (int i = 0; i < outport_count; i++) {
				std::vector<int> outports;
				std::vector<int> switch_ids;
		
				switch_ids.insert(switch_ids.end(), curr_switch);
				outports.insert(outports.end(), i);
				
				//outports[0] = i;
				// CASE 1: neighbor happens to be an Optical Switch
				//print_switch_link(curr_switch_outport_connection[i]);
				if (is_optical_switch(curr_switch_outport_connection[i]->dest_sid, 
										num_groups, 
										switches_per_group)) {
					int optical_switch = curr_switch_outport_connection[i]->dest_sid;

					auto optical_incoming_links_iter = topology_physical_connection_inport.find(optical_switch);
					const std::vector<switch_link *>& incoming_links = optical_incoming_links_iter->second;
					int optical_switch_radix = incoming_links.size();
					//std::cout << "Optical switch has " + std::to_string(optical_switch_radix) + " outports" << std::endl;
					int optical_switch_inport = 0;
					for (optical_switch_inport = 0; optical_switch_inport < optical_switch_radix; optical_switch_inport++) {
						

						switch_link *curr_incoming_link = incoming_links[optical_switch_inport];

						

						int src_switch = curr_incoming_link->src_sid;

						

						if (src_switch == curr_switch) {
							optical_switch_inport = optical_switch_inport;
							break;
						}
					}
					optical_switch = optical_switch - num_groups*switches_per_group;
					//std::cout << "Optical switch has value of : " << std::to_string(optical_switch) << std::endl;
					int optical_switch_outport = optical_inout_configuration[optical_switch][optical_switch_inport];
					optical_switch = optical_switch + num_groups*switches_per_group;
					auto optical_switch_outport_iter = topology_physical_connection_outport.find(optical_switch);
					const switch_link* outgoing_link = optical_switch_outport_iter->second[optical_switch_outport];
					//assert(outgoing_link->src_sid == optical_switch); // the src switch of this outgoing link has to be the same 
																	// as the id of the optical switch
					int dst_elec_switch = outgoing_link->dest_sid;
					//for (int j = 0; j < optical_switch_radix; j++) {
						// one distance is due to optical switch
					//outports
					std::cout << std::to_string(dst_elec_switch) << std::endl;
					if (distance_vector[curr_switch] + 2 < distance_vector[dst_elec_switch]) {
						outports.insert(outports.end(), optical_switch_outport);
						switch_ids.insert(switch_ids.end(), optical_switch);
						distance_vector[dst_elec_switch] = distance_vector[curr_switch] + 2;
						clear_path_collections(path_collection, dst_elec_switch);
						add_to_path(path_collection, dst_elec_switch, path_collection[curr_switch], outports, switch_ids, false);
						if (!visited_nodes[dst_elec_switch])
							switch_queue.push(dst_elec_switch);
					} else if (distance_vector[curr_switch] + 2 == distance_vector[dst_elec_switch]) {
						outports.insert(outports.end(), optical_switch_outport);
						switch_ids.insert(switch_ids.end(), optical_switch);
						add_to_path(path_collection, dst_elec_switch, path_collection[curr_switch], outports, switch_ids, true);
						if (!visited_nodes[dst_elec_switch])
							switch_queue.push(dst_elec_switch);
					} 
					
					//NOTE: Don't push the optical switch into the not visited stack because they work differently.
					//if (!visited_nodes[dst_elec_switch])
						//switch_queue.push(dst_elec_switch);
				} else {
					//auto elec_switch_outport_iter = topology_physical_connection_outport.find()
					const int dst_elec_switch = curr_switch_outport_connection[i]->dest_sid;
					// CASE 2: neighbor is an Electrical Switch
				
					if (distance_vector[curr_switch] + 1 < distance_vector[dst_elec_switch]) {
						distance_vector[dst_elec_switch] = distance_vector[curr_switch] + 1;
						clear_path_collections(path_collection, dst_elec_switch);
						add_to_path(path_collection, dst_elec_switch, path_collection[curr_switch], outports, switch_ids, false);
						if (!visited_nodes[dst_elec_switch])
							switch_queue.push(dst_elec_switch);
					} else if (distance_vector[dst_elec_switch] == distance_vector[curr_switch] + 1) {
						add_to_path(path_collection, dst_elec_switch, path_collection[curr_switch], outports, switch_ids, true);
						if (!visited_nodes[dst_elec_switch])
							switch_queue.push(dst_elec_switch);
					}

					// if the neighbor is not optical switch and have yet been visited, then push it into stack
					
				}
			}
		} 
		//visited_nodes[curr_switch] = true;
	}
	std::cout << "The distance vector: " << std::endl;
	for (int i = 0; i < total_switches; i++) {
		std::cout << " | " + std::to_string(distance_vector[i]);
	}
	std::cout << std::endl;
	//std::cout << "\n\n\n" << std::endl;	
};




}
}
