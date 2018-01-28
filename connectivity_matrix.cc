#include "connectivity_matrix.h"
#include <iostream>
#include <cassert>
#include <climits>

using namespace std;

/**
 * searches through an array with size @size and finds the minimum element and returns the index 
 * of the minimum value, and stores the acutal value of the minimum in solution.
 **/
int search_min(int array[], int size, int& solution) {
	if (size == 1) {
		solution = array[0];
		return 0;
	}

	int index = -1;
	int current_min = INT_MAX;
	for (int i = 0; i < size; i++) {
		if (array[i] == -1)
			continue;
		if (current_min > array[i] && array[i] != -1) {
			index = i;
			current_min = array[i];
		}
	}
	array[index] = -1;
	return index;
}

void target_groups(int num_groups, 
					int src_group, 
					std::vector<int>& incoming_switches_left, 
					std::vector<int>& groups) {

	int target_groups[num_groups-1];
	int index = 0;
	// given a src_group, forms an array of all
	for (int i = 0; i < num_groups; i++) {
		if (src_group == i) 
			continue;
		target_groups[index] = i;
		index++;
	}


	int size = num_groups - 1;
	int array[size];
	for (int i = 0; i < size; i++) {
		array[i] = incoming_switches_left[target_groups[i]];
	}
	for (int i = 0; i < size; i++) {
		int solution;
		int index = search_min(array, size, solution);
		groups.push_back(target_groups[index]);
	}
}

void generate_butterfly(int num_groups,  std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	int num_optical_switches = num_groups - 1;
	for (int i = 0; i < num_groups; i++) {
		groups[i] = new node(i, num_optical_switches, false);
		assert(groups[i]);
	}
	for (int i = 0; i < num_optical_switches; i++) {
		optical_switches[i] = new node(i, num_groups, true);
		assert(optical_switches[i]);
	}
};

void link_butterfly(int num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	for (int i = 0; i < num_groups; i++) {
		node* curr_group_node = groups[i];
		assert(curr_group_node);
		assert(!curr_group_node->is_optical_switch());
		for (int j = 0; j < num_groups - 1; j++) {
			node* curr_optical_switch = optical_switches[j];
			assert(optical_switches[j]);
			curr_group_node->set_child(j, curr_optical_switch);
			printf(" The pointer for optical_switch: %d is %p\n", j, curr_optical_switch);
		}
	}
	for (int i = 0; i < num_groups - 1; i++) {
		node* current_optical_switch = optical_switches[i];
		assert(current_optical_switch);
		assert(current_optical_switch->is_optical_switch());
		for (int j = 0; j < num_groups; j++) {
			node* current_group = groups[j];
			assert(current_group);
			printf(" The pointer for group: %d is %p\n", j, current_group);
			current_optical_switch->set_child(j, current_group);
		}
	}
};

/**
 * Butterfly is the canonical Dragonfly case in which every group should be able to get 
 * to every other group in the interconnect using just 1 inter-group link. This function 
 * both generates the group nodes and the optical switch nodes and then subsequently links 
 * them together. 
 * NOTE: This function has to always be followed up using a delete_buttefly function
 * 			call in order to free all the nodes which were malloc-ed earlier using form butterfly
 **/
void form_butterfly(int num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	generate_butterfly(num_groups, groups, optical_switches);
	link_butterfly(num_groups, groups,optical_switches);
	for (int i = 0; i < num_groups; i++) {
		node* group = groups[i];
		assert(!group->is_optical_switch());
		std::cout << "current group is: " << group->get_id() << std::endl;
		for (int j = 0; j < num_groups - 1; j++) {
			printf("port: %d is connected to switch: %d with pointer %p\n",j,optical_switches[j]->get_id(),optical_switches[j]);
			//std::cout << "port: " << std::to_string(j) << " connects to switch: " << std::to_string((group->get_child(j))->get_id()) << std::endl;
		} 
	}

	for (int i = 0; i < num_groups - 1; i++) {
		node* optical_switch = optical_switches[i];
		std::cout << "current optical_switch is: " << optical_switch->get_id() << std::endl;
		assert(optical_switch->is_optical_switch());
		for (int j = 0; j < num_groups; j++) {
			printf("port: %d is connected to group: %d with pointer %p\n",j,groups[j]->get_id(),groups[j]);
		}
	}
	//for (int i)
};

void delete_butterfly(int num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	for (int i = 0; i < groups.size(); i++) {
		delete groups[i];
	}	
	for (int i = 0; i < optical_switches.size(); i++) {
		delete optical_switches[i];
	}
};

/**
 * The algorithm that calls the DFS function to search the butterfly graph 
 * that is the Flexfly topology to figure out what the configuration should 
 * be 
 **/
void canonical_dragonfly_config_greedy(int num_groups, std::vector<node*>& groups, 
										std::vector<node*>& optical_switches,  
										std::vector<std::vector<int> >& optical_inout_connections) {
	// for all switches, find the optical switch that gets them to every other switch
	std::vector<int> incoming_switches_left (num_groups);
	std::fill(incoming_switches_left.begin(), incoming_switches_left.end(), num_groups - 1);
	for (int i = 0; i < num_groups; i++) {
		node* curr_group = groups[i];
		assert(curr_group);
		int curr_group_id = curr_group->get_id();
		std::vector<int> in_order_group;
		
		// prioritize the 
		target_groups(num_groups, i, incoming_switches_left, in_order_group);
		for (int j = 0; j < num_groups - 1; j++) {	
			int target_group = in_order_group[j];
			// NOTE: j is the target group's id
			node* optical_switch = dfs(curr_group, 3, target_group);
			assert(curr_group->get_id() != target_group);
			if (optical_switch == nullptr) {
				std::cout << "the failed switch: " << std::to_string(curr_group->get_id()) << " and target group: " << std::to_string(target_group) << std::endl;
			} else {
				incoming_switches_left[target_group]--; // NOTE: keep this in mind
				assert(optical_switch->is_optical_switch());
				std::vector<int>& inout_connection_for_optical_switch = optical_inout_connections[optical_switch->get_id()];
				if (inout_connection_for_optical_switch.size() != num_groups) {
					inout_connection_for_optical_switch.resize(num_groups);
				}
				inout_connection_for_optical_switch[i] = target_group; 
			}
		}
	}
};

/**
 * short for depth-first search
 */ 
node* dfs(node* curr_node, int depth, int target_group_id) {
	if (curr_node == nullptr) 
		return nullptr;
	
	node* found_node = nullptr;
	if (depth == 1) {
		assert(!curr_node->is_optical_switch());
		found_node = (target_group_id == curr_node->get_id()) ? curr_node : nullptr;
	} else {
		//cout << to_string(curr_node->get_num_child()) << endl;
		int iter = curr_node->get_num_child();
		for (int i = 0; i < iter; i++) {
			found_node = dfs(curr_node->get_child(i), depth - 1, target_group_id);
			if (found_node == nullptr)
				continue;
			
				// return curr node so that the group will know which 
				// optical switch can be connected to the destination group
			curr_node->set_child(i, nullptr);
			if (curr_node->is_optical_switch()) {
				assert(depth == 2);
				found_node = curr_node;
			}
			break;
		}
	}
	return found_node;
};


void check_switch_config(std::vector< std::vector <int> >& optical_inout_connections) {
	for (int i = 0; i < optical_inout_connections.size(); i++) {
		std::vector<int>& current_config = optical_inout_connections[i];
		std::cout << "Switch: " << std::to_string(i) << std::endl;
		for (int j = 0; j < current_config.size(); j++) {
			std::cout << "inport: " << std::to_string(j) << " -> outport: " << std::to_string(current_config[j]) << std::endl; 
		}
		std::cout << std::endl;
	}	
}

/**
 * Combines all the above the algorithms above and into one. It begins by generating the physical
 * Butterfly topology connection where middle layer switches are all optical-switches. Next it will
 * run a DFS algorithm to find the canonical configuration (dfly with 1 inter-group link to every group)
 * and finally frees all the nodes that were allocated when calling form_butterfly() function.
 *
 * NOTE: This is the function that should really be called by external codes as it abstracts away a lot of 
 * 			the internal implementation of finding the optical switch inport-outport configuration needed to
 * 			achieve the canonical dragonfly topology.
 **/
void configure_optical_switches_canonical(int num_groups, std::vector< std::vector<int> >& optical_inout_connections) {
	std::vector<node*> groups;
	std::vector<node*> optical_switches;
	groups.resize(num_groups);
	groups.reserve(num_groups);
	optical_switches.resize(num_groups - 1);
	optical_switches.reserve(num_groups - 1);
	// first form (generate and link) butterfly physical connection
	form_butterfly(num_groups, groups, optical_switches);
	optical_inout_connections.resize(num_groups - 1);
	optical_inout_connections.reserve(num_groups - 1);
	// traverse the butterfly graph to find the canonical dfly optical configuration
	canonical_dragonfly_config_greedy(num_groups, groups, optical_switches, optical_inout_connections);
	// finally delete all the groups and optical switch nodes;
	delete_butterfly(num_groups, groups, optical_switches);
}


inline int wrapped_increment(int index, int max_num) {
	return (index + 1) % max_num;
}

/**
 * Forms the default group-level all-to-all connectivity to emulate the Dragonfly topology
 * Also the model assumes that the number of optical switches equals to the number of switches per group
 * and that the radix of each optical switch equals the number of groups in the Flexfly topology
 *
 **/
void configure_default_simpler_model(int num_groups, std::unordered_map<int, std::vector<int>>& optical_inout_connections) {
	int num_optical_switches = num_groups - 1;
	for (int i = 0; i < num_groups; i++) {
		// i equals the current group
		int tmp = i; // stores the current index
		for (int j = 0; j < num_groups - 1; j++) {
			std::vector<int>& current_switch_config = optical_inout_connections[j];
			if (current_switch_config.size() == 0) {
				//current_switch_config.reserve(num_groups);
				current_switch_config.resize(num_groups);
				current_switch_config.reserve(num_groups);
			}
			tmp = wrapped_increment(tmp, num_groups);
			current_switch_config[i] = tmp;
		}
		
	}
}

