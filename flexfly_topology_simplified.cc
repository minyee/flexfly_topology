/**
 * Author: Min Yee Teh
 */
#include <iostream>	
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <sstmac/common/node_address.h>
#include <sprockit/sim_parameters.h>
#include <stdlib.h>
#include <sstmac/hardware/topology/topology.h>
#include "flexfly_topology_simplified.h"
#include "connectivity_matrix.h"
#include "flexfly_routing_algo.h"
#include "link_stealing.h"

namespace sstmac {
namespace hw {

 // NOTE: the switch with the same id as the max_switch_id_ is the optical switch.
 flexfly_topology_simplified::flexfly_topology_simplified(sprockit::sim_parameters* params) : 
                              structured_topology(params,InitMaxPortsIntra::I_Remembered, 
                                                  InitGeomEjectID::I_Remembered) {
                                //std::cout << "1 " << std::endl;
  num_groups_ = params->get_int_param("groups"); // controls g
 	switches_per_group_ = params->get_int_param("switches_per_group"); // controls a
  switches_per_group_ = num_groups_ - 1;
 	nodes_per_switch_ = params->get_int_param("nodes_per_switch");
  optical_switch_radix_ = params->get_optional_int_param("optical_switch_radix", switches_per_group_);
  is_simplified_model_ = true;
 	num_optical_switches_ = 1;
  optical_switch_radix_ = num_groups_;
  max_switch_id_ = num_optical_switches_ + (num_groups_ * switches_per_group_) - 1;
  group_connectivity_matrix_.resize(num_groups_);
  for (int i = 0; i < num_groups_; i++) {
    group_connectivity_matrix_[i].resize(num_groups_);
    for (int j = 0; j < num_groups_; j++) {
      if (i == j) {
        group_connectivity_matrix_[i][j] = 0;
      } else {
        group_connectivity_matrix_[i][j] = 1;
      }
    }
  }

  routing_table_.resize(num_groups_ * switches_per_group_);
  distance_matrix_.resize(num_groups_ * switches_per_group_);
  routing_table_2_.resize(num_groups_ * switches_per_group_);
  for (int i = 0; i < num_groups_ * switches_per_group_; i++) {
    distance_matrix_[i].resize(num_groups_ * switches_per_group_);
    routing_table_[i].resize(num_groups_ * switches_per_group_);
    routing_table_2_[i].resize(num_groups_ * switches_per_group_);
    
  }
 	setup_flexfly_topology_simplified();

  /*
   * link_stealking algorithm is applied here to group_connectivity_matrix;
   * begin
   */
  if (params->get_bool_param("link_steal")) {
    matrix_float float_mat;
    std::string filename = "lulesh_512r_32g_res.txt";

    read_float_matrix(filename, float_mat);
    simple_rounding_down(float_mat, group_connectivity_matrix_);
    std::cout << "hello from flexfly simplified topology" << std::endl;
    vector_int row_constraints;
    vector_int col_constraints;
    row_constraints.resize(num_groups_);
    col_constraints.resize(num_groups_);
    std::fill(row_constraints.begin(), row_constraints.end(), num_groups_ - 1);
    std::fill(col_constraints.begin(), col_constraints.end(), num_groups_ - 1);
    massage_random(group_connectivity_matrix_, row_constraints, col_constraints);
  }
  /*
   * end
   */

  route_topology_minimal();
 }

 flexfly_topology_simplified::~flexfly_topology_simplified() {
   for (const std::pair<switch_id, std::vector<switch_link*>> elem : switch_outport_connection_map_) {
     const std::vector<switch_link*>& conn_vector = elem.second;  
     for (auto const&  switch_link_ptr : conn_vector) {
        if (switch_link_ptr)
          free(switch_link_ptr);
     }
   }
   //std::cout << "2 " << std::endl;
   std::cout << "Deconstructor is called" << std::endl;
   for (int i = 0; i < num_groups_ * switches_per_group_; i++) {
    for (int j = 0; j < num_groups_ * switches_per_group_; j++) {
      clear_path_collections(routing_table_[i], j);
    }
   }
 }

 /**
  * the implementation is the one where the optical radix is the same as the # of 
  * groups, and there will be switches_per_group_ number of optical switches
  */
 void flexfly_topology_simplified::setup_flexfly_topology_simplified() {
   // first step: connect all the switches within the same group together
  //std::cout << "3 " << std::endl;
  int last_used_outport[max_switch_id_ + 1];
  int last_used_inport[max_switch_id_ + 1];
  std::memset(&last_used_outport, 0, sizeof(int) * (max_switch_id_ + 1));
  std::memset(&last_used_inport, 0, sizeof(int) * (max_switch_id_ + 1));
   for (int group = 0; group < num_groups_; group++) {
     switch_id group_offset = group * switches_per_group_;
     for (int index = 0; index < switches_per_group_ - 1; index++) {
       switch_id swid = group_offset + index;
       for (int target_index = index + 1; target_index < switches_per_group_; target_index++) {
         switch_id target_swid = group_offset + target_index;
         connect_switches(swid, last_used_outport[swid], target_swid, last_used_inport[target_swid], Electrical);
         last_used_outport[swid]++;
         last_used_inport[target_swid]++;
         connect_switches(target_swid, last_used_outport[target_swid], swid, last_used_inport[swid], Electrical);
         last_used_outport[target_swid]++;
         last_used_inport[swid]++;
       }
     }  
   }

   // second step: connect all the switches within groups to the optical switches
   switch_id optical_swid = max_switch_id_;
   for (int group = 0; group < num_groups_; group++) {
     switch_id group_offset = group * switches_per_group_;
     for (int index = 0; index < switches_per_group_; index++) {
       switch_id swid = group_offset + index;
       assert(last_used_outport[optical_swid] == swid);
       assert(last_used_inport[optical_swid] == swid);
       assert(last_used_outport[swid]==switches_per_group_ - 1);
       assert(last_used_inport[swid]==switches_per_group_ - 1);
       connect_switches(swid, last_used_outport[swid], optical_swid, last_used_inport[optical_swid], Optical);
       
       last_used_outport[swid]++;
       last_used_inport[optical_swid]++;

       connect_switches(optical_swid, last_used_outport[optical_swid], swid, last_used_inport[swid], Optical);
       last_used_outport[optical_swid]++;
       last_used_inport[swid]++;
     }
   }
 }

 /**
  * Called only once at the start of the constructor. Finds all the paths 
  * from 1 switch to the other given the initial optical switch inout port
  * connectivity configurations. Uses minimal shortest path routing for this.
  **/
 void flexfly_topology_simplified::setup_routing_table() {
  //std::cout << "4 " << std::endl;
  int total_switches = num_optical_switches_ + num_groups_ * switches_per_group_;
  if (routing_table_.size() == 0) 
    routing_table_.resize(total_switches);
  for (int i = 0; i < total_switches; i++) {
    routing_table_.resize(total_switches);
  }
 }

 /**
  * IMPORTANT: This function will route the minimal path
  **/
 void flexfly_topology_simplified::minimal_route_to_switch(switch_id src_switch_addr, 
 												switch_id dst_switch_addr, 
 												routable::path& path) const {
  //std::cout << "6 " << std::endl;
  int src_group = src_switch_addr/switches_per_group_;
  int dst_group = src_switch_addr/switches_per_group_;
  // both the source and dest switch cannot be the same, because the caller 
  // should not call this function if this were the case
  assert(src_switch_addr != dst_switch_addr);
  // if both nodes are connected to the same switch
  flexfly_path* fpath = routing_table_[src_switch_addr][dst_switch_addr];
  assert(fpath != nullptr);
  switch_port_pair* spp = fpath->path[0];
  assert(spp->switch_id == src_switch_addr);
  path.set_global_outport(spp->outport);
 };

 int flexfly_topology_simplified::get_output_port(int src_switch, int dst_switch) const {
  //std::cout << "7 " << std::endl;
  auto sl_iter = switch_outport_connection_map_.find(src_switch);
  assert(sl_iter != switch_outport_connection_map_.end());
  auto slv = sl_iter->second;
  //int i = 0;
  for (auto sl : slv) {
    if (sl->dest_sid == dst_switch)
      return sl->src_outport;
    //i++; 
  }
  spkt_abort_printf("Get output port function in flexfly_topology_simplified has failed ");
  return -1;
 };

 /**
  * checks if a given switch id in the topology is an optical switch or not
  */
 bool flexfly_topology_simplified::is_optical_switch(switch_id swid) const {
  //std::cout << "8 " << std::endl;
  if (!valid_switch_id(swid)) {
 		return false;
 	} 
  //std::cout << "got out of 8" << std::endl;
 	return swid == max_switch_id_;
 };

 /**
  * checks if a given switch id in the topology is an electrical switch or not
  */
 bool flexfly_topology_simplified::is_electrical_switch(switch_id swid) const {
  //std::cout << "9" << std::endl;
  return valid_switch_id(swid) && !is_optical_switch(swid);
 };

 // DONE (RECHECK)
 void flexfly_topology_simplified::configure_optical_or_electrical_port_params(switch_id swid, 
                                                                    std::string& str, 
                                                                    sprockit::sim_parameters* switch_params) const {
 	//std::cout << "10 " << std::endl;
  if (!switch_params) {
 		return;
 	}
 	// refers to either optical or electrical
  sprockit::sim_parameters* link_params = switch_params->get_namespace("link");
 	double bandwidth = switch_params->get_optional_bandwidth_param("bandwidth",1000); // in units of bytes/sec I think
 	long buf_space = switch_params->get_optional_byte_length_param("buffer_size", 1024);
  int switch_link_redundancy = switch_params->get_optional_int_param("switch_link_redundancy", 1);
  int credits = ((int) buf_space)*switch_link_redundancy;
  int port_count = is_optical_switch(swid) ? this->num_groups_ : this->switches_per_group_ + this->nodes_per_switch_;
 	for (int i = 0; i < port_count; i++) {
 		topology::setup_port_params(i, credits, bandwidth, link_params, switch_params);
 	}
  //std::cout << "got out of 10" << std::endl;
 };

 // DONE
 void flexfly_topology_simplified::configure_individual_port_params(switch_id src,
          												 sprockit::sim_parameters* switch_params) const {
  //std::cout << "11 " << std::endl;
 	if (valid_switch_id(src)) {
 		std::string str;
 		if (is_optical_switch(src)) {
 			str = "optical";
 		} else {
 			str = "electrical";
 		}
 		configure_optical_or_electrical_port_params(src, str, switch_params);
 	} 	
 };

 /**
  * @Brief Given a source switch and a destination switch, connects the source switch with the 
  * dest switch
  * NOTE: This member function should form a single switch_link
  */
  void flexfly_topology_simplified::connect_switches(switch_id src, int src_outport, switch_id dst, int dst_inport, Link_Type ltype) {
    std::vector<switch_link*>& src_outport_connection_vector = switch_outport_connection_map_[src];
    std::vector<switch_link*>& dst_inport_connection_vector = switch_inport_connection_map_[dst];
    switch_link *conns = new switch_link();
    conns->src_sid = src;
    conns->dest_sid = dst;
    conns->dest_inport = dst_inport; 
    conns->src_outport = src_outport;
    conns->type = ltype;
    src_outport_connection_vector.push_back(conns);
    dst_inport_connection_vector.push_back(conns);
    assert(conns->src_outport == src_outport_connection_vector.size() - 1);
    assert(conns->dest_inport == dst_inport_connection_vector.size() - 1);
    return;
 }


void flexfly_topology_simplified::connected_outports(const switch_id src, 
                                            std::vector<topology::connection>& conns) const {
  //std::cout << "13" << std::endl;
  std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator got = switch_outport_connection_map_.find(src);
  int cidx = 0;
  if (got != switch_outport_connection_map_.end()) {
    const std::vector<switch_link*>& switch_link_vectors = got->second;
    for (switch_link* current_switch_link : switch_link_vectors) {
      conns.push_back(topology::connection());
      conns[cidx].src = src;
      conns[cidx].dst = current_switch_link->dest_sid;
      conns[cidx].src_outport = current_switch_link->src_outport; 
      conns[cidx].dst_inport = current_switch_link->dest_inport;
      //conns[cidx].link_type = current_switch_link->type;
      cidx++;
    }
  }
  //std::cout << "got out of 13" << std::endl;
  //std::cout << "cidx is: " << std::to_string(cidx) << std::endl;
}


bool flexfly_topology_simplified::switch_id_slot_filled(switch_id sid) const {
  //std::cout << "14 " << std::endl;
  return (sid <= max_switch_id_);
}

 void flexfly_topology_simplified::configure_vc_routing(std::map<routing::algorithm_t, int>& m) const {
  //std::cout << "15 " << std::endl;
  m.insert({routing::minimal, 3});
  m.insert({routing::minimal_adaptive, 3});
  m.insert({routing::valiant, 3});
  m.insert({routing::ugal, 3});
  return;
 };

  switch_id flexfly_topology_simplified::node_to_ejection_switch(node_id addr, uint16_t& port) const {
    //std::cout << "16 " << std::endl;
    switch_id swid = addr / nodes_per_switch_; // this gives us the switch id of the switch node addr is connected to
    std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_iter = switch_outport_connection_map_.find(swid);
    const std::vector<switch_link*>& conn_vector = tmp_iter->second;
    port = std::max((int) (conn_vector.size() - 1), 0) + ((int) swid) * nodes_per_switch_; // CHECK THIS AGAIN
    return swid;
  };
  
  switch_id flexfly_topology_simplified::node_to_injection_switch(node_id addr, uint16_t& port) const {
    //std::cout << "17 " << std::endl;
    return node_to_ejection_switch(addr, port);
  };

  /**
   * NOTE: This method does not include the hop to an optical switch
   **/
  int flexfly_topology_simplified::minimal_distance(switch_id src, switch_id dst) const {
    //std::cout << "src switch id: " << std::to_string(src) << " dst switch id: " << std::to_string(dst) << std::endl;
    //std::cout << "18 " << std::endl;
    int src_group = group_from_swid(src);
    int dst_group = group_from_swid(dst);
    if (src == dst) { // same switch
      return 0;
    } else if (src_group == dst_group) { // same group
      return 1;
    } else { // different group but can reach either by 1 global and 1 local or 1 local and then 1 global
      std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_iter = switch_outport_connection_map_.find(src);
      const std::vector<switch_link*>& conn_vector = tmp_iter->second;
      bool two_or_three = true; 
      
      for (switch_link* tmp_link : conn_vector ) {
        if (tmp_link->type == Electrical)
          continue;
      }
      return two_or_three ? 2 : 3;
    }
  };

  int flexfly_topology_simplified::num_hops_to_node(node_id src, node_id dst) const {
    //std::cout << "19 " << std::endl;
    int src_swid = src / (nodes_per_switch_);
    int dst_swid = dst / (nodes_per_switch_);
    int min_dist = distance_matrix_[src_swid][dst_swid];
    return min_dist + 2; // added by 2 because each node is 1 hop away from it's switch
  };

  void flexfly_topology_simplified::nodes_connected_to_injection_switch(switch_id swid, 
                                                              std::vector<injection_port>& nodes) const {
    //std::cout << "20 " << std::endl;
    int i = 0;
    int port_offset = switches_per_group_;
    nodes.resize(nodes_per_switch_);
    for (int i = 0; i < nodes_per_switch_; i++) {
      int node_id = swid * nodes_per_switch_ + i;
      int port_ind = port_offset + i;
      nodes[i].nid = node_id;
      nodes[i].port = port_ind;
    }
    //std::cout << "did i get out " << std::endl;
  };

  void flexfly_topology_simplified::nodes_connected_to_ejection_switch(switch_id swid, 
                                                              std::vector<injection_port>& nodes) const { 
    flexfly_topology_simplified::nodes_connected_to_injection_switch(swid, nodes);
  };

  // returns the group id of a given switch
  int flexfly_topology_simplified::group_from_swid(switch_id swid) const {
    return swid / (switches_per_group_);
  };

  /**
   * prints out the all the connections for each switch
   */
  void flexfly_topology_simplified::print_port_connection_for_switch(switch_id swid) const {
    //std::cout << "1000" << std::endl;
    std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_iter = 
                                        switch_outport_connection_map_.find(swid);
    
    if (tmp_iter == switch_outport_connection_map_.end()) {
      return;
    }

    const std::vector<switch_link*>& connection_vector = tmp_iter->second;

    std::string message;
    int i = 0;
    for (switch_link* sl_ptr : connection_vector) {
      // check if null, if null have to throw an error 
      if (sl_ptr) {
        message += ("Dest switch_id: " + std::to_string(sl_ptr->dest_sid));
        message += (" Dest inport: " + std::to_string(sl_ptr->dest_inport));
        message += " Link type: ";
        if (sl_ptr->type == Electrical) {
          message += ("ELECTRICAL\n");
        } else {
          message += ("OPTICAL\n");
        }
      } else {
        spkt_abort_printf("A switch link with swid: %d is null\n", swid);
      }
      i++;
    }
    std::cout << message; 
  };


  /**
   * @brief num_endpoints To be distinguished slightly from nodes.
   * Multiple nodes can be grouped together with a netlink.  The netlink
   * is then the network endpoint that injects to the switch topology
   * @return
   */
  int flexfly_topology_simplified::num_netlinks() const {
    //std::cout << "22 " << std::endl;
    //std::cout << "num_netlinks?" << std::endl;
    return 1;
  }; 

  switch_id flexfly_topology_simplified::max_netlink_id() const {
    //std::cout << "23 " << std::endl;
    //std::cout << "max_netlink_id?" << std::endl;
    return max_switch_id_;
  };

  bool flexfly_topology_simplified::netlink_id_slot_filled(node_id nid) const {
    //std::cout << "24 " << std::endl;
    //std::cout << "netlink_id_slot_filled?" << std::endl;
    return true;
  };

  /**
     For a given node, determine the injection switch
     All messages from this node inject into the network
     through this switch
     @param nodeaddr The node to inject to
     @param switch_port [inout] The port on the switch the node injects on
     @return The switch that injects from the node
  */
  switch_id flexfly_topology_simplified::netlink_to_injection_switch(
        netlink_id nodeaddr, uint16_t& switch_port) const {
    //std::cout << "25 " << std::endl;
    return max_switch_id_;
  };

  /**
     For a given node, determine the ejection switch
     All messages to this node eject into the network
     through this switch
     @param nodeaddr The node to eject from
     @param switch_port [inout] The port on the switch the node ejects on
     @return The switch that ejects` into the node
  */
  switch_id flexfly_topology_simplified::netlink_to_ejection_switch(
        netlink_id nodeaddr, uint16_t& switch_port) const {
    //std::cout << "26 " << std::endl;
    return max_switch_id_;
  };

  
  bool flexfly_topology_simplified::node_to_netlink(node_id nid, node_id& net_id, int& offset) const {
    //std::cout << "27 " << std::endl;
    return true;
  };

  switch_id flexfly_topology_simplified::node_to_switch(node_id nid) const {
    //std::cout << "28 " << std::endl;
    switch_id swid = nid / (nodes_per_switch_);
    return swid;
  };

  /**
   * given two group id's, group1 and group2, returns the number of intergroup links 
   * 
   */
  int flexfly_topology_simplified::num_links_between_groups(int group1, int group2) const {
    //std::cout << "29 " << std::endl;
    return group_connectivity_matrix_[group1][group2];
  };

  /**
   * Checks that each in and out port of a given index of all optical switches
   * get connected to the exact same switches.
   **/
  void flexfly_topology_simplified::check_intergroup_connection() const {
    //std::cout << "30 " << std::endl;
    for (int i = 0; i < num_optical_switches_; i++) {
      int index = i + num_groups_ * switches_per_group_;
      std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_in_iter = switch_inport_connection_map_.find(index);
      std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_out_iter = switch_outport_connection_map_.find(index);
      const std::vector<switch_link*>& conn_in_vector = tmp_in_iter->second;
      const std::vector<switch_link*>& conn_out_vector = tmp_out_iter->second;
      for (int i = 0; i < conn_in_vector.size(); i++) {
        if (conn_in_vector[i]->src_sid != conn_out_vector[i]->dest_sid) {
          spkt_abort_printf("Optical Switch id: %d has different in-out connectivity\n",index);
        }
        switch_id swid = conn_in_vector[i]->src_sid;
        if (group_from_swid(swid) != i) {
          spkt_abort_printf("Optical Switch id: %d doesn't have properly-connected inout ports to group\n",index);  
        }
      }
    }
  };




  void flexfly_topology_simplified::check_routing_table() const {
    std::cout << "Routing table has " + std::to_string(routing_table_.size()) + " entries." << std::endl;
    assert((num_groups_  * switches_per_group_) == routing_table_.size());
    for (int i = 0; i < routing_table_.size(); i++) {
      for (int j = 0; j < routing_table_[i].size(); j++) {
        if (routing_table_[i][j])
          std::cout << " 1 ";
        else 
          std::cout << " 0 ";
      }
      std::cout << std::endl;
    }

    for (int i = 0; i < routing_table_.size(); i++) {
      for (int j = 0; j < routing_table_[i].size(); j++) {
        if (i == j) {
          continue;
        }

        assert(routing_table_[i][j]);
        std::cout << "Path from switch: " + std::to_string(i) + " to " << std::to_string(j) << std::endl;
        std::cout << "Path size: " + std::to_string(routing_table_[i][j]->path.size()) + " and path length: " + std::to_string(routing_table_[i][j]->path_length) << std::endl;
        for (int k = 0; k < routing_table_[i][j]->path_length; k++) {
          std::cout << "      switch id: " + std::to_string(routing_table_[i][j]->path[k]->switch_id) + " outport: "  + std::to_string(routing_table_[i][j]->path[k]->outport)<< std::endl;
        }
        assert(routing_table_[i][j]->path.size() == routing_table_[i][j]->path_length);
      }
      std::cout << std::endl;
    }
  };

  // NEWWWWWWW
  // Written 12/16/2017, used by optical_network to figure out which outports to direct a packet coming from group i to group j
  // From the optical network side, should use some decision algrotihm to at least figure out if the outgoing port can lead directly to
  // target switch
  void flexfly_topology_simplified::configure_optical_network(std::vector<std::vector<std::vector<int>>>& outport_options) const {
    //std::cout << "wsdihgairghierhi" << std::endl;
    auto connection_vector_iter = switch_outport_connection_map_.find(num_groups_ * switches_per_group_);
    assert(connection_vector_iter != switch_outport_connection_map_.end());
    auto connection_vector = connection_vector_iter->second;
    std::vector<std::vector<int>> target_group_set(num_groups_);
    uint32_t seed = 12;
    std::srand(seed);
    //std::cout << "31 " << std::endl;
    for (auto sl : connection_vector) {
      int dst_group = group_from_swid(sl->dest_sid);
      target_group_set[dst_group].push_back(sl->src_outport);
    }
    
    for (int i = 0; i < num_groups_; i++) {
      for (int j = 0; j < num_groups_; j++) {
        int count = 0;
        while (group_connectivity_matrix_[i][j] > count) {
          int ssize = target_group_set[j].size();
          int rand_num = rand() % ssize;
          outport_options[i][j].push_back(target_group_set[j][rand_num]);
          target_group_set[j].erase(target_group_set[j].begin() + rand_num);
          count++; 
        }
      }
    }
  };

  // FIGURE OUT HOW TO HAVE A GLOBAL ROUTING FOR ALL NODES. FIRST, TAKE A LOOK AT GROUP TO GROUP CONNECTIVITY MATRIX
  // THEN BASED ON THAT RANDOMLY CHOOSE INTERGROUP SWITCH PAIRS TO CONNECT, ASSUMING THAT THE GROUP CONNECTIVITY MATRIX
  // DOES NOT VIOLATE THE CONSTRAINTS
  // THEN PERFORM DIJKSTRA
  void flexfly_topology_simplified::route_single_switch_minimal(switch_id src, std::vector<std::vector<switch_id>>& adjacency_list) {
    //std::cout << "cibaikia" << std::endl;
    if (src > max_switch_id_)
      spkt_abort_printf("The switch for which we want to route should be within the topology");
    int distance_vector[max_switch_id_ + 1];
    bool visited_switches[max_switch_id_ + 1];
    switch_id parent_switch[max_switch_id_ + 1];

    // this is to make sure routing is equalized for all switches get used evenly for routing
    // and no 1 switch gets too much favoritism
    uint8_t used_for_routing[max_switch_id_ + 1];  

    // Initialization phase BEGIN
    for (int id = 0; id <= max_switch_id_; id++) {
      distance_vector[id] = INT_MAX;
      visited_switches[id] = false;
    }
    distance_vector[src] = 0;
    // Initialization phase END
    
    std::queue<switch_id> queue;
    queue.push(src);

    // Breadth first search portion
    while (!queue.empty()) {
      switch_id curr_switch = queue.front();
      queue.pop();
      visited_switches[curr_switch] = true;
      auto outgoing_edges = adjacency_list[curr_switch];
      // CONTINUE HERE
      for (auto neighbor_switch : outgoing_edges) {
        
        if ((distance_vector[neighbor_switch] == distance_vector[curr_switch] + 1) &&
            (curr_switch != src) &&
            used_for_routing[curr_switch] < used_for_routing[parent_switch[neighbor_switch]]) {
          
          parent_switch[neighbor_switch] = curr_switch;

        }

        if (distance_vector[neighbor_switch] > distance_vector[curr_switch] + 1) {
          distance_vector[neighbor_switch] = distance_vector[curr_switch] + 1;
          parent_switch[neighbor_switch] = curr_switch;
        }

        

        if (!visited_switches[neighbor_switch])
          queue.push(neighbor_switch);
      }
    }
    
    // Now that Dijkstra's is done, put routing information
    for (switch_id id = 0; id <= max_switch_id_ - 1; id++) {
      if (id == src) continue;
      
      if (routing_table_[src][id] != nullptr) delete routing_table_[src][id];
      routing_table_[src][id] = new flexfly_path();
      switch_id parent = parent_switch[id];
      switch_id current = id;
      
      while (current != src)  {
        if (group_from_swid(parent) != group_from_swid(current)) {
          switch_port_pair* spp1 = new switch_port_pair();
          switch_port_pair* spp2 = new switch_port_pair();
          spp1->switch_id = max_switch_id_;
          spp1->outport = current; // NOTE that for the optical network, the port number leading to a switch is the same as the switch_id
          spp2->switch_id = parent;
          spp2->outport = switches_per_group_ - 1 ;
          (routing_table_[src][id]->path).insert(routing_table_[src][id]->path.begin(), spp1);
          (routing_table_[src][id]->path).insert(routing_table_[src][id]->path.begin(), spp2);
          current = parent;
          parent = parent_switch[parent];
        } else {
          switch_port_pair* spp = new switch_port_pair();
          spp->switch_id = parent;
          auto it = switch_inport_connection_map_.find(current);
          //int port = 0;
          for (auto incoming_switch : it->second) {
            if (incoming_switch->src_sid == parent) {
              spp->outport = incoming_switch->src_outport;
              break;
            }
            //port++;
          }
          (routing_table_[src][id]->path).insert(routing_table_[src][id]->path.begin(), spp);
          current = parent;
          parent = parent_switch[parent];
        }
        
      }
    }
  };
/*
  void flexfly_topology_simplified::route_single_switch_minimal(switch_id src, std::vector<std::vector<switch_id>>& adjacency_list) {
    if (src > max_switch_id_)
      spkt_abort_printf("The switch for which we want to route should be within the topology");
    int distance_vector[max_switch_id_ + 1];
    bool visited_switches[max_switch_id_ + 1];
    std::vector<switch_id>& parent_switch = routing_table_2_[src];

    // this is to make sure routing is equalized for all switches get used evenly for routing
    // and no 1 switch gets too much favoritism
    uint8_t used_for_routing[max_switch_id_ + 1];  

    // Initialization phase BEGIN
    for (int id = 0; id <= max_switch_id_; id++) {
      distance_vector[id] = INT_MAX;
      visited_switches[id] = false;
    }
    distance_vector[src] = 0;
    // Initialization phase END
    
    std::queue<switch_id> queue;
    queue.push(src);

    // Breadth first search portion
    while (!queue.empty()) {
      switch_id curr_switch = queue.front();
      queue.pop();
      visited_switches[curr_switch] = true;
      auto outgoing_edges = adjacency_list[curr_switch];
      // CONTINUE HERE
      for (auto neighbor_switch : outgoing_edges) {
        
        if ((distance_vector[neighbor_switch] == distance_vector[curr_switch] + 1) &&
            (curr_switch != src) &&
            used_for_routing[curr_switch] < used_for_routing[parent_switch[neighbor_switch]]) {
          
          parent_switch[neighbor_switch] = curr_switch;

        }

        if (distance_vector[neighbor_switch] > distance_vector[curr_switch] + 1) {
          distance_vector[neighbor_switch] = distance_vector[curr_switch] + 1;
          parent_switch[neighbor_switch] = curr_switch;
        }

        if (!visited_switches[neighbor_switch])
          queue.push(neighbor_switch);
      }
    }
    
    // Now that Dijkstra's is done, put routing information

  };
  */
  void flexfly_topology_simplified::route_topology_minimal() {
    // Massage the virtual adjacency list of the graph 
    //std::cout << "siaoe" << std::endl;
    std::vector<std::vector<switch_id>> adjacency_list; // NOTE THAT FLEXFLY IS FUNDAMENTALLY A DIRECTED TOPOLOGY
                                                        // So the adjacency represents a directed topology
    adjacency_list.resize(max_switch_id_);
    for (auto row : routing_table_) {
      for (auto f_path : row) {
        if (nullptr != f_path) delete f_path;
      }
    }
    form_virtual_intergroup_topology(adjacency_list);
    for (int swid = 0; swid < max_switch_id_; swid++) {
      route_single_switch_minimal(swid, adjacency_list);
    }
    return;
  };

  /**
   * This function forms the virtual inter-group switch-pair connection
   * so that route_topology_minimal can actually use this to do some routing
   **/
  void flexfly_topology_simplified::form_virtual_intergroup_topology(std::vector<std::vector<switch_id>>& adjacency_list) {
    //std::cout << "tamade" << std::endl;
    std::vector<std::vector<switch_id>> groups_available_switches_outports;
    std::vector<std::vector<switch_id>> groups_available_switches_inports;
    groups_available_switches_outports.resize(num_groups_);
    groups_available_switches_inports.resize(num_groups_);
    std::srand(100);
    // 1) first we want to form the sets of switches that each group can use to
    //    to connect
    for (int k = 0; k < num_groups_; k++) {
      int id_offset = k * switches_per_group_;
      for (int m = 0; m < switches_per_group_; m++) {
        groups_available_switches_outports[k].push_back(id_offset + m);
        groups_available_switches_inports[k].push_back(id_offset + m);
      }
      for (int i = 0; i < switches_per_group_ - 1; i++) {
        switch_id sid_src = k * switches_per_group_ + i;
        for (int j = i + 1; j < switches_per_group_; j++) {
          switch_id sid_dst = k * switches_per_group_ + j;
          adjacency_list[sid_src].push_back(sid_dst);
          adjacency_list[sid_dst].push_back(sid_src);
        }
      }
    }

    for (int i = 0; i < num_groups_; i++) {
      for (int j = 0; j < num_groups_; j++) {
        int formed_pairs = 0;
        while (formed_pairs < group_connectivity_matrix_[i][j]) {
          switch_id src_id = groups_available_switches_outports[i][rand() % groups_available_switches_outports[i].size()];
          switch_id dst_id = groups_available_switches_inports[j][rand() % groups_available_switches_inports[j].size()];
          adjacency_list[src_id].push_back(dst_id);
          // Delete the switches from the sets of available switches from either groups
          auto src_it = std::find(groups_available_switches_outports[i].begin(), 
                                  groups_available_switches_outports[i].end(), 
                                  src_id);
          auto dst_it = std::find(groups_available_switches_inports[j].begin(), 
                                  groups_available_switches_inports[j].end(), 
                                  dst_id);
          assert(src_it != groups_available_switches_outports[i].end());
          assert(dst_it != groups_available_switches_inports[j].end());
          groups_available_switches_outports[i].erase(src_it);
          groups_available_switches_inports[j].erase(dst_it);
          formed_pairs++;
        }
      }
    }
  };

  /**
   * This 
   **/
  void flexfly_topology_simplified::minimal_route_to_switch_optical(switch_id src, switch_id dst, int& outport_arg) const {
    //std::cout << "weeeeeeeee" << std::endl;
    flexfly_path* fpath = routing_table_[src][dst];
    assert(fpath != nullptr);
    for (auto spp : fpath->path) {
      if (spp->switch_id == max_switch_id_) {
        outport_arg = spp->outport; 
      }
    }
  };

  bool flexfly_topology_simplified::minimal_route_special_flexfly(switch_id src, switch_id dst, switch_id curr_switch, int& outport) const {
    //std::cout << "yuuuuuuuuu" << std::endl;
    flexfly_path* fp = routing_table_[src][dst];
    for (auto spp : fp->path) {
      if (curr_switch == spp->switch_id) {
        outport = spp->outport;
        return true;
      }
    }
    return false;
  };
/*
  bool flexfly_topology_simplified::minimal_route_special_flexfly_2(switch_id src, switch_id dst, switch_id curr_switch, int& outport) const {
    switch_id curr = dst;
    switch_id parent = routing_table_2_[src][dst];
    while (parent != curr_switch) {
      curr = routing_table
    }
    for (auto spp : fp->path) {
      if (curr_switch == spp->switch_id) {
        outport = spp->outport;
        return true;
      }
    }
    return false;
  };
  */
}
}

