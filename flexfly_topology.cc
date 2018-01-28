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
#include "flexfly_topology.h"
#include "connectivity_matrix.h"
#include "flexfly_routing_algo.h"

namespace sstmac {
namespace hw {

 
 flexfly_topology::flexfly_topology(sprockit::sim_parameters* params) : 
                              structured_topology(params,InitMaxPortsIntra::I_Remembered, 
                                                  InitGeomEjectID::I_Remembered) {
 	num_groups_ = params->get_optional_int_param("groups", 7); // controls g
 	switches_per_group_ = params->get_optional_int_param("switches_per_group", 6); // controls a
  switches_per_group_ = num_groups_ - 1;
 	nodes_per_switch_ = params->get_optional_int_param("nodes_per_switch", 4);
  optical_switch_radix_ = params->get_optional_int_param("optical_switch_radix", switches_per_group_);
  is_simplified_model_ = false;
 	num_optical_switches_ = switches_per_group_;
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
  


  distance_matrix_.resize(num_groups_ * switches_per_group_);
  for (int i = 0; i < num_groups_ * switches_per_group_; i++) {
    distance_matrix_[i].resize(num_groups_ * switches_per_group_);
  }
 	setup_flexfly_topology();
  
  check_intergroup_connection();
  configure_default_simpler_model(num_groups_, optical_inout_connectivity_); 
  
  routing_table_.resize(num_groups_ * switches_per_group_);
  for (int i = 0; i < num_groups_ * switches_per_group_; i++) {
    routing_table_[i].resize(num_groups_ * switches_per_group_);
  }
  for (int i = 0; i < num_groups_ * switches_per_group_; i++) {
    dijkstra_minimal_route(i, 
                            switch_outport_connection_map_, 
                            switch_inport_connection_map_, 
                            distance_matrix_[i],
                            routing_table_[i],
                            optical_inout_connectivity_,
                            num_groups_,
                            switches_per_group_);
  }
  std::cout << "Routing succesful" <<std::endl;
  check_routing_table();
  
  updated_routing_table_ = false; // turn this to false first.

  /****************************************************
   ************ Testing for configure_optical ********* 
   ****************************************************/
  
  /*
  for (int optical_index = num_groups_ * switches_per_group_; 
        optical_index < num_optical_switches_ + num_groups_ * switches_per_group_; 
        optical_index++) {
    std::pair<int, std::vector<int>> curr_pair (optical_index, std::vector<int> (optical_switch_radix_));
    optical_inout_connectivity_.insert(curr_pair);
  }
  
  configure_optical_switches_general(group_connectivity_matrix_, optical_inout_connectivity_);
  */
  /****************************************************
   ************ Testing for configure_optical ********* 
   ****************************************************/

 }

 flexfly_topology::~flexfly_topology() {
   for (const std::pair<switch_id, std::vector<switch_link*>> elem : switch_outport_connection_map_) {
     const std::vector<switch_link*>& conn_vector = elem.second;  
     for (auto const&  switch_link_ptr : conn_vector) {
        if (switch_link_ptr)
          free(switch_link_ptr);
     }
   }

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
 void flexfly_topology::setup_flexfly_topology() {
   // first step: connect all the switches within the same group together
   for (int group = 0; group < num_groups_; group++) {
     switch_id group_offset = group * switches_per_group_;
     for (int index = 0; index < switches_per_group_; index++) {
       switch_id swid = group_offset + index;
       for (int target_index = index + 1; target_index < switches_per_group_; target_index++) {
         switch_id target_swid = group_offset + target_index;
         connect_switches(swid, target_swid, Link_Type::electrical);
         connect_switches(target_swid, swid, Link_Type::electrical);
       }
     }  
   }

   // second step: connect all the switches within groups to the optical switches
   for (int group = 0; group < num_groups_; group++) {
     switch_id group_offset = group * switches_per_group_;
     for (int index = 0; index < switches_per_group_; index++) {
       switch_id swid = group_offset + index;
       switch_id optical_swid = num_groups_ * switches_per_group_ + index;
       connect_switches(swid, optical_swid, Link_Type::optical);
       connect_switches(optical_swid, swid, Link_Type::optical);
     }
   }
 }

 /**
  * Called only once at the start of the constructor. Finds all the paths 
  * from 1 switch to the other given the initial optical switch inout port
  * connectivity configurations. Uses minimal shortest path routing for this.
  **/
 void flexfly_topology::setup_routing_table() {
  int total_switches = num_optical_switches_ + num_groups_ * switches_per_group_;
  if (routing_table_.size() == 0) 
    routing_table_.resize(total_switches);
  for (int i = 0; i < total_switches; i++) {
    routing_table_.resize(total_switches);
  }
 }

 void flexfly_topology::configure_metis(metis_config* configuration) const {
	if (!configuration) {
		return;
	}
	configuration->nvtxs = (idx_t) num_groups_ * switches_per_group_;
 }

 /**
  * IMPORTANT: This function will route the minimal path
  **/
 void flexfly_topology::minimal_route_to_switch(switch_id src_switch_addr, 
 												switch_id dst_switch_addr, 
 												routable::path& path) const {
  int src_group = src_switch_addr/switches_per_group_;
  int dst_group = src_switch_addr/switches_per_group_;
  // both the source and dest switch cannot be the same, because the caller 
  // should not call this function if this were the case
  assert(src_switch_addr != dst_switch_addr);
  // if both nodes are connected to the same switch
  if (src_group == dst_group) {
    std::vector<topology::connection> conns;  
    connected_outports(src_switch_addr, conns);
    for (int i = 0; i < conns.size(); i++) {
      if (conns[i].dst == dst_switch_addr) {
        path.set_outport(i);
        return;
      }
    }
  } else {
  }

 };

 int flexfly_topology::get_output_port(int src_switch, int dst_switch) const {
  auto sl_iter = switch_outport_connection_map_.find(src_switch);
  assert(sl_iter != switch_outport_connection_map_.end());
  auto slv = sl_iter->second;
  int i = 0;
  for (auto sl : slv) {
    if (sl->dest_sid == dst_switch)
      return i;
    i++; 
  }
  return -1;
 };

 /**
  * checks if a given switch id in the topology is an optical switch or not
  */
 bool flexfly_topology::is_optical_switch(switch_id swid) const {
  if (!valid_switch_id(swid)) {
 		return false;
 	} 
 	return swid >= num_groups_ * switches_per_group_;
 };

 /**
  * checks if a given switch id in the topology is an electrical switch or not
  */
 bool flexfly_topology::is_electrical_switch(switch_id swid) const {
  return valid_switch_id(swid) && !is_optical_switch(swid);
 };

 // DONE (RECHECK)
 void flexfly_topology::configure_optical_or_electrical_port_params(switch_id swid, 
                                                                    std::string& str, 
                                                                    sprockit::sim_parameters* switch_params) const {
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
 };

 // DONE
 void flexfly_topology::configure_individual_port_params(switch_id src,
          												 sprockit::sim_parameters* switch_params) const {
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
  * NOTE: This member function should form a bidirectional switch_link
  */
  void flexfly_topology::connect_switches(switch_id src, switch_id dst, Link_Type ltype) {
    std::vector<switch_link*>& src_outport_connection_vector = switch_outport_connection_map_[src];
    std::vector<switch_link*>& dst_inport_connection_vector = switch_inport_connection_map_[dst];
    switch_link *conns = new switch_link();
    conns->src_sid = src;
    conns->dest_sid = dst;
    conns->dest_inport = dst_inport_connection_vector.size(); 
    conns->src_outport = src_outport_connection_vector.size();
    conns->type = ltype;
    int src_port = src_outport_connection_vector.size();
    src_outport_connection_vector.push_back(conns);
    dst_inport_connection_vector.push_back(conns);
    return;
 }


void flexfly_topology::connected_outports(const switch_id src, 
                                            std::vector<topology::connection>& conns) const {
  std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator got = switch_outport_connection_map_.find(src);
  int cidx = 0;
  if (got != switch_outport_connection_map_.end()) {
    const std::vector<switch_link*>& switch_link_vectors = got->second;
    for (switch_link* current_switch_link : switch_link_vectors) {
      conns.push_back(topology::connection());
      conns[cidx].src = src;
      conns[cidx].dst = current_switch_link->dest_sid;
      conns[cidx].src_outport = cidx; 
      conns[cidx].dst_inport = current_switch_link->dest_inport;
      conns[cidx].link_type = current_switch_link->type;
      cidx++;
    }
  }
}


bool flexfly_topology::switch_id_slot_filled(switch_id sid) const {
  return (sid < max_switch_id_);
}

 void flexfly_topology::configure_vc_routing(std::map<routing::algorithm_t, int>& m) const {
  m.insert({routing::minimal, 3});
  m.insert({routing::minimal_adaptive, 3});
  m.insert({routing::valiant, 3});
  m.insert({routing::ugal, 3});
  return;
 };

  switch_id flexfly_topology::node_to_ejection_switch(node_id addr, uint16_t& port) const {
    switch_id swid = addr / nodes_per_switch_; // this gives us the switch id of the switch node addr is connected to
    std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_iter = switch_outport_connection_map_.find(swid);
    const std::vector<switch_link*>& conn_vector = tmp_iter->second;
    port = std::max((int) (conn_vector.size() - 1), 0) + ((int) swid) * nodes_per_switch_; // CHECK THIS AGAIN
    return swid;
  };
  
  switch_id flexfly_topology::node_to_injection_switch(node_id addr, uint16_t& port) const {
    return node_to_ejection_switch(addr, port);
  };

  /**
   * NOTE: This method does not include the hop to an optical switch
   **/
  int flexfly_topology::minimal_distance(switch_id src, switch_id dst) const {
    //std::cout << "src switch id: " << std::to_string(src) << " dst switch id: " << std::to_string(dst) << std::endl;
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
        if (tmp_link->type == electrical)
          continue;
      }
      return two_or_three ? 2 : 3;
    }
  };

  int flexfly_topology::num_hops_to_node(node_id src, node_id dst) const {
    int src_swid = src / (nodes_per_switch_);
    int dst_swid = dst / (nodes_per_switch_);
    int min_dist = distance_matrix_[src_swid][dst_swid];
    return min_dist + 2; // added by 2 because each node is 1 hop away from it's switch
  };

  void flexfly_topology::nodes_connected_to_injection_switch(switch_id swid, 
                                                              std::vector<injection_port>& nodes) const {
    int i = 0;
    int port_offset = switches_per_group_;
    for (int i = 0; i < nodes_per_switch_; i++) {
      int node_id = swid * nodes_per_switch_ + i;
      int port_ind = port_offset + i;
      injection_port ijp;
      ijp.nid = node_id;
      ijp.port = port_ind;
      nodes.push_back(ijp);
    }

  };

  void flexfly_topology::nodes_connected_to_ejection_switch(switch_id swid, 
                                                              std::vector<injection_port>& nodes) const { 
    flexfly_topology::nodes_connected_to_injection_switch(swid, nodes);
  };

  // returns the group id of a given switch
  int flexfly_topology::group_from_swid(switch_id swid) const {
    return swid / (switches_per_group_);
  };

  /**
   * Used to print the connectivity matrix of the entire topology 
   * Doesn't need to take in any sort of arguments and does not return any
   * argument
   */
  void flexfly_topology::print_topology() const {
    for (std::pair<switch_id, std::vector<switch_link*>> key_val_pair : 
                switch_outport_connection_map_ ) {
      switch_id swid = key_val_pair.first;
      print_port_connection_for_switch(swid);
    }
  };

  /**
   * prints out the all the connections for each switch
   */
  void flexfly_topology::print_port_connection_for_switch(switch_id swid) const {
    std::unordered_map<switch_id, std::vector<switch_link*>>::const_iterator tmp_iter = 
                                        switch_outport_connection_map_.find(swid);
    if (tmp_iter == switch_outport_connection_map_.end()) {
      return;
    }

    const std::vector<switch_link*>& connection_vector = tmp_iter->second;

    std::stringstream message;
    int i = 0;
    for (switch_link* sl_ptr : connection_vector) {
      // check if null, if null have to throw an error 
      if (sl_ptr) {
        message << "Dest switch_id: " << std::to_string(sl_ptr->dest_sid);
        message << " Dest inport: " << std::to_string(sl_ptr->dest_inport);
        message << " Link type: ";
        if (sl_ptr->type == Link_Type::electrical) {
          message << "ELECTRICAL" << std::endl;
        } else {
          message << "OPTICAL" << std::endl;
        }
      } else {
        spkt_abort_printf("A switch link with swid: %d is null\n", swid);
      }
      i++;
    }
    std::cout << message.str(); 
  };


  /**
   * @brief num_endpoints To be distinguished slightly from nodes.
   * Multiple nodes can be grouped together with a netlink.  The netlink
   * is then the network endpoint that injects to the switch topology
   * @return
   */
  int flexfly_topology::num_netlinks() const {
    std::cout << "num_netlinks?" << std::endl;
    return 1;
  }; 

  switch_id flexfly_topology::max_netlink_id() const {
    std::cout << "max_netlink_id?" << std::endl;
    return max_switch_id_;
  };

  bool flexfly_topology::netlink_id_slot_filled(node_id nid) const {
    std::cout << "netlink_id_slot_filled?" << std::endl;
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
  switch_id flexfly_topology::netlink_to_injection_switch(
        netlink_id nodeaddr, uint16_t& switch_port) const {
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
  switch_id flexfly_topology::netlink_to_ejection_switch(
        netlink_id nodeaddr, uint16_t& switch_port) const {
    return max_switch_id_;
  };

  
  bool flexfly_topology::node_to_netlink(node_id nid, node_id& net_id, int& offset) const {
    return true;
  };

  switch_id flexfly_topology::node_to_switch(node_id nid) const {
    switch_id swid = nid / (nodes_per_switch_);
    return swid;
  };

  /**
   * given two group id's, group1 and group2, returns the number of intergroup links 
   * 
   */
  int flexfly_topology::num_links_between_groups(int group1, int group2) const {
    return group_connectivity_matrix_[group1][group2];
  };

  /**
   * Checks that each in and out port of a given index of all optical switches
   * get connected to the exact same switches.
   **/
  void flexfly_topology::check_intergroup_connection() const {
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

  /**
   * This is the key function that will be called by the electrical switches
   **/
  void flexfly_topology::route_minimal(int src_switch, int dst_switch, flexfly_packet* f_packet) {
    assert(src_switch != dst_switch);
    if (!updated_routing_table_) {
      // if haven't updated yet, just use the routing table
      flexfly_path *flexpath = routing_table_[src_switch][dst_switch];
      f_packet->set_path(flexpath);
      return;
    } else {
      // reroute the network here
    }
  };

  /**
   * Given an optical switch id, and the input-output connection vector, 
   * this function determines reupdates the inout connectivity vector for
   * this switch, with each index being the input port number.
   * NOTE: this function is only called by optical switches which will
   *       update the topology's inout connection table
   **/
  void flexfly_topology::optical_switch_update_inout(int optical_swid, std::vector<int>& inout_configuration) {
    assert(is_optical_switch(optical_swid));
    int port_count = inout_configuration.size();
    std::vector<int> &inout_vector = optical_inout_connectivity_[optical_swid];
    for (int i = 0; i < port_count; i++) {
      inout_vector[i] = inout_configuration[i];
    }
    updated_routing_table_ = true; // now update this to true
    return;
  };

  void flexfly_topology::check_routing_table() const {
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
     /*
    }
    for (const std::vector<flexfly_path*>& routing_entry_iter : routing_table_) {
      //std::cout << " got in here for: " + std::to_string(i) << std::endl;
      //const std::vector<flexfly_path *>& switch_path = routing_entry_iter.second; 
      std::string str = "Source Switch id: " + std::to_string(i) + " ";
      int j = 0;
      for (flexfly_path* f_path : routing_entry_iter) {
        if (i == j) continue;

        //const flexfly_path* f_path = flexfly_path_elem;
        str += ("Dest Switch id: " + std::to_string(j) + " \n");
        assert(f_path);
        for (int i = 0; i < f_path->path.size(); i++) {
          switch_port_pair* spp = f_path->path[i];
          str += ("Switch id: " + std::to_string(spp->switch_id) + " outport: " + std::to_string(spp->outport) + "\n"); 
        }
        str += "\n";
        j++;
        
      }
      std::cout << str << std::endl;
      i++;
    }
    std::cout << "Getting out of check_routing_table" << std::endl;
    */
  };

  /**
   * Note: total_switches contains the number of optical switches as well
   * Note: unordered_map will have to be initialized by the caller
   **/
  void flexfly_topology::configure_optical_switches_general(std::vector<std::vector<int>>& connectivity_requirement,
                                                            std::unordered_map<int, std::vector<int>>& optical_switch_inout_config) {
    /**********************************************************************************************************
    ******************************** WORKING SECTION **********************************************************
    **********************************************************************************************************/
    // vector or lists for a group, each element in list is an 
    // optical switch's id that this group is connected to
    std::cout << "In configure_optical_switches_general function" << std::endl;
    int total_switches = num_optical_switches_ + num_groups_ * switches_per_group_;
    
    // potential number of ways to connect 1 group from another
    std::vector<std::vector<int>> group_physical_connectivity(num_groups_); 
    
    /*
     * Initialize the group_connectivity_matrix to all zeros
     */
    for (int i = 0 ; i < num_groups_; i++) {
      group_physical_connectivity[i].resize(num_groups_);
      for (int j = 0; j < num_groups_; j++) {
        group_physical_connectivity[i][j] = 0;
      }
    }

    std::vector<std::vector<int>> group_nodes(num_groups_);
    std::unordered_map<int, std::vector<int>> optical_switch_nodes;
  

    /*
     * Use vectors/list to emulate a graph node approach to solving this problem
     * First initiate the flexfly topology as a 3 layer graph. For each group node,
     * figure out all the all the 
     */
    for (int group = 0; group < num_groups_; group++) {
      int group_offset = group * switches_per_group_;
      for (int switch_index = 0; switch_index < switches_per_group_; switch_index++) {
        int curr_switch = group_offset + switch_index;
        auto switch_link_vector_iter = switch_outport_connection_map_.find(curr_switch);
        assert(switch_link_vector_iter != switch_outport_connection_map_.end());
        const std::vector<switch_link *>& switch_link_vector = switch_link_vector_iter->second;
        for (auto sl : switch_link_vector) {
          if (sl->type == Link_Type::electrical) continue;
          // dest_sid has to belong to an optical switch if topology wiring is done correctly
          assert(is_optical_switch(sl->dest_sid)); 
          group_nodes[group].insert(group_nodes[group].end(),  sl->dest_sid);
        }
      }
    }
    /*
     * Now do the same for the second layer of optical switches
     */
    int optical_switch_offset = num_groups_ * switches_per_group_;
    for (int opt_switch = 0; opt_switch < num_optical_switches_; opt_switch++) {
      int curr_switch = optical_switch_offset + opt_switch; // CURRENT OPTICAL SWITCH
      auto switch_link_vector_iter = switch_outport_connection_map_.find(curr_switch);
      const std::vector<switch_link *>& switch_link_vector = switch_link_vector_iter->second;
      optical_switch_nodes.insert(std::pair<int, std::vector<int>> (curr_switch, std::vector<int>()));
      for (auto sl : switch_link_vector) {
        optical_switch_nodes[curr_switch].insert(optical_switch_nodes[curr_switch].end(), group_from_swid(sl->dest_sid));
      }
    }

    

    for (int k = 0; k < total_switches; k++) {
      if (is_optical_switch(k)) continue;
      auto connection_vector_iter = switch_outport_connection_map_.find(k);
      const std::vector<switch_link*>& switch_outport_vector = connection_vector_iter->second;
      int curr_switch_group = group_from_swid(k);
      for (auto sl : switch_outport_vector) {
        if (!is_optical_switch(sl->dest_sid)) continue;
        auto optical_connection_vector_iter = switch_outport_connection_map_.find(sl->dest_sid);
        const std::vector<switch_link*>& optical_connection_vector = optical_connection_vector_iter->second;
        for (auto op_sl : optical_connection_vector) {
          int dest_switch_group = group_from_swid(op_sl->dest_sid);
          group_physical_connectivity[curr_switch_group][dest_switch_group]++;
        }
      }
    }
    /**
     * FOR DEBUGGING PURPOSE
     **/
    std::cout << "Connectivity_requirement" << std::endl;
    for (auto x : connectivity_requirement) {
      for (auto y : x) {
        std::cout << std::to_string(y) + " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    /**
     * FOR DEBUGGING PURPOSE
     **/
    /**
     * FOR DEBUGGING PURPOSE
     **/
    for (auto x : group_physical_connectivity) {
      for (auto y : x) {
        std::cout << std::to_string(y) + " ";
      }
      std::cout << std::endl;
    }
    /**
     * FOR DEBUGGING PURPOSE
     **/

    int hack = 6;
    for (int i = 0; i < num_optical_switches_; i++) {
      std::cout << "Optical Switch: " + std::to_string(i) << std::endl;
      //std::cout << std::to_string(optical_switch_nodes[].size())
      for (auto dummy_target_group : optical_switch_nodes[hack]) {
        std::cout << " " + std::to_string(dummy_target_group) + " ";
      }
      hack++;
      std::cout << std::endl;
    }
    /*
     * Main body of the algorithm
     */
    for (int i = 0; i < num_groups_; i++) {
      for (int j = 0; j < num_groups_; j++) {
        while (connectivity_requirement[i][j] > 0) {
          bool found_switch = false;
          int group_node_port = 0;
          //int used_optical_switch_id = -1;
          for (int optical_switch : group_nodes[i]) {
            
            
            if (optical_switch < 0) {
              group_node_port++;
              continue;
            }
            //used_optical_switch_id = optical_switch;
            std::cout << "here once " + std::to_string(optical_switch) + " for pair " + std::to_string(i) + " to "+ std::to_string(j) << std::endl;
            int optical_node_port = 0;
            for (auto target_group : optical_switch_nodes[optical_switch]) {
              if (target_group < 0 || target_group != j) {
                optical_node_port++;
                continue;
              }
              auto switch_inport_vector_iter = switch_inport_connection_map_.find(optical_switch);
              const std::vector<switch_link *>& switch_inport_vector = switch_inport_vector_iter->second;
              bool pass_check = true;
              assert(target_group == j);
              for (auto sl : switch_inport_vector) {
                int src_group = group_from_swid(sl->src_sid);
                //if (src_group == i || ) continue;
                if (src_group == i || std::find(group_nodes[src_group].begin(), group_nodes[src_group].end(), optical_switch) == group_nodes[src_group].end()) continue;
                if (group_physical_connectivity[src_group][j] - 1 < connectivity_requirement[src_group][j]) {
                  pass_check = false;
                  break;
                } 
              } 
              // assume only each group's connection will go to the same switch ONCE
              if (pass_check) {
                for (auto outgoing_group : optical_switch_nodes[optical_switch]) {
                  if (outgoing_group < 0 || outgoing_group == j) continue;
                  
                  if (group_physical_connectivity[i][outgoing_group] - 1 < connectivity_requirement[i][outgoing_group]) {
                    std::cout << std::to_string(outgoing_group) << std::endl;
                    pass_check = false;
                    break;
                  }
                }
              }

              if (pass_check) {
                found_switch = true;
                int optical_switch_inport = -1;
                for (auto sl : switch_inport_vector) { 
                  int src_group = group_from_swid(sl->src_sid);
                  for (auto y : group_nodes[src_group]) {
                    //if (y < 0) continue;
                    if (y == optical_switch)
                      group_physical_connectivity[src_group][j]--;
                  }
                }
                for (auto x : optical_switch_nodes[optical_switch]) {
                  if (x < 0 || x == j) continue; // prevents the double counting
                  group_physical_connectivity[i][x]--;
                }
                auto switch_inout_connection_iter = optical_switch_inout_config.find(optical_switch);
                assert(switch_inout_connection_iter != optical_switch_inout_config.end());
                std::vector<int>& switch_inout_config = switch_inout_connection_iter->second;
                switch_inout_config[optical_switch_inport] = optical_node_port;
    
                connectivity_requirement[i][j]--;
                assert(optical_switch_nodes[optical_switch][optical_node_port] == j);
                optical_switch_nodes[optical_switch][optical_node_port] = -1; 
                assert(group_nodes[i][group_node_port] == optical_switch);
                group_nodes[i][group_node_port] = -1;
                break;
              }
            }
            group_node_port++;
            if (found_switch) {
              std::cout << "Connectivity_requirement" << std::endl;
              for (auto x : connectivity_requirement) {
                for (auto y : x) {
                  std::cout << std::to_string(y) + " ";
                }
                std::cout << std::endl;
              }
              std::cout << std::endl;

              std::cout << "current form of group_physical_connectivity" << std::endl;
              for (auto x : group_physical_connectivity) {
                for (auto y : x) {
                  std::cout << std::to_string(y) + " ";
                }
              std::cout << std::endl;
              }
              break;
            }
            
          }

          if (!found_switch) {
            std::cout << "Failed at " << std::to_string(i) <<  ", " << std::to_string(j) << std::endl;
            spkt_abort_printf("Cannot possibly find the proper configuration for optical switches")
          } 
          //connectivity_requirement[i][j]--;
        }
      }
    }
  /**********************************************************************************************************
  ******************************** WORKING SECTION **********************************************************
  **********************************************************************************************************/
    return;
  };


/**
   * Note: total_switches contains the number of optical switches as well
   * Note: unordered_map will have to be initialized by the caller
   **/
  void flexfly_topology::configure_optical_switches_general3(std::vector<std::vector<int>>& connectivity_requirement,
                                                            std::unordered_map<int, std::vector<int>>& optical_switch_inout_config) {
    /**********************************************************************************************************
    ******************************** WORKING SECTION **********************************************************
    **********************************************************************************************************/
    // vector or lists for a group, each element in list is an 
    // optical switch's id that this group is connected to
    int total_switches = num_optical_switches_ + num_groups_ * switches_per_group_;
    
    // entry i,j is the set of all optical switches that can still connect group i to j
    std::vector<std::vector<std::vector<int>>> group_physical_connectivity(num_groups_); 
    
    /*
     * Initialize the group_connectivity_matrix to all zeros
     */
    for (int i = 0 ; i < num_groups_; i++) {
      group_physical_connectivity[i].resize(num_groups_);
    }

    std::vector<std::vector<int>> group_nodes(num_groups_);
    std::unordered_map<int, std::vector<int>> optical_switch_nodes;  
    /*
     * Use vectors/list to emulate a graph node approach to solving this problem
     * First initiate the flexfly topology as a 3 layer graph. For each group node,
     * figure out all the all the 
     */
    for (int group = 0; group < num_groups_; group++) {
      int group_offset = group * switches_per_group_;
      for (int switch_index = 0; switch_index < switches_per_group_; switch_index++) {
        int curr_switch = group_offset + switch_index;
        auto switch_link_vector_iter = switch_outport_connection_map_.find(curr_switch);
        assert(switch_link_vector_iter != switch_outport_connection_map_.end());
        const std::vector<switch_link *>& switch_link_vector = switch_link_vector_iter->second;
        for (auto sl : switch_link_vector) {
          if (sl->type == Link_Type::electrical) continue;
          // dest_sid has to belong to an optical switch if topology wiring is done correctly
          assert(is_optical_switch(sl->dest_sid)); 
          group_nodes[group].insert(group_nodes[group].end(),  sl->dest_sid);
        }
      }
    }
    /*
     * Now do the same for the second layer of optical switches
     */
    int optical_switch_offset = num_groups_ * switches_per_group_;
    for (int opt_switch = 0; opt_switch < num_optical_switches_; opt_switch++) {
      int curr_switch = optical_switch_offset + opt_switch; // CURRENT OPTICAL SWITCH
      auto switch_link_vector_iter = switch_outport_connection_map_.find(curr_switch);
      const std::vector<switch_link *>& switch_link_vector = switch_link_vector_iter->second;
      optical_switch_nodes.insert(std::pair<int, std::vector<int>> (curr_switch, std::vector<int>()));
      for (auto sl : switch_link_vector) {
        optical_switch_nodes[curr_switch].push_back(group_from_swid(sl->dest_sid));
      }
    }

    /*
     * Fill up the group_physical_connectivity matrix of sets
     */
    for (int switch_id = switches_per_group_ * num_groups_; switch_id < total_switches; switch_id++) {
      auto outport_vector_iter = switch_outport_connection_map_.find(switch_id);
      auto inport_vector_iter = switch_inport_connection_map_.find(switch_id);
      assert(outport_vector_iter != switch_outport_connection_map_.end());
      auto outport_vector = outport_vector_iter->second;
      auto inport_vector = inport_vector_iter->second;
      for (auto sl1 : outport_vector) {
        assert(!is_optical_switch(sl1->dest_sid));
        int j = group_from_swid(sl1->dest_sid);
        for (auto sl2 : inport_vector) {
          int i = group_from_swid(sl2->src_sid);
          if (std::find(group_physical_connectivity[i][j].begin(), 
                        group_physical_connectivity[i][j].end(), 
                        switch_id) == group_physical_connectivity[i][j].end()) {
            group_physical_connectivity[i][j].push_back(switch_id);
          }
        }
      }
    }

    std::cout << "Printing set" << std::endl;
    for (auto x : group_physical_connectivity[0][1]) {
      std::cout << std::to_string(x) + " ";
    }
    std::cout << std::endl;


    /*
     * Main body of the algorithm
     */
    for (int i = 0; i < num_groups_; i++) {
      for (int j = 0; j < num_groups_; j++) {
        std::cout << "i, j is: " + std::to_string(i) + " " + std::to_string(j) << std::endl;
        while (connectivity_requirement[i][j] > 0) {
          // Find the potential optical switch to use
          bool use_potential_optical_switch = false;

          std::cout << "Printing set" << std::endl;
          for (auto x : group_physical_connectivity[i][j]) {
            std::cout << std::to_string(x) + " ";
          }
          std::cout << std::endl;
          std::cout << "Printing set" << std::endl;
          for (auto x : group_physical_connectivity[2][4]) {
            std::cout << std::to_string(x) + " ";
          }
          std::cout << std::endl;

          for (auto potential_optical_switch : group_physical_connectivity[i][j]) {
            if (use_potential_optical_switch) break;
            
            std::cout << "potential_optical_switch is: " + std::to_string(potential_optical_switch) << std::endl;
            assert(std::find(optical_switch_nodes[potential_optical_switch].begin(), 
                            optical_switch_nodes[potential_optical_switch].end(), 
                            j) != optical_switch_nodes[potential_optical_switch].end());
            assert(std::find(group_nodes[i].begin(), 
                            group_nodes[i].end(), 
                            potential_optical_switch) != group_nodes[i].end());

            std::vector<std::pair<int, int>> entries_affected; // the (i,j) entries that are affected
            bool pass_check = true;
            for (auto sl : switch_inport_connection_map_[potential_optical_switch]) {
              int src_group = group_from_swid(sl->src_sid);
              // check if this potential optical_switch is still connected to the source group
              if (std::find(group_nodes[src_group].begin(), 
                            group_nodes[src_group].end(),
                            potential_optical_switch) != group_nodes[src_group].end() && src_group != i) {

                // means that the group still might have an option for using this switch
                assert(!is_optical_switch(sl->src_sid));
                assert(delete_item(group_physical_connectivity[src_group][j], potential_optical_switch));
                entries_affected.push_back(std::make_pair(src_group, j));
               // std::cout << "union set size: " + std::to_string(union_set_size(group_physical_connectivity[src_group], j, num_groups_)) << std::endl;
                if (union_set_size(group_physical_connectivity[src_group], 0, num_groups_) < remaining_row_requirement(connectivity_requirement[src_group])
                      || !check_remaining(group_physical_connectivity[src_group], connectivity_requirement[src_group])) {
                  std::cout << "cibai here la for dshGI;Herwgiuaerhugi src_group: " + std::to_string(src_group) << std::endl ;
                  pass_check = false;
                  //break;
                } 
              }
            }

            bool second_checkpoint = true;

            if (pass_check) {
              for (auto target_group : optical_switch_nodes[potential_optical_switch]) {
                if (target_group == j) continue;
                assert(delete_item(group_physical_connectivity[i][target_group], potential_optical_switch));
                entries_affected.push_back(std::make_pair(i, target_group));
                if (union_set_size(group_physical_connectivity[i], 0, num_groups_) < remaining_row_requirement(connectivity_requirement[i])
                    || !check_remaining(group_physical_connectivity[i], connectivity_requirement[i])) {
                  std::cout << "cibai here la for dst_group: " + std::to_string(target_group) << std::endl;
                  pass_check = false;
                  second_checkpoint = false;
                  //break;
                } 
              }
            } else {
              for (auto p : entries_affected) {
                int row = p.first;
                int col = p.second;
                add_item(group_physical_connectivity[row][col], potential_optical_switch);
              }
            }



            if (pass_check) {
              use_potential_optical_switch = true;
              assert(delete_item(group_physical_connectivity[i][j], potential_optical_switch));
              assert(delete_item(group_nodes[i], potential_optical_switch));
              assert(delete_item(optical_switch_nodes[potential_optical_switch], j));
              connectivity_requirement[i][j]--;
              break;
            } else if (!pass_check && !second_checkpoint) {
              for (auto p : entries_affected) {
                int row = p.first;
                int col = p.second;
                add_item(group_physical_connectivity[row][col], potential_optical_switch);
              }
            }
            
          }

          //assert(use_potential_optical_switch);
          if (!use_potential_optical_switch) {
            std::cout << "FAILED FOR i: " + std::to_string(i) + " j: " + std::to_string(j) << std::endl; 
            connectivity_requirement[i][j]--;
            //assert(use_potential_optical_switch); 
          }
        }
      }
    }
  /**********************************************************************************************************
  ******************************** WORKING SECTION **********************************************************
  **********************************************************************************************************/
    return;
  };

  bool flexfly_topology::delete_item(std::vector<int>& set, int optical_switch) {
    bool to_ret = false;
    int i;
    for (i = 0; i < set.size(); i++) {
      if (optical_switch == set[i]) {
        to_ret = true;
        set.erase(set.begin() + i);
        break;
      }
    }
    //if (to_ret) set.erase(set.begin() + i);
    return to_ret;
  };

  void flexfly_topology::add_item(std::vector<int>& set, int optical_switch) {
    if (std::find(set.begin(), set.end(), optical_switch) == set.end()) 
      set.insert(set.end(), optical_switch);
  };

  int flexfly_topology::remaining_row_requirement(std::vector<int>& comm_requirement) {
    int requirement = 0;
    for (int i = 0; i < comm_requirement.size(); i++) {
      requirement += comm_requirement[i];
    }
    return requirement;
  };


  int flexfly_topology::union_set_size(std::vector<std::vector<int>>& group_connectivity_vector, int start, int end) {
    bool appearance[num_optical_switches_];
    for (int i = 0; i < num_optical_switches_; i++) {
      appearance[i] = false;
    }
    int offset = num_groups_ * switches_per_group_;
    for (int index = start; index < end; index++) {
      for (auto optical_switch_id : group_connectivity_vector[index]) {
        appearance[optical_switch_id - offset] = true;
      }
    }
    int final_size = 0;
    for (int i = 0; i < num_optical_switches_; i++) {
      if (appearance[i]) {
        final_size++;
      }
    }
    return final_size;
  };

  bool flexfly_topology::check_remaining(std::vector<std::vector<int>>& group_connectivity_vector, std::vector<int>& requirement) {
    bool result = true;
    for (int i = 0 ; i < group_connectivity_vector.size(); i++) {
      if (group_connectivity_vector[i].size() < requirement[i]) {
        result = false;
        break;
      }
    }
    return result;
  };



  /**
   * Note: total_switches contains the number of optical switches as well
   * Note: unordered_map will have to be initialized by the caller
   **/
  void flexfly_topology::configure_optical_switches_general4(std::vector<std::vector<int>>& connectivity_requirement,
                                                            std::unordered_map<int, std::vector<int>>& optical_switch_inout_config) {
    /**********************************************************************************************************
    ******************************** WORKING SECTION **********************************************************
    **********************************************************************************************************/
    // vector or lists for a group, each element in list is an 
    // optical switch's id that this group is connected to
    int total_switches = num_optical_switches_ + num_groups_ * switches_per_group_;
    
    // entry i,j is the set of all optical switches that can still connect group i to j
    std::vector<std::vector<std::vector<int>>> group_physical_connectivity(num_groups_); 
    
    /*
     * Initialize the group_connectivity_matrix to all zeros
     */
    for (int i = 0 ; i < num_groups_; i++) {
      group_physical_connectivity[i].resize(num_groups_);
    }

    std::vector<std::vector<int>> group_nodes(num_groups_);
    std::unordered_map<int, std::vector<int>> optical_switch_nodes;  
    /*
     * Use vectors/list to emulate a graph node approach to solving this problem
     * First initiate the flexfly topology as a 3 layer graph. For each group node,
     * figure out all the all the 
     */
    for (int group = 0; group < num_groups_; group++) {
      int group_offset = group * switches_per_group_;
      for (int switch_index = 0; switch_index < switches_per_group_; switch_index++) {
        int curr_switch = group_offset + switch_index;
        auto switch_link_vector_iter = switch_outport_connection_map_.find(curr_switch);
        assert(switch_link_vector_iter != switch_outport_connection_map_.end());
        const std::vector<switch_link *>& switch_link_vector = switch_link_vector_iter->second;
        for (auto sl : switch_link_vector) {
          if (sl->type == Link_Type::electrical) continue;
          // dest_sid has to belong to an optical switch if topology wiring is done correctly
          assert(is_optical_switch(sl->dest_sid)); 
          group_nodes[group].insert(group_nodes[group].end(),  sl->dest_sid);
        }
      }
    }
    /*
     * Now do the same for the second layer of optical switches
     */
    int optical_switch_offset = num_groups_ * switches_per_group_;
    for (int opt_switch = 0; opt_switch < num_optical_switches_; opt_switch++) {
      int curr_switch = optical_switch_offset + opt_switch; // CURRENT OPTICAL SWITCH
      auto switch_link_vector_iter = switch_outport_connection_map_.find(curr_switch);
      const std::vector<switch_link *>& switch_link_vector = switch_link_vector_iter->second;
      optical_switch_nodes.insert(std::pair<int, std::vector<int>> (curr_switch, std::vector<int>()));
      for (auto sl : switch_link_vector) {
        optical_switch_nodes[curr_switch].push_back(group_from_swid(sl->dest_sid));
      }
    }

    /*
     * Fill up the group_physical_connectivity matrix of sets
     */
    for (int switch_id = switches_per_group_ * num_groups_; switch_id < total_switches; switch_id++) {
      auto outport_vector_iter = switch_outport_connection_map_.find(switch_id);
      auto inport_vector_iter = switch_inport_connection_map_.find(switch_id);
      assert(outport_vector_iter != switch_outport_connection_map_.end());
      auto outport_vector = outport_vector_iter->second;
      auto inport_vector = inport_vector_iter->second;
      for (auto sl1 : outport_vector) {
        assert(!is_optical_switch(sl1->dest_sid));
        int j = group_from_swid(sl1->dest_sid);
        for (auto sl2 : inport_vector) {
          int i = group_from_swid(sl2->src_sid);
          if (std::find(group_physical_connectivity[i][j].begin(), 
                        group_physical_connectivity[i][j].end(), 
                        switch_id) == group_physical_connectivity[i][j].end()) {
            group_physical_connectivity[i][j].push_back(switch_id);
          }
        }
      }
    }

    std::cout << "Printing set" << std::endl;
    for (auto x : group_physical_connectivity[0][1]) {
      std::cout << std::to_string(x) + " ";
    }
    std::cout << std::endl;


    /*
     * Main body of the algorithm
     */
    for (int i = 0; i < num_groups_; i++) {
      for (int j = 0; j < num_groups_; j++) {
        std::cout << "i, j is: " + std::to_string(i) + " " + std::to_string(j) << std::endl;
        while (connectivity_requirement[i][j] > 0) {
          // Find the potential optical switch to use
          bool use_potential_optical_switch = false;

          std::cout << "Printing set" << std::endl;
          for (auto x : group_physical_connectivity[i][j]) {
            std::cout << std::to_string(x) + " ";
          }
          std::cout << std::endl;

          for (auto potential_optical_switch : group_physical_connectivity[i][j]) {
            if (use_potential_optical_switch) break;
            
            std::cout << "potential_optical_switch is: " + std::to_string(potential_optical_switch) << std::endl;
            assert(std::find(optical_switch_nodes[potential_optical_switch].begin(), 
                            optical_switch_nodes[potential_optical_switch].end(), 
                            j) != optical_switch_nodes[potential_optical_switch].end());
            assert(std::find(group_nodes[i].begin(), 
                            group_nodes[i].end(), 
                            potential_optical_switch) != group_nodes[i].end());

            std::vector<std::pair<int, int>> entries_affected; // the (i,j) entries that are affected
            bool pass_check = true;


            for (int ii = 0; ii < num_groups_; ii++) {
              delete_item(group_physical_connectivity[ii][j], potential_optical_switch);
              entries_affected.push_back(std::make_pair(ii,j));
            }
            for (int ii = 0; ii < num_groups_; ii++) {
              delete_item(group_physical_connectivity[i][ii], potential_optical_switch);
              entries_affected.push_back(std::make_pair(i,ii));
            }
            connectivity_requirement[i][j]--;
            
            for (auto sl : switch_inport_connection_map_[potential_optical_switch]) {
              int src_group = group_from_swid(sl->src_sid);
              if (src_group == i) continue;
              if (std::find(group_nodes[src_group].begin(), 
                            group_nodes[src_group].end(),
                            potential_optical_switch) != group_nodes[src_group].end() && src_group != i) {
                if (union_set_size(group_physical_connectivity[src_group], 0, num_groups_) < remaining_row_requirement(connectivity_requirement[src_group])
                      || !check_remaining(group_physical_connectivity[src_group], connectivity_requirement[src_group])) {
                  std::cout << "cibai here la for dshGI;Herwgiuaerhugi src_group: " + std::to_string(src_group) << std::endl ;
                  pass_check = false;
                  break;
                } 
              }
            }

            if (pass_check) {
              
            } else {
              for (auto p : entries_affected) {
                int row = p.first;
                int col = p.second;
                add_item(group_physical_connectivity[row][col], potential_optical_switch);
              }
            }
            bool second_checkpoint = true;
            if (union_set_size(group_physical_connectivity[i], 0, num_groups_) < remaining_row_requirement(connectivity_requirement[i])
                    || !check_remaining(group_physical_connectivity[i], connectivity_requirement[i])) {
              //std::cout << "cibai here la for dst_group: " + std::to_string(target_group) << std::endl;
              pass_check = false;
              second_checkpoint = false;
              //break;
            }

            if (pass_check) {
              use_potential_optical_switch = true;
              assert(delete_item(group_physical_connectivity[i][j], potential_optical_switch));
              assert(delete_item(group_nodes[i], potential_optical_switch));
              assert(delete_item(optical_switch_nodes[potential_optical_switch], j));
              for (int ii = 0; ii < num_groups_; ii++) {
                delete_item(group_physical_connectivity[i][ii], potential_optical_switch);
                
              }
              connectivity_requirement[i][j]--;
              break;
            } else if (!pass_check && !second_checkpoint) {
              for (auto p : entries_affected) {
                int row = p.first;
                int col = p.second;
                add_item(group_physical_connectivity[row][col], potential_optical_switch);
              }
            }
            
          }

          //assert(use_potential_optical_switch);
          if (!use_potential_optical_switch) {
            std::cout << "FAILED FOR i: " + std::to_string(i) + " j: " + std::to_string(j) << std::endl; 
            connectivity_requirement[i][j]--;
            assert(use_potential_optical_switch); 
          }
        }
      }
    }
  /**********************************************************************************************************
  ******************************** WORKING SECTION **********************************************************
  **********************************************************************************************************/
    return;
  };


}
}

