
/*
 *  This file is part of SST/macroscale:
 *               The macroscale architecture simulator from the SST suite.
 *  Copyright (c) 2009 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the LICENSE file in the top
 *  SST/macroscale directory.
 */

#include <sstmac/hardware/switch/network_switch.h>
//#include <sstmac/hardware/logp/my_logp_switch.h>
#include "my_logp_switch.h"
#include <sstmac/hardware/nic/nic.h>
#include <sstmac/hardware/switch/dist_dummyswitch.h>
#include <sstmac/hardware/topology/topology.h>
#include <sstmac/hardware/interconnect/interconnect.h>
#include <sstmac/common/event_manager.h>
#include <sprockit/util.h>
#include <sprockit/sim_parameters.h>
#include <sstmac/software/launch/launch_event.h>

namespace sstmac {
namespace hw {

my_logp_switch::my_logp_switch(sprockit::sim_parameters *params, uint64_t id, event_manager *mgr) :
  network_switch(params, id, mgr, device_id::logp_overlay)
{
  //sprockit::sim_parameters* link_params = params->get_namespace("link");
  //sprockit::sim_parameters* ej_params = params->get_namespace("ejection"); 
  int my_id = params->get_int_param("id");
  electrical_bw_ = params->get_bandwidth_param("electrical_bandwidth");
  optical_bw_ = params->get_bandwidth_param("optical_bandwidth");

  inv_optical_bw_ = 1.0/optical_bw_;
  inv_electrical_bw_ = 1.0/electrical_bw_;

  //double inj_bw = ej_params->get_optional_bandwidth_param("bandwidth", net_bw);
  double inj_bw = electrical_bw_;
  inj_bw_inverse_ = 1.0/inj_bw;
  /*
  if (ej_params->has_param("send_latency")){
    inj_lat_ = ej_params->get_time_param("send_latency");
  } else {
    inj_lat_ = ej_params->get_time_param("latency");
  }
  */
  dbl_inj_lat_ = 2 * inj_lat_;

  interconn_ = interconnect::static_interconnect(params, mgr);

  nics_.resize(top_->num_nodes());
  neighbors_.reserve(1000); //nproc - just reserve a large block for now
  ftop_ = dynamic_cast<flexfly_topology *>(topology::static_topology(params));
  if (ftop_ == nullptr) {
    ftop_simplified_ = safe_cast(flexfly_topology_simplified, topology::static_topology(params));
  } else {
    ftop_simplified_ = nullptr;
  }

  if (ftop_) {
    switches_per_group_ = ftop_->switches_per_group();
    nodes_per_switch_ = ftop_->nodes_per_switch();
  } else {
    switches_per_group_ = ftop_simplified_->switches_per_group();
    nodes_per_switch_ = ftop_simplified_->nodes_per_switch();
  }
  init_links(params);
  std::cout << "Exited my_logp_switch constructor" << std::endl;
}

my_logp_switch::~my_logp_switch()
{
#if !SSTMAC_INTEGRATED_SST_CORE
  delete mtl_handler_;
#endif
}

link_handler*
my_logp_switch::payload_handler(int port) const
{
  return new SST::Event::Handler<my_logp_switch>(
        const_cast<my_logp_switch*>(this), &my_logp_switch::handle);
}

void
my_logp_switch::connect_output(sprockit::sim_parameters *params,
                              int src_outport, int dst_inport,
                              event_handler *mod) {
    node_id nid = src_outport;
    //std::cout << "What about now" << std::endl;
    switch_debug("Connecting my_LogP to NIC %d", nid);
    nics_[nid] = mod;
    switch_id sid = src_outport;
    switch_debug("Connecting to my_LogP switch %d", sid);
    if (sid >= neighbors_.size()){
      neighbors_.resize(sid+1);
    }
    neighbors_[sid] = mod;
}

void
my_logp_switch::connect_input(sprockit::sim_parameters *params,
                              int src_outport, int dst_inport,
                              event_handler *mod)
{
  return;
}

void
my_logp_switch::incoming_message(message* msg, node_id src, node_id dst)
{
  bool skip = false;
  bool local_src = nics_[src];
  long num_bytes = msg->byte_length();
  timestamp delay = 2 * inj_bw_inverse_ * num_bytes; //bw term
  switch_id src_switch = src / nodes_per_switch_;
  switch_id dst_switch = dst / nodes_per_switch_;
  int src_group = src_switch / switches_per_group_;
  int dst_group = dst_switch / switches_per_group_;
  if (!skip) {
    if (src_switch == dst_switch) {
    } else if (src_group == dst_group) {
      delay += (inv_electrical_bw_) * num_bytes; 
    } else {
      int num_links = 0;
      if (ftop_ == nullptr) {
        num_links = ftop_simplified_->num_links_between_groups(src_group, dst_group);
      } else {
        num_links = ftop_->num_links_between_groups(src_group, dst_group);
      } 
      double net_intergroup_bw = num_links * optical_bw_;
      delay += (num_links * inv_optical_bw_ + 2 * inv_electrical_bw_) * num_bytes;
    }
    send_delayed_to_link(delay, nics_[dst], msg);
  } else {
    /*
        switch_debug("incoming message over %d hops with extra delay %12.8e and inj lat %12.8e: %s",
               num_hops, delay.sec(), dbl_inj_lat_.sec(), msg->to_string().c_str());
               */
    send_to_link(nics_[dst], msg);
  }
  //
}

// TODO: Need to modify this, the delayed to link function has to include bandwidth from group to group
void
my_logp_switch::outgoing_message(message* msg, node_id src, node_id dst)
{
  //local going remote - just accumulate latency delay
  //int num_hops = ftop_->num_hops_to_node(src, dst);
  switch_id src_switch = src / nodes_per_switch_;
  switch_id dst_switch = dst / nodes_per_switch_;
  int src_group = src_switch / switches_per_group_;
  int dst_group = dst_switch / switches_per_group_;
  timestamp delay; //factor of 2 for in-out
  if (src_switch == dst_switch) {
    delay = 2 * inj_bw_inverse_;
  } else if (src_group == dst_group) {
    delay = 2 * inj_bw_inverse_ + inv_electrical_bw_; 
  } else {
    int num_links = 0;
    if (ftop_ == nullptr) {
      num_links = ftop_simplified_->num_links_between_groups(src_group, dst_group);
    } else {
      num_links = ftop_->num_links_between_groups(src_group, dst_group);
    } 
    double net_intergroup_bw = num_links * optical_bw_;
    delay = 2 * inj_bw_inverse_ + num_links * inv_optical_bw_ + 2 * inv_electrical_bw_;
  }
  


  //int dst_switch = interconn_->node_to_logp_switch(dst);
  /*
  switch_debug("outgoing message over %d hops with extra delay %12.8e and inj lat %12.8e: %s",
               num_hops, delay.sec(), dbl_inj_lat_.sec(), msg->to_string().c_str());
               */
  send_delayed_to_link(delay, dbl_inj_lat_, neighbors_[dst_switch], msg);
}

void
my_logp_switch::bcast_local_message(message* msg, node_id src)
{
  std::cout << "myLogP Broadcasting message" << std::endl;
  sw::start_app_event* lev = safe_cast(sw::start_app_event, msg);
  sw::task_mapping::ptr mapping = lev->mapping();
  int num_ranks = mapping->num_ranks();
  //std::cout << "Broadcasting message" << std::endl;
  long num_bytes = msg->byte_length();
  for (int i=0; i < num_ranks; ++i){
    node_id dst_node = mapping->rank_to_node(i);
    auto dst_nic = nics_[dst_node];
    bool local_dst = dst_nic;
    bool local_src = nics_[src];

    if (local_dst){
      sw::start_app_event* new_lev = lev->clone(i, src, dst_node);
      //int num_hops = top_->num_hops_to_node(src, dst_node);
      switch_id src_switch = src / nodes_per_switch_;
      switch_id dst_switch = dst_node / nodes_per_switch_;
      int src_group = src_switch / switches_per_group_;
      int dst_group = dst_switch / switches_per_group_;

      timestamp delay= 2 * inj_bw_inverse_ * num_bytes; //factor of 2 for in-out
      
      if (src_switch == dst_switch) {
        //delay = 2 * inj_bw_inverse_;
      } else if (src_group == dst_group) {
        delay += inv_electrical_bw_ * num_bytes; 
      } else {
        int num_links = 0;
        if (ftop_ == nullptr) {
          num_links = ftop_simplified_->num_links_between_groups(src_group, dst_group);
        } else {
          num_links = ftop_->num_links_between_groups(src_group, dst_group);
        } 
        double net_intergroup_bw = num_links * optical_bw_;
        delay += (num_links * inv_optical_bw_ + 2 * inv_electrical_bw_) * num_bytes;
      }
      //send_to_link(dst_nic, new_lev);
      send_delayed_to_link(delay, dst_nic, new_lev);
    }
  }
  //this one not needed anymore
  delete lev;
}

void
my_logp_switch::forward_bcast_message(message* msg, node_id dst)
{
  int dst_switch = interconn_->node_to_logp_switch(dst);
  //only accumulate inj lat - exact hop latency gets added on the other side
  send_delayed_to_link(timestamp(), dbl_inj_lat_, neighbors_[dst_switch], msg);
}

void
my_logp_switch::handle(event* ev)
{
  //this should only handle messages
  message* msg = safe_cast(message, ev);
  node_id dst = msg->toaddr();
  node_id src = msg->fromaddr();
  bool local_dst = nics_[dst];

  switch_debug("handling message %d->%d of size %ld: %s",
               src, dst, msg->byte_length(),
               msg->to_string().c_str());

  if (msg->is_bcast()){
    if (local_dst){
      bcast_local_message(msg, src);
    } else {
      forward_bcast_message(msg, dst);
    }
  } else if (local_dst){
    incoming_message(msg, src, dst);
  } else { //locall going remote
    outgoing_message(msg, src, dst);
  }
}



}
}