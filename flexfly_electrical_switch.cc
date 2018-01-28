#include "flexfly_electrical_switch.h"
//#include "flexfly_packet.h"
#include <sstmac/hardware/common/connection.h>
#include <sstmac/hardware/topology/topology.h>
#include <sprockit/sim_parameters.h>
#include <sprockit/util.h>
#include <sstmac/hardware/router/router_fwd.h>
#include <sstmac/hardware/switch/network_switch.h>
#include <sprockit/factories/factory.h>
#include <sstmac/hardware/router/minimal_routing.h>
#include <sstmac/hardware/pisces/pisces.h>
#include <sstmac/software/launch/launch_event.h>
#include <sstmac/hardware/network/network_message.h>

namespace sstmac {
namespace hw {

	flexfly_electrical_switch::flexfly_electrical_switch(sprockit::sim_parameters* params,
    														uint64_t id,
    														event_manager* mgr,
    														device_id::type_t ty) : 
																network_switch(params, 
																				id,
																				mgr, 
																				device_id::logp_overlay) {
		my_addr_ = params->get_int_param("id");
		radix_ = params->get_int_param("total_radix");
		switches_per_group_ = params->get_int_param("switches_per_group");
		num_groups_ = params->get_int_param("num_groups");
		nodes_per_switch_ = radix_ - switches_per_group_;
		inport_handlers_.reserve(radix_);
		outport_handlers_.reserve(radix_);
		queue_length_ = new int[radix_];

		/*
		 * Grab the link parameters like electrical and optical link bandwidths
		 * and even send and credit latencies
		 */
		get_link_params(params);
		/*
		 * Grab the link parameters like electrical and optical link bandwidths
		 * and even send and credit latencies
		 */

		sprockit::sim_parameters* rtr_params = params->get_optional_namespace("router");
		rtr_params->add_param_override_recursive("id", int(my_addr_));
		router_ = router::factory::get_param("name", rtr_params, top_, this);
		ftop_ = dynamic_cast<flexfly_topology *>(topology::static_topology(params));
		if (ftop_ == nullptr) {
			ftop_simplified_ = safe_cast(flexfly_topology_simplified, topology::static_topology(params));
		} else {
			ftop_simplified_ = nullptr;
		}
		init_links(params);
	}

	flexfly_electrical_switch::~flexfly_electrical_switch() {
		delete [] queue_length_;
	}

	/**
	 * Obtain the link parameters
	 **/
	void flexfly_electrical_switch::get_link_params(sprockit::sim_parameters* switch_params) {
		sprockit::sim_parameters* link_params = switch_params->get_namespace("link");
		if (link_params == nullptr) {
			return;
		}
		send_latency_ = link_params->get_time_param("send_packet_latency");
		credit_latency_ = link_params->get_time_param("credit_latency");
		electrical_link_bw_ = link_params->get_bandwidth_param("electrical_link_bandwidth");
		optical_link_bw_ = link_params->get_bandwidth_param("optical_link_bandwidth");
		inv_electrical_link_bw_ = 1 / electrical_link_bw_;
		inv_optical_link_bw_ = 1 / optical_link_bw_;
	}

	int flexfly_electrical_switch::queue_length(int port) const {
		//std::cout << "queue length function is called" << std::endl;
		return queue_length_[port];
	}

	void flexfly_electrical_switch::connect_input(sprockit::sim_parameters* params, 
                              						int src_outport, 
                              						int dst_inport,
                              						event_handler* credit_handler) {
		//if (dst_inport < 0 || dst_inport >= radix_)
		//	spkt_abort_printf("Invalid inport %d in flexfly_electrical_switch::connect_input", dst_inport);
		
		inport_handlers_[dst_inport] = credit_handler;
	}

	void flexfly_electrical_switch::connect_output(sprockit::sim_parameters* params, 
                              						int src_outport, 
                              						int dst_inport,
                              						event_handler* payload_handler) {
		//if (src_outport < 0 || src_outport >= radix_)
		//	spkt_abort_printf("Invalid inport %d in flexfly_electrical_switch::connect_output", src_outport);

		outport_handlers_[src_outport] = payload_handler;
	}

	link_handler* flexfly_electrical_switch::credit_handler(int port) const {
		if (port < switches_per_group_) {
			return new_link_handler(this, &flexfly_electrical_switch::recv_credit);
		} else {
			return new_link_handler(this, &flexfly_electrical_switch::recv_nodal_credit);
		}
	}

	link_handler* flexfly_electrical_switch::payload_handler(int port) const {
		if (port < switches_per_group_) {
			return new_link_handler(this, &flexfly_electrical_switch::recv_payload);
		} else {
			return new_link_handler(this, &flexfly_electrical_switch::recv_nodal_payload);
		}
	}

	/**
	 * Receiving a packet or message from another switch, not from a node.
	 **/
	void flexfly_electrical_switch::recv_payload(event* ev) {
		
		if (ftop_ != nullptr) {
			flexfly_packet* fpacket = safe_cast(flexfly_packet, ev);
			node_id dst = fpacket->get_pisces_packet()->toaddr();
  			node_id src = fpacket->get_pisces_packet()->fromaddr();
  			// Case 1: route it to a connecting node
  			//std::cout << "electrical switch received a payload" << std::endl;
  			switch_id dst_swid = ftop_->node_to_switch(dst);
  			if (dst_swid == my_addr_) {
  				int outport = dst - (my_addr_ * nodes_per_switch_) + switches_per_group_;
  				send_to_link(outport_handlers_[outport], fpacket->get_pisces_packet());
  				return;
  			}
  			send_to_link(outport_handlers_[fpacket->next_outport()], fpacket);
  		} else {
  			pisces_default_packet* msg = safe_cast(pisces_default_packet, ev);
  			int num_bits = msg->num_bytes();
  			timestamp electrical_delay = num_bits * inv_electrical_link_bw_ + send_latency_;
			timestamp optical_delay = num_bits * inv_optical_link_bw_ + send_latency_;
  			int dst = msg->toaddr();
			int src = msg->fromaddr();
			int dst_switch = ftop_simplified_->node_to_switch(dst);
			int src_switch = ftop_simplified_->node_to_switch(src);
			int dst_group = ftop_simplified_->group_from_swid(dst_switch);
			int src_group = ftop_simplified_->group_from_swid(src_switch);
			if (dst_switch == my_addr_) {
				int offset = dst % nodes_per_switch_;
				send_delayed_to_link(electrical_delay, outport_handlers_[offset + switches_per_group_], msg);
				//std::cout << "Some message actually got delivered to its destination" << std::endl;
			} else if (dst_group == ftop_simplified_->group_from_swid(my_addr_)) {
				int outport = ftop_simplified_->get_output_port(my_addr_, dst_switch);
				assert(outport >= 0);
				send_delayed_to_link(electrical_delay, outport_handlers_[outport], ev);
			} else {
				assert(dst_group != ftop_simplified_->group_from_swid(my_addr_));
				send_delayed_to_link(optical_delay, outport_handlers_[switches_per_group_ - 1], msg);
			}
  		}
	}

	/**
	 * Receiving a packet or message from a node, not from a switch.
	 * This can be said to be an injection packet. Must be routed through
	 * the network 
	 **/
	void flexfly_electrical_switch::recv_nodal_payload(event* ev) {
		// need to generate a new flexfly_packet here
		pisces_default_packet* msg = safe_cast(pisces_default_packet, ev);
		//std::cout << "received_nodal_payload?" << std::endl;
		//std::cout << "From node: " << std::to_string(msg->fromaddr()) << " to node: " + std::to_string(msg->toaddr()) << std::endl;
		int dst = msg->toaddr();
		int src = msg->fromaddr();
		int dst_switch = 0;
		int src_switch = 0;
		//std::cout << "received an message at electrical_switch" << std::endl;
		int num_bits = 8 * msg->num_bytes();
		timestamp electrical_delay = num_bits * inv_electrical_link_bw_ + send_latency_;
		timestamp optical_delay = num_bits * inv_optical_link_bw_ + send_latency_;

		if (ftop_ != nullptr) {
			dst_switch = ftop_->node_to_switch(dst);
			src_switch = ftop_->node_to_switch(src);
		} else {
			dst_switch = ftop_simplified_->node_to_switch(dst);
			src_switch = ftop_simplified_->node_to_switch(src);
		}
		
		/*
		 * Case 1: When the packet has reached its final switch
		 */
		if (src_switch == dst_switch) {
			//send_credit_to_node(msg->num_bytes(), src);
			send_credit_to_node(1000000, src);
			send_packet_to_node(dst, electrical_delay, ev);
			return; 
		}

		// If we are here, that means that we do need to route properly
		
		if (ftop_ != nullptr) {
			flexfly_packet* fpacket = new flexfly_packet(msg);
			ftop_->route_minimal(my_addr_, ftop_->node_to_switch(msg->toaddr()), fpacket);
			int next_port = fpacket->next_outport();
			//std::cout << "Next port is: " + std::to_string(next_port) << std::endl;
			send_delayed_to_link(electrical_delay, outport_handlers_[next_port], fpacket);
		} else { 
			// In this case, we are in the simplified
			int src_group = ftop_simplified_->group_from_swid(src_switch);
			int dst_group = ftop_simplified_->group_from_swid(dst_switch);
			routable::path new_path;
			ftop_simplified_->minimal_route_to_switch(my_addr_, dst_switch, new_path);
			int outport = new_path.outport();
			send_delayed_to_link(electrical_delay, outport_handlers_[outport], ev);
			return;
			/*
			if (src_group == dst_group) {
				int outgoing_port = ftop_simplified_->get_output_port(src_switch, dst_switch);
				assert(outgoing_port >= 0);
				send_delayed_to_link(electrical_delay, outport_handlers_[outgoing_port], ev);
			} else {
				// If not in the same group, then just send it over to the optical_network
				send_delayed_to_link(optical_delay, outport_handlers_[switches_per_group_ - 1], ev);
			}
			*/
		}
		
		
	};

	void flexfly_electrical_switch::recv_nodal_credit(event* ev) {
		
		// need to generate a new flexfly_packet here
		pisces_credit* msg = safe_cast(pisces_credit, ev);
		//std::cout << "From node: " << std::to_string(msg->fromaddr()) << " to node: " + std::to_string(msg->toaddr()) << std::endl;
		//std::cout << "The credits port is " << std::to_string(msg->port()) << std::endl;
		//std::cout << "The num credit is " << std::to_string(msg->num_credits()) << std::endl;
		//std::cout << "received an message at electrical_switch" << std::endl;
		send_delayed_to_link(credit_latency_, inport_handlers_[msg->port()], ev);
		//delete msg;
		return;
	};

	void flexfly_electrical_switch::recv_credit(event* ev) {
		pisces_credit* msg = safe_cast(pisces_credit, ev);
		msg->port();
		//std::cout << "received an message at electrical_switch" << std::endl;
		send_delayed_to_link(credit_latency_, inport_handlers_[msg->port()], ev);
		return;
	};

	//void flexfly_electrical_switch::send_credit(event_handler* ev, pisces_credit* pc) {

	//}

	/**
	 * This function sends a delayed packet to the node connected to this switch with node_id 
	 * as the node id
	 **/
	void flexfly_electrical_switch::send_packet_to_node(int node_id, timestamp delay, event* ev) {
		assert(ftop_simplified_->node_to_switch(node_id) == my_addr_);
		int offset = node_id % nodes_per_switch_;
		int outport = offset + switches_per_group_;
		send_delayed_to_link(delay, outport_handlers_[outport], ev);
	}

	void flexfly_electrical_switch::send_credit_to_node(uint32_t credit_amount, int node_id) {
		assert(ftop_simplified_->node_to_switch(node_id) == my_addr_);
		int port = switches_per_group_ + (node_id % nodes_per_switch_);
		pisces_credit* pc = new pisces_credit(port, 0, credit_amount);
		send_delayed_to_link(send_latency_ , inport_handlers_[port], pc);
	}
}
}