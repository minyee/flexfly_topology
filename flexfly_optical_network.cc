#include "flexfly_optical_network.h"
#include "flexfly_events.h"
#include <sstmac/hardware/common/connection.h>
#include <iostream>
#include <stdio.h>
#include <sstmac/hardware/pisces/pisces.h>
#include <sprockit/util.h>

//#include <sstmac/software/launch/launch_event.h>
#include <sstmac/hardware/network/network_message.h>


namespace sstmac {
namespace hw {
	flexfly_optical_network::flexfly_optical_network(sprockit::sim_parameters* params,
    												uint64_t id,
    												event_manager* mgr) : optical_switch(params, id, mgr) {
		my_addr_ = params->get_int_param("id");
		num_electrical_switches_ = params->get_int_param("num_electrical_switches");
		num_groups_ = params->get_int_param("num_groups");
		num_ports_ = num_groups_ * num_electrical_switches_;
		switches_per_group_ = params->get_int_param("switches_per_group");
		inport_handler_.reserve(num_ports_);
		outport_handler_.reserve(num_ports_);
		//ftop_ = safe_cast(flexfly_topology_simplified, topology::static_topology(params));
		ftop_ = dynamic_cast<flexfly_topology *>(topology::static_topology(params));
		if (ftop_ == nullptr) {
			ftop_simplified_ = safe_cast(flexfly_topology_simplified, topology::static_topology(params));
			outport_options_.resize(num_groups_);
			for (int i = 0; i < num_groups_; i++) {
				outport_options_[i].resize(num_groups_);
			}
			ftop_simplified_->configure_optical_network(outport_options_);
		} else {
			ftop_simplified_ = nullptr;
		}

		init_links(params); // this has to be called upon class initialization
	};

	flexfly_optical_network::~flexfly_optical_network() { 
		// Do nothing for now since all private member variables are allocated on the stack
	}

	void flexfly_optical_network::setup() {
		return; 
	};

	/**
	 * NOTE: This method is only called once, and it hooks the input port to a handler 
	 * 			so that when the switch receives an 
	 */
	void flexfly_optical_network::connect_input(sprockit::sim_parameters* params, 
                              					int src_outport, 
                              					int dst_inport,
                              					event_handler* credit_handler) {
 		int input_port = dst_inport;
 		inport_handler_[input_port] = credit_handler;
 		//std::cout << "is flexfly_optical_network's connect_input function ever called? " << std::endl;
		return;
		
	};

	void flexfly_optical_network::connect_output(sprockit::sim_parameters* params, 
                              					int src_outport, 
                              					int dst_inport,
                              					event_handler* payload_handler) {
		//std::cout << "src_outport: " << std::to_string(src_outport) << " and dst_inport: " << std::to_string(dst_inport) << std::endl;
		//std::cout << "src_outport: " << std::to_string(src_outport) << " and dst_inport: " << std::to_string(dst_inport) << std::endl;
 		if (src_outport < 0 || src_outport >= num_ports_)
			spkt_abort_printf("Invalid inport %d in flexfly_optical_network::connect_output", src_outport);
 		int output_port = src_outport;
 		outport_handler_[output_port] = payload_handler;
		return;
	};

	link_handler* flexfly_optical_network::payload_handler(int port) const {
		return new_link_handler(this, &flexfly_optical_network::recv_payload);
	};

	link_handler* flexfly_optical_network::credit_handler(int port) const {
		return new_link_handler(this, &flexfly_optical_network::recv_credit);
	};

	void flexfly_optical_network::deadlock_check() {
		return;
	};

	void flexfly_optical_network::deadlock_check(event* ev) {
		return;
	};

	void flexfly_optical_network::recv_credit(event* ev) {
		std::cout << "Optical optical network received a credit" << std::endl;
		return;
	};

	void flexfly_optical_network::recv_payload(event* ev) {
		//message* msg = safe_cast(message, ev);
		//std::cout << "Optical Network received a payload" << std::endl;

		pisces_default_packet *packet = safe_cast(pisces_default_packet, ev);
		assert(packet);
		//pisces_payload* msg = fpacket->get_pisces_packet();
		int dst_node = packet->toaddr();
		int src_node = packet->fromaddr();
		//ftop_->group_from_swid(dst_switch);
		if (ftop_ == nullptr) {
			int dst_switch = ftop_simplified_->node_to_switch(dst_node);
			int src_switch = ftop_simplified_->node_to_switch(src_node);
			int outport;
			ftop_simplified_->minimal_route_to_switch_optical(src_switch, dst_switch, outport);
			send_to_link(outport_handler_[outport], ev);
		}
		
		//send_to_link(outport_handler_[], fpacket);

	};

	void flexfly_optical_network::recv_payload_smart(event* ev) {
		pisces_default_packet *packet = safe_cast(pisces_default_packet, ev);
		assert(packet);
		//pisces_payload* msg = fpacket->get_pisces_packet();
		int dst_node = packet->toaddr();
		int src_node = packet->fromaddr();
		//std::cout << "This packet goes from source node id of: " << std::to_string(src_node) << " to node id: " << std::to_string(dst_node) << std::endl;
		if (ftop_ == nullptr) {
			int dst_switch = ftop_simplified_->node_to_switch(dst_node);
			int src_switch = ftop_simplified_->node_to_switch(src_node);
			int dst_group = ftop_simplified_->group_from_swid(dst_switch);
			int src_group = ftop_simplified_->group_from_swid(src_switch);
			int outport = ftop_simplified_->get_output_port(my_addr_, dst_switch);	
			if (std::find(outport_options_[src_group][dst_group].begin(), outport_options_[src_group][dst_group].end(), outport) == outport_options_[src_group][dst_group].end()) {
				int size = outport_options_[src_group][dst_group].size();
				uint32_t seed = 12;
   				std::srand(seed);
   				int choice = rand() % size;
   				outport = outport_options_[src_group][dst_group][choice];
			}
			send_to_link(outport_handler_[outport], ev);
		}
	};
	/*
	void flexfly_optical_network::recv_config_msg(event* ev) {
		optical_configuration* config = safe_cast(optical_configuration, ev);

	}
	*/
}
}