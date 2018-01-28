#include "flexfly_optical_switch.h"
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
	flexfly_optical_switch::flexfly_optical_switch(sprockit::sim_parameters* params,
    												uint64_t id,
    												event_manager* mgr) : optical_switch(params, id, mgr) {
		my_addr_ = params->get_int_param("id");
		num_ports_ = params->get_int_param("optical_switch_radix");
		num_electrical_switches_ = params->get_int_param("num_electrical_switches");
		inport_handler_.reserve(num_ports_);
		outport_handler_.reserve(num_ports_);
		inout_connection_ = new int[num_ports_];
		for (int i = 0; i < num_ports_; i++) {
			inout_connection_[i] = i;
		}
		top_ = safe_cast(flexfly_topology, topology::static_topology(params));
		init_links(params); // this has to be called upon class initialization
		std::vector<int> tmp;
		top_->optical_switch_configuration(my_addr_ - num_electrical_switches_, tmp);
		for (int i = 0; i < num_ports_; i++) {
			inout_connection_[i] = tmp[i];
		}

		std::cout << "switch: " << std::to_string(my_addr_) << std::endl;

		for (int i = 0; i < num_ports_; i++) {
			std::cout << "inport: " << std::to_string(i) << " -> outport: " << std::to_string(inout_connection_[i]) << std::endl;
		}
	};

	flexfly_optical_switch::~flexfly_optical_switch() {
		delete [] inout_connection_;
	}

	void flexfly_optical_switch::setup() {
		return; 
	};

	/**
	 * NOTE: This method is only called once, and it hooks the input port to a handler 
	 * 			so that when the switch receives an 
	 */
	void flexfly_optical_switch::connect_input(sprockit::sim_parameters* params, 
                              					int src_outport, 
                              					int dst_inport,
                              					event_handler* credit_handler) {
		if (dst_inport < 0 || dst_inport >= num_ports_)
			spkt_abort_printf("Invalid inport %d in flexfly_optical_switch::connect_input", dst_inport);
 		inport_connections_[dst_inport] = src_outport;
 		int input_port = dst_inport;
 		inport_handler_[input_port] = credit_handler;
		return;
		
	};

	void flexfly_optical_switch::connect_output(sprockit::sim_parameters* params, 
                              					int src_outport, 
                              					int dst_inport,
                              					event_handler* payload_handler) {
		//std::cout << "src_outport: " << std::to_string(src_outport) << " and dst_inport: " << std::to_string(dst_inport) << std::endl;
		//std::cout << "src_outport: " << std::to_string(src_outport) << " and dst_inport: " << std::to_string(dst_inport) << std::endl;
 		if (src_outport < 0 || src_outport >= num_ports_)
			spkt_abort_printf("Invalid inport %d in flexfly_optical_switch::connect_output", src_outport);
 		outport_connections_[src_outport] = dst_inport;
 		int output_port = src_outport;
 		outport_handler_[output_port] = payload_handler;
		return;
	};

	link_handler* flexfly_optical_switch::payload_handler(int port) const {
		return new_link_handler(this, &flexfly_optical_switch::recv_payload);
	};

	link_handler* flexfly_optical_switch::credit_handler(int port) const {
		return new_link_handler(this, &flexfly_optical_switch::recv_credit);
	};

	void flexfly_optical_switch::deadlock_check() {
		return;
	};

	void flexfly_optical_switch::deadlock_check(event* ev) {
		return;
	};

	void flexfly_optical_switch::recv_credit(event* ev) {
		std::cout << "Optical Switch received a credit" << std::endl;
		return;
	};

	void flexfly_optical_switch::recv_payload(event* ev) {
		//message* msg = safe_cast(message, ev);
		//std::cout << "Optical Switch received a payload" << std::endl;

		flexfly_packet *fpacket = safe_cast(flexfly_packet, ev);
		assert(fpacket);
		pisces_payload* msg = fpacket->get_pisces_packet();
		int curr_outport = fpacket->next_outport();
		send_to_link(outport_handler_[curr_outport], fpacket);
		/*
		pisces_payload* msg = safe_cast(pisces_payload, ev);
		node_id dst_id = msg->toaddr();
		node_id src_id = msg->fromaddr();
		switch_id src_swid = top_->node_to_switch(src_id);
		
		int src_group = top_->group_from_swid(src_swid);
		int inport = src_group;
		int outport = inout_connection_[inport];
		//std::cout << "inport is: " << std::to_string(inport) << " and outport is: " << std::to_string(outport) << std::endl; 
		//std::cout << "src_node is: " << std::to_string(src_id) << " and dst_node is: " << std::to_string(dst_id) << std::endl; 
		//std::cout << "src_swid is: " << std::to_string(src_swid) << std::endl;
		send_to_link(outport_handler_[outport],ev);
		*/
	};

	bool flexfly_optical_switch::outport_connected(int outport) const {
		for (int i = 0; i < num_ports_; i++) {
			if (inout_connection_[i] == outport) {
				return true;
			}
		}
		return false;
	}

	bool flexfly_optical_switch::setup_inout_connection(int inport, int outport) {
		// check if the current outport is already connected to any of the inports
		if (outport_connected(outport)) {
			return false;
		}
		// schedule some sort of delay?
		teardown_outport_connection(outport);
		inout_connection_[inport] = outport; 
		return true;
	}

	void flexfly_optical_switch::teardown_outport_connection(int outport) {
		for (int i = 0; i < num_ports_; i++) {
			if (inout_connection_[i] == outport) {
				inout_connection_[i] = -1;
			}	
		}
	}
	/*
	void flexfly_optical_switch::recv_config_msg(event* ev) {
		optical_configuration* config = safe_cast(optical_configuration, ev);

	}
	*/
}
}