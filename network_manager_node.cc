#include "network_manager_node.h"
#include "connectivity_matrix.h"
#include <sstmac/software/launch/launch_event.h>
#include <sstmac/common/event_scheduler.h>
namespace sstmac {
namespace hw {
	network_manager_node::network_manager_node(sprockit::sim_parameters* params,
												 uint64_t id, event_manager* mgr) : node(params, id, mgr){
		id_ = params->get_int_param("id");
		num_groups_ = params->get_int_param("num_groups");
		num_optical_switches_ = num_groups_ - 1;
		outport_handler_.reserve(num_optical_switches_);
		inport_handler_.reserve(num_optical_switches_);
		init_links(params);
	};

	network_manager_node::~network_manager_node() {};

	void network_manager_node::connect_input(sprockit::sim_parameters* params,
												int src_outport,
												int dst_inport,
												event_handler* payload_handler) {
		inport_handler_[dst_inport] = payload_handler;
	};

	void network_manager_node::connect_output(sprockit::sim_parameters* params,
												int src_outport,
												int dst_inport,
												event_handler* payload_handler) {
		outport_handler_[src_outport] = payload_handler;
	};

	link_handler* network_manager_node::credit_handler(int port) const {
		return nullptr;
	};

	link_handler* network_manager_node::payload_handler(int port) const {
		return new_link_handler(this, &network_manager_node::receive_payload);
	};

	void network_manager_node::receive_payload(event* ev) {
		sw::start_app_event* lev = dynamic_cast<sw::start_app_event*>(ev);
		// ignore everything but the start events in logp
		if (lev == nullptr) {
			return;
		}

		sw::app_id aid = lev->aid();
		std::cout << "Application: " << std::to_string(aid) << " received a start event" << std::endl;
		// if we are here, it means that we did indeed receive a start event for an app
		
		/**
		* Run the configuration algorithm here and then get the proper inout port configuration
		* for each optical switch, and then send each one of them out.
		*
		**/
		std::vector<std::vector<int>> connectivity_matrix;
		//configure_simpler_model(num_groups_, connectivity_matrix);
		for (int i = 0 ; i < num_optical_switches_; i++) {
			std::vector<int>& optical_inout_connectivity = connectivity_matrix[i];
			optical_configuration_event* configuration = new optical_configuration_event();
			configuration->set_optical_configurations(optical_inout_connectivity);
			send_to_link(outport_handler_[i], configuration);
		}
	};

	void network_manager_node::execute(ami::COMP_FUNC func, event* data, callback* cb) {};

	optical_configuration_event::optical_configuration_event() {

	};
}	
}