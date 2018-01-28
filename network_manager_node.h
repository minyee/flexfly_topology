#include "flexfly_events.h"
#include <sstmac/hardware/node/node.h>
#include <vector>
namespace sstmac {
namespace hw{

class network_manager_node : public node {

RegisterComponent("",node,network_manager_node,
	"macro", COMPONENT_CATEGORY_NETWORK,
	"A network manager node that is used to manage the optical link configurations in Flexfly topologies");

public:
	network_manager_node(sprockit::sim_parameters* params, uint64_t id, event_manager* mgr);
	
	~network_manager_node();

	virtual void connect_output(sprockit::sim_parameters* params,
								int src_outport,
								int dst_inport,
								event_handler* payload_handler) override;

	virtual void connect_input(sprockit::sim_parameters* params,
								int src_outport,
								int dst_inport,
								event_handler* credit_handler) override;

	virtual link_handler* credit_handler(int port) const override;

	virtual link_handler* payload_handler(int port) const override;

	void receive_payload(event* ev); 

	virtual void execute(ami::COMP_FUNC func, event* data, callback* cb) override;

private:
	std::vector<event_handler*> outport_handler_;

	std::vector<event_handler*> inport_handler_;
	
	int id_;
	
	int num_groups_;
	
	int num_optical_switches_;
};

class optical_configuration_event : public flexfly_event {
public:
	optical_configuration_event();

	~optical_configuration_event() {};

	// copies the entries from one switch configuration vector into the other 
	void set_optical_configurations(std::vector<int> inout_connection) {
		inout_connection_ = inout_connection;
	};

private:
	std::vector<int> inout_connection_;
};

}
}