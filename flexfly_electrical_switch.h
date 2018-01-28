#include <sstmac/hardware/switch/network_switch.h>
#include <sstmac/hardware/topology/topology.h>
#include <vector> 
#include <sstmac/hardware/router/router_fwd.h>
#include <sstmac/hardware/router/router.h>
#include <sstmac/hardware/router/minimal_routing.h>
#include "flexfly_topology.h"
#include "flexfly_topology_simplified.h"

//RegisterNamespaces("switch");
namespace sstmac {
namespace hw {

//typedef int port_id;

class flexfly_electrical_switch : public network_switch {
public:

RegisterComponent("flexfly_electrical_switch | flexfly_electrical", network_switch, flexfly_electrical_switch,
         			"macro", COMPONENT_CATEGORY_NETWORK,
         			"An electrical switch which has electrical components to it which is used in the Flexfly project");
 
flexfly_electrical_switch(sprockit::sim_parameters* params,
							uint64_t id,
							event_manager* mgr,
							device_id::type_t ty = device_id::router);// : network_switch(params, id, mgr, ty);
 
~flexfly_electrical_switch();

virtual std::string to_string() const override {
	return "Flexfly Electrical Switch for the Flexfly project";
};

/**
 * @brief connect_output
 * @param params
 * @param src_outport
 * @param dst_inport
 * @param payload_handler
 */
virtual void connect_output(sprockit::sim_parameters* params, 
                                int src_outport, 
                                int dst_inport, 
                                event_handler* payload_handler) override;


/**
 * @brief connect_input
 * @param params
 * @param src_outport
 * @param dst_inport
 * @param credit_handler Can be null, if no credits are ever sent
 */
virtual void connect_input(sprockit::sim_parameters* params, 
                              int src_outport, 
                              int dst_inport,
                              event_handler* credit_handler) override;

/**
 * @brief credit_handler
 * @param port
 * @return Can be null, if no credits are ever to be received
 */
virtual link_handler* credit_handler(int port) const override;

/**
 * @brief payload_handler
 * @param port
 *  @return
 */
virtual link_handler* payload_handler(int port) const override;

virtual int queue_length(int port) const override;

protected:

void recv_payload(event* ev);

void recv_credit(event* ev);

void recv_nodal_payload(event* ev);

void recv_nodal_credit(event* ev);

private:
void send_packet_to_node(int node_id, timestamp delay, event* ev);

void send_credit_to_node(uint32_t credit_amount, int node_id);

void get_link_params(sprockit::sim_parameters* switch_param);

protected:
router* router_;
std::vector<event_handler*> inport_handlers_;
std::vector<event_handler*> outport_handlers_;

private:
//uint64_t my_id_;
int radix_;
int *queue_length_;
switch_id my_addr_;
int port_cnt_;
int switches_per_group_;
int nodes_per_switch_;
flexfly_topology* ftop_;
flexfly_topology_simplified* ftop_simplified_;
int num_groups_;

 //void (std::vector<>)
/*
 * IMPORTANT: ADDED (01/17/2018 by minyee)
 * BEGIN
 * Note: These are all the bandwidths and send and credit lantecies 
 *
 */
double optical_link_bw_;
double electrical_link_bw_;
timestamp inv_electrical_link_bw_;
timestamp inv_optical_link_bw_;
timestamp send_latency_;
timestamp credit_latency_;
/*
 * IMPORTANT: ADDED (01/17/2018 by minyee)
 * END
 * Note: These are all the bandwidths and send and credit lantecies 
 *
 */
};
 
}
}