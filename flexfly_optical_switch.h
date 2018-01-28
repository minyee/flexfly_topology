/**
 * Author: Min Yee Teh
 *
 */
#include <sprockit/debug.h>                   
#include <sstmac/hardware/router/router_fwd.h>
#include <sstmac/hardware/topology/topology.h>
#include <sprockit/sim_parameters_fwd.h>
#include <sstmac/hardware/switch/network_switch.h>
//#include <sstmac/hardware/pisces/pisces.h>
#include "flexfly_topology.h"
#if SSTMAC_INTEGRATED_SST_CORE
#include <sstmac/sst_core/integrated_component.h>
#endif

#include <unordered_map>

//DeclareDebugSlot(network_switch)
#define switch_debug(...) \
  debug_printf(sprockit::dbg::network_switch, "Switch %d: %s", \
    int(addr()), sprockit::printf(__VA_ARGS__).c_str())


namespace sstmac {
namespace hw {

class optical_switch : public connectable_component {
  public:
  DeclareFactory(optical_switch,uint64_t,event_manager*)

  optical_switch(sprockit::sim_parameters* params, uint64_t id, event_manager* mgr) :
    connectable_component(params,id,
      device_id(params->get_int_param("id"),device_id::router),mgr) {
    };
};

/**
 * @brief The network_switch class
 * A class encapsulating a network switch that packets must traverse on the network.
 * The network switch performs both routing computations and congestion modeling.
 */
class flexfly_optical_switch :
  public optical_switch {

RegisterComponent("flexfly_optical_switch | flexfly_optical", optical_switch, flexfly_optical_switch,
        "macro", COMPONENT_CATEGORY_NETWORK,
        "An optical switch used in the Flexfly project");

 public:
  flexfly_optical_switch(sprockit::sim_parameters* params, uint64_t id, event_manager* mgr);

  virtual ~flexfly_optical_switch();

  virtual std::string to_string() const override {
    return "Flexfly Optical Switch";
  };

  virtual void deadlock_check() override;

  virtual void deadlock_check(event* ev) override;

  virtual void setup() override; //needed for SST core compatibility

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
   * @return
   */
  virtual link_handler* payload_handler(int port) const override;


  switch_id addr() const {
    return my_addr_;
  }; 

  bool outport_connected(int outport) const;


public:
  void recv_payload(event* ev);
  
  void recv_credit(event* ev);

  void recv_config_msg(event* ev){};

private:
  // sooner or later we need to simulate the fact that there will be configuration messages 
  

  bool setup_inout_connection(int inport, int outport);

private: 
  switch_id my_addr_;
  
  int num_ports_;
  
  std::unordered_map<int,event_handler*> outport_handler_;
  
  std::unordered_map<int,event_handler*> inport_handler_;
  // given an index, the value of the entry is the output port that said inport is currently connected to
  
  int* inout_connection_; 
  
  std::unordered_map<int, int> inport_connections_; // maps this switch's inport to the source switch's outport
  
  std::unordered_map<int, int> outport_connections_; // maps this switch's outport to the destination's inport
  
  void teardown_outport_connection(int outport);
  
  flexfly_topology* top_;

  int num_electrical_switches_;
};

}
}

