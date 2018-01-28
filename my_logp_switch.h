#include <sstmac/common/event_handler.h>
#include <sstmac/common/event_scheduler.h>
#include <sstmac/common/messages/sst_message_fwd.h>
#include <sstmac/hardware/nic/nic_fwd.h>
#include <sstmac/hardware/switch/network_switch.h>
#include <sstmac/hardware/interconnect/interconnect_fwd.h>
#include <sprockit/unordered.h>
#include "flexfly_topology.h"
#include "flexfly_topology_simplified.h"

namespace sstmac {
namespace hw {

/**
 * @brief Implements a switch that does very basic congestion modeling
 *        using the LogGP model.  See "LogGP in Theory and Practice"
 *        by Hoefler and Schneider.
 */
class my_logp_switch :
  public network_switch
{
 public:
  RegisterComponent("my_logP | my_simple | my_LogP | my_logp", network_switch, my_logp_switch,
         "macro", COMPONENT_CATEGORY_NETWORK,
         "A switch that implements a basic delay model with no congestion modeling")

 public:

  my_logp_switch(sprockit::sim_parameters* params, uint64_t id, event_manager* mgr);

  void handle(event* ev);

  std::string to_string() const override {
    return "my_LogP switch";
  }

  int queue_length(int port) const override {
    return 0;
  }

  virtual ~my_logp_switch();

  void connect_output(sprockit::sim_parameters *params,
                      int src_outport, int dst_inport,
                      event_handler* handler) override;

  void connect_input(sprockit::sim_parameters *params,
                     int src_outport, int dst_inport,
                     event_handler* handler) override;

  link_handler* payload_handler(int port) const override;

  link_handler* credit_handler(int port) const override {
    return nullptr;
  }

 private:
  void incoming_message(message* msg, node_id src, node_id dst);
  void outgoing_message(message* msg, node_id src, node_id dst);
  void bcast_local_message(message* msg, node_id src);
  void forward_bcast_message(message* msg, node_id dst);

 private:
  double inj_bw_inverse_;

  timestamp inj_lat_;

  timestamp dbl_inj_lat_;

  double electrical_bw_;

  double optical_bw_;

  timestamp inv_electrical_bw_;

  timestamp inv_optical_bw_;

  interconnect* interconn_;

  int switches_per_group_;

  int nodes_per_switch_;

  std::vector<event_handler*> neighbors_;
  std::vector<event_handler*> nics_;

  flexfly_topology* ftop_;
  flexfly_topology_simplified* ftop_simplified_;

#if !SSTMAC_INTEGRATED_SST_CORE
  link_handler* mtl_handler_;
#endif

};

}
}

