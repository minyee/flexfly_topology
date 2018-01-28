#include "data_structures.h"
#include <vector>
#include <unordered_map>

void generate_butterfly(int num_groups,  std::vector<node*>& groups, std::vector<node*>& optical_switches);

void link_butterfly(int num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches);

void form_butterfly(int num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches);

void delete_butterfly(int num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches);

void canonical_dragonfly_config_greedy(int num_groups, std::vector<node*>& groups, 
										std::vector<node*>& optical_switches,  
										std::vector<std::vector<int> >& optical_inout_connections);

void configure_optical_switches_canonical(int num_groups, std::vector< std::vector<int> >& optical_inout_connections);

node* dfs(node* curr_node, int depth, int target_group_id);

void configure_default_simpler_model(int num_groups, std::unordered_map<int, std::vector<int>>& optical_inout_connections);
