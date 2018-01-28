#include <sstmac/hardware/router/routable.h>
#include <sstmac/hardware/pisces/pisces.h>
#include <queue>
#include <sprockit/factories/factory.h>
#include <sprockit/debug.h>

namespace sstmac{
namespace hw {

#ifndef FLEXFLY_PACKET
#define FLEXFLY_PACKET
/** Contains information on the switch target
  *
 **/
struct switch_port_pair {
	int switch_id;
	int outport;
};

struct flexfly_path {
	std::vector<switch_port_pair*> path;
	int path_length;
	int curr_index; 
};

class flexfly_packet : public pisces_default_packet {

//ImplementSerializable(flexfly_packet)

public:
	flexfly_packet(pisces_default_packet* pisces_packet_arg) {
		uint32_t default_packet_size = 1024;// in bytes
		set_packet_data_size(default_packet_size);
		pisces_packet_ = pisces_packet_arg;
		fpath_ = nullptr;
		counter_along_path_ = 0;
	};
	
	~flexfly_packet() {
		if (fpath_)
			delete fpath_;
	}; 

	/**
	* Returns the next outport in the switch. Usually called by switches to figure out which outport to use
	*
	**/

	int next_outport() {
		assert(fpath_);
		switch_port_pair* spp = fpath_->path[counter_along_path_];
		//std::cout << "in next_outport - counter_along_path is: " << std::to_string(counter_along_path_) <<std::endl;
		//std::cout << "size of path is: " << std::to_string(fpath_->path.size()) <<std::endl;
		assert(spp);
		//std::cout << "      switch id: " + std::to_string(spp->switch_id) + " outport: "  + std::to_string(spp->outport)<< std::endl;
		counter_along_path_++;
		if (counter_along_path_ >= fpath_->path_length) 
			counter_along_path_ = fpath_->path_length - 1;
		return spp->outport;
	}

	void reset_path_counter() {
		counter_along_path_ = 0;
	}
	/**
	 * pops the pointer to the next switch port pair. returns a nullptr if there are
	 * no more switches in the path. Note that this will decrement the count of the 
	 * switch__port_pair queue by 1.
	 **/
	switch_port_pair* next_switch() {
		if (!fpath_)
			return nullptr;
		switch_port_pair* to_ret = fpath_->path.front();
		//fpath_->path.pop();
		counter_along_path_++;
		return to_ret;
	};

	void set_path(flexfly_path* fpath) {
		fpath_ = fpath;
		return;
	};

	uint32_t get_packet_data_size() const {
		return packet_data_size_;	
	}

	void set_packet_data_size(uint32_t size) {
		packet_data_size_ = size;
		return;
	};

	pisces_default_packet *get_pisces_packet() const {
		return pisces_packet_;
	}

	virtual int next_vc() const override {
		return 0;
	};

private: 
	uint32_t packet_data_size_;
	flexfly_path* fpath_;
	uint8_t counter_along_path_;
	pisces_default_packet *pisces_packet_; 
	//route* r_;
};




#endif
}
}
