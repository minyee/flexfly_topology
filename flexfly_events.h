/**
 * These are just essentially classes that represent the packets in the Flexfly 
 * interconnect network
 */
#include <sstmac/hardware/common/connection.h>
#include <sstmac/common/sstmac_config.h>
#include <sstmac/common/sst_event.h>
typedef uint32_t byte_size;
namespace sstmac {


class flexfly_event : public event {
public:
ImplementSerializable(flexfly_event)
};

class flexfly_credit_event : public flexfly_event {

public:
	flexfly_credit_event(int src_id, int src_outport, int dst_id, int dst_inport) : flexfly_event() {
		src_id_ = src_id;
		src_outport_ = src_outport;
		dst_id_ = dst_id;
		dst_inport_ = dst_inport;
	};

	~flexfly_credit_event() {};
	//ImplementSerializable(flexfly_credit_event);

	int src_outport() const {
		return src_outport_;
	};

	int dst_inport() const {
		return dst_inport_;
	}

	int src_id() const {
		return src_id_;
	}

	int dst_id() const {
		return dst_id_;
	}

	void set_src_outport(int port) {
		src_outport_ = port;
	};

	void set_dst_inport(int port) {
		dst_inport_ = port;
	};

	void set_src_id(int id) {
		src_id_ = id;
	};

	void set_dst_id(int id) {
		dst_id_ = id;
	};

private:
	int src_id_;
	int src_outport_;
	int dst_inport_;
	int dst_id_;
	byte_size msg_size_;
	
};

class flexfly_payload_event : public flexfly_event {

public:
	flexfly_payload_event(int src_id, int src_outport, int dst_id, int dst_inport) : flexfly_event() {
		src_id_ = src_id;
		src_outport_ = src_outport;
		dst_id_ = dst_id;
		dst_inport_ = dst_inport;
	};

	~flexfly_payload_event() {};
	//ImplementSerializable(flexfly_credit_event);

	int src_outport() const {
		return src_outport_;
	};

	int dst_inport() const {
		return dst_inport_;
	}

	int src_id() const {
		return src_id_;
	}

	int dst_id() const {
		return dst_id_;
	}

	void set_src_outport(int port) {
		src_outport_ = port;
	};

	void set_dst_inport(int port) {
		dst_inport_ = port;
	};

	void set_src_id(int id) {
		src_id_ = id;
	};

	void set_dst_id(int id) {
		dst_id_ = id;
	};

private:
	int src_id_;
	int src_outport_;
	int dst_inport_;
	int dst_id_;
	byte_size msg_size_;
};


}