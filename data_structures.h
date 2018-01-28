#include <vector>
#include <iostream>
#include <cassert>

#ifndef DATA_STRUCTURES_MINE
#define DATA_STRUCTURES_MINE

class node {
public:
	node(int id) : id_(id), num_child_(0) {
		//set_null_children();
		id_ = id;
		num_child_ = 0;
	}

	node(int id, size_t num_child, bool is_optical_switch) : id_(id), num_child_(num_child) {
		children_.reserve(num_child_);
		children_.resize(num_child_);
		id_ = id;
		num_child_ = num_child;
		is_optical_switch_ = is_optical_switch;
		set_null_children();
	};

	~node() {};

	int get_num_child() const {
		int return_val = num_child_;
		assert(num_child_==children_.size());
		return return_val;
	};

	void set_num_child(int num_child) {
		if (num_child > 0) {
			num_child_ = num_child;
			children_.reserve(num_child_);
			children_.resize(num_child_);
		}
	};

	void set_child(int index, node* child_node) {
 		assert(index >= 0 && index < num_child_);
 		//assert(child_node);
		if (index < num_child_ && index >= 0) {
			children_[index] = child_node;
		}
	};

	node* get_child(int index) const {

		node* child = nullptr;
		assert(index < num_child_ && index >= 0);
		if (index < num_child_ && index >= 0) {
			child = children_[index];
		}
		return children_[index];
		//return child;
	};

	int get_id() const {
		return id_;
	};

	bool is_optical_switch() const {
		return is_optical_switch_;
	}

private:
	void set_null_children() {
		if (num_child_ == 0) {
			return;
		}
		for (int i = 0; i < num_child_; i++) {
			children_[i] = nullptr;
		}
	};

private:
	int num_child_;
	
	int id_;
	
	std::vector<node*> children_;

	bool is_optical_switch_;
};
#endif
