#include <iostream>

int flattened_index(int switch, int row, int col) {
	
}

void formulate_constraint_matrix(int num_groups) {
	int num_rows = num_groups;
	int num_cols = num_groups;
	for (int switch_id = 0; switch_id < num_groups - 1; switch_id++) {
		int switch_offset = switch_id * num_rows * num_cols;
		for (int row = 0; row < num_rows; row++) {
			int row_offset = switch_offset + row * num_cols;
			for (int col = 0; col < num_cols; col++ ) {
				int col_offset += row_offset + col;

			}
		}
	}

};

void formualte() {

};

void formu
