#include <iostream>
#include <fstream>

#include "storage.h"
#include "config.h"



int main(int argc, char* argv[]) {
	Storage storage;
	parse_yaml_file(storage);
	storage.checkfile();
	storage.readfile();

	for (auto& column : storage.columns) {
		column.normalize();
	}
	storage.coeffs.init_coeffs(storage.n_vars);

	std::ofstream output;
	output.open("cost.dat");
	for (int i = 0; i < 10000; i ++) {
		storage.coeffs.print();
		storage.update_coeffs();
		std::vector<double> cost = storage.total_cost();
		for (auto& elem : cost) {
			output << elem << "\t";
		}	
		output << "\n";
	}
	output.close();
	return 0;
}