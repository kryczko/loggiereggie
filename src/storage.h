#ifndef _STORAGE_H_
#define _STORAGE_H_

#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
#include <limits>
#include <cmath>

class Column {
public:
	std::vector<double> data;
	std::vector<double> mean_norm_data;
	std::vector<double> norm_data;
	void add_data(double d) {
		this->data.push_back(d);
	}

	void print_data() {
		std::vector<double>& d = this->data;
		if (d.size() > 10) {
			std::cout << d[0] << "\n";
			std::cout << d[1] << "\n";
			std::cout << "." << "\n";
			std::cout << "." << "\n";
			std::cout << "." << "\n";
			std::cout << d[d.size() - 1] << "\n";
		} else {
			for (auto val : d) {
				std::cout << val << "\n";
			}
		}
	}

	void print_mnd_data() {
		std::vector<double>& d = this->mean_norm_data;
		if (d.size() > 10) {
			std::cout << d[0] << "\n";
			std::cout << d[1] << "\n";
			std::cout << "." << "\n";
			std::cout << "." << "\n";
			std::cout << "." << "\n";
		} else {
			for (auto val : d) {
				std::cout << val << "\n";
			}
		}
	}

	double max() {
		std::vector<double>& d = this->data;
		double max_val = std::numeric_limits<int>::min();
		for (auto& val : d ){
			if (val > max_val) {
				max_val = val;
			}
		}
		return max_val;
	}
	
	double min() {
		std::vector<double>& d = this->data;
		double max_val = std::numeric_limits<int>::max();
		for (auto& val : d ){
			if (val < max_val) {
				max_val = val;
			}
		}
		return max_val;
	}

	double mean() {
		std::vector<double>& d = this->data;
		double sum = 0;
		#pragma omp parallel for reduction(+:sum)
		for (int i = 0 ; i < d.size(); i ++) {
			sum = sum + d[i];
		}
		return sum / d.size();
	}

	void normalize() {
		this->mean_norm_data.resize(this->data.size());
		this->norm_data.resize(this->data.size());
		double avg = mean();
		double range = max() - min();
		#pragma omp parallel for
		for (int i = 0; i < this->data.size(); i ++) {
			this->mean_norm_data[i] = (this->data[i] - avg) / range;
			this->norm_data[i] = this->data[i] / range;
		}
	}
};

class Coeffs {
public:
	int n_coeffs;
	std::vector<double> coeffs;
	
	void init_coeffs(int n) {
		this->n_coeffs = n;
		for (int i = 0; i < n; i++) {
			this->coeffs.push_back(1.0);
		}

	}

	void print() {
		std::cout << "Coefficients:\n";
		for (auto& val : this->coeffs) {
			std::cout << val << "\t";
		}
		std::cout<< "\n";
	}
};
class Storage {
	std::string filename;
	int cat_column;
	

public:
	double alpha;
	int n_rows;
	int n_vars;

	std::vector<double> costs;
	std::vector<Column> columns;
	Coeffs coeffs;

	void set_filename(std::string name) {
		this->filename = name;
	}
	void set_alpha(double a) {
		this->alpha = a;
	}

	void checkfile() {
		std::vector<int> comma_count;
		std::ifstream datafile;
		datafile.open(filename.c_str());
		int count = 0;
		while (!datafile.eof()) {
			std::string stuff;
			std::getline(datafile, stuff);
			size_t n = std::count(stuff.begin(), stuff.end(), ',');
			comma_count.push_back(n+1);
			if (count > 1) {
				if (comma_count[count] != comma_count[count -1]) {
					std::cout << "Input file is corrupt, exiting...\n";
					exit(-1);
				}
			}
			count ++;
		}
		this->n_vars = comma_count[0];
		this->n_rows = count;
		this->cat_column = comma_count[0];
		this->columns.resize(n_vars);
	}

	void readfile() {
		std::ifstream datafile;
		datafile.open(filename.c_str());
		while(!datafile.eof()) {
			std::string stuff;
			std::getline(datafile, stuff);
			char* pch;
			pch = strtok (const_cast<char*>(stuff.c_str()), ",");
			int ccount = 0;
			while (pch != NULL) {
			    this->columns[ccount].add_data(std::stod(pch));
			    ccount ++;
			    if (ccount == this->n_vars) {
			    	ccount = 0;
			    }
			    pch = strtok (NULL, ",");
			  }
		}
	}

	double hypothesis(double x) {
		return 1.0 / ( 1.0 + exp(-x) );
	}


	double cost(std::vector<Column>& c, int row) {
		double dot = 0;
		for (int i = 1; i < this->coeffs.coeffs.size(); i ++) {
			dot += this->coeffs.coeffs[i] * c[i - 1].mean_norm_data[row];
		}
		dot += this->coeffs.coeffs[0];
		return hypothesis(dot);
	} 

	void update_coeffs() {
		for (int i = 1; i < this->coeffs.coeffs.size(); i ++) {
			double sum = 0;
			#pragma omp parallel for reduction(+:sum)
			for (int j = 0; j < this->n_rows; j ++) {
				std::vector<double>& data = this->columns[i - 1].mean_norm_data;
				double val = this->cost(this->columns, j);
				sum = sum + (val - this->columns[this->columns.size() - 1].data[j]) * data[j]; 
			}
			this->coeffs.coeffs[i] = this->coeffs.coeffs[i] - (this->alpha / this->n_rows) * sum; 
		}
	}
	std::vector<double> total_cost() {
		std::vector<double> a;
		for (int i = 0; i < this->n_rows; i ++) {
			double val = this->cost(this->columns, i);
			double total = -this->columns[this->columns.size() - 1].data[i] * log(val) -(1-this->columns[this->columns.size() - 1].data[i]) * log(1 - val);
			a.push_back(val);
		}
		return a;
	}

	double avg_cost() {
		double total_cost = 0;
		for (int i = 0; i < this->n_rows; i ++) {
			double val = this->cost(this->columns, i);
			double total = -this->columns[this->columns.size() - 1].data[i] * log(val) -(1-this->columns[this->columns.size() - 1].data[i]) * log(1 - val);
			total_cost += total;
		}
		return total_cost / this->n_rows;
	}

};




#endif