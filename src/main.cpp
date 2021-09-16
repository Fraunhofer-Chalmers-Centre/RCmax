/*
MIT License

Copyright (c) 2021 Fraunhofer-Chalmers Centre

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.																	*/

#include "IoUtils.h"
#include "RCmax.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
/*
	A command line interface to the branch-and-bound RCmax solver 
	use ? or help to see the usage or see readme.md
	Can be seen as example useage of the rcmax::bnb_solver function
*/
int main(int argc, char** argv)
{
	if (argc <= 1 || (argc == 2 && (argv[1][0] == '?' || std::string(argv[1]) == "help"))) {
		std::cout << "No instance file provided" << std::endl;
		std::cout << "Call is: RCmax.exe \"path_to_instance.dat\" time_limit_in_seconds flag1 flag2..." << std::endl;
		std::cout << "Potential flags are (sa. sensitivity analysis in paper):" << std::endl;
		std::cout << "\t-nV\t\tnot use variable fixing" << std::endl;
		std::cout << "\t-nS\t\tuse few subgradient iterations" << std::endl;
		std::cout << "\t-nD\t\tnot use deflected subgradient" << std::endl;
		std::cout << "\t-nL\t\tnot use large neighbourhood search" << std::endl;
		std::cout << "\t-nE\t\tnot use ergodic sequence heuristic" << std::endl;
		std::cout << "\t-o path_to_result_file\n\t\t\tdirect the result to the output_file" << std::endl;
		std::cout << "\t-v\t\tprint progress" << std::endl;
		std::cout << "\t-g m n t\tgenerate a new instance with m machines and n tasks and write to the instance path, where t = u,j,m, (unrelated, correlated jobs, correlated machines)" << std::endl;
		return 1;
	}
	auto get_option = [argv, argc] (std::string opt) {return std::find_if(argv, argv + argc, [opt] (char* cmd) { return std::string(cmd) == opt; }); };
	auto has_option = [argv, argc] (std::string opt) {return std::find_if(argv, argv + argc, [opt] (char* cmd) { return std::string(cmd) == opt; }) != argv + argc; };

	// Determine output
	std::ofstream out_file;

	auto write_header = true;
	auto& out = [argv, argc, &out_file, get_option,&write_header] () ->std::ostream& {
		auto it = get_option("-o");
		if (it != argv + argc && it + 1 != argv + argc) {
			auto out_file_name = *(it + 1);
			if (std::filesystem::exists(out_file_name)) {
				write_header = false;
			}
			out_file.open(out_file_name, std::ofstream::app);
			return out_file;
		}
		return std::cout;
	}();

	// read or write?
	auto instance_file = argv[1];
	const bool generate = has_option("-g");
	if (generate) {
		int mg, ng;
		try {
			mg = std::stoi(*(get_option("-g") + 1));
			ng = std::stoi(*(get_option("-g") + 2));
		}
		catch(const std::exception&){
			std::cout << "Error: cannot read integers m = " << *(get_option("-g") + 1) << ", n = " << *(get_option("-g") + 2);
			return 1;
		};
		auto tg = (*(get_option("-g") + 3))[0];
		if (mg <= 0 || mg > 20 || ng <= 0 || ng > 200) {
			std::cout << "Error: Too large problem, this algorithm is only tested for m <= 20, n <= 200, m = " << mg << ", n = " << ng << std::endl;
			return 1;
		}
		auto type = -1;
		if (tg == 'u') {
			type = 0;
		}
		if (tg == 'j') {
			type = 1;
		}
		if (tg == 'm') {
			type = 2;
		}
		if (type < 0) {
			std::cout << "Error: type = " << tg << " is not supported" << std::endl;
			return 1;
		}
		auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
		rcmax::io::write_rcmax(instance_file, rcmax::io::generate_rcmax(mg, ng, type, seed));
	}

	const auto c = rcmax::io::read_rcmax(instance_file);
	if (c.empty()) {
		std::cout << "Cannot read file: " << instance_file << std::endl;
		return 1;
	}

	auto m = c.size();
	auto n = c.front().size();
	if (!std::all_of(c.begin(), c.end(), [n] (const auto& v) {return v.size() == n; })) {
		std::cout << "Error: Corrput dimensions: not same number of task for each machine. "<< std::endl;
		return 1;
	}
	//if (m > 20 || n > 200) {
	//	std::cout << "Error: Too large problem, this algorithm is only tested for m <= 20, n <= 200, m = " << m << ", n = " << n;
	//	return 1;
	//}
	if (n == 0) {
		std::cout << "Error: No tasks found";
		return 1;
	}
	//if (!std::all_of(c.begin(), c.end(), [](const auto& v) {return std::all_of(v.begin(), v.end(), [](auto c) {return (c >= 1 && c <= 200) || c == 5000; }); }))
	//{
	//	std::cout << "Error: Invalid processing times, assume 1 <= p_ij <= 200, use p_ij = 5000, to indicate an unavailable machine for task j." << std::endl;
	//	return 1;
	//}
	//for (size_t j = 0; j < n; ++j) {
	//	if (!std::any_of(c.begin(), c.end(), [j] (const auto& v) {return v[j] <= 200; })) {
	//		std::cout << "Error: Task " << j + 1 << "has no available machine.";
	//		return 1;
	//	}
	//}
	auto time_limit = 120.0;
	if (argc >= 3) {
		try {
			time_limit = std::stod(argv[2]);
		}catch (const std::exception&) {
			std::cout << "Error: cannot interpret as a time limit " << argv[2] << std::endl;
			return 1;
		};
	}
	
	if (!std::isfinite(time_limit) || time_limit < 0.0) {
		std::cout << "Invalid time_limit!" << time_limit <<  std::endl;
		return 1;
	}
	rcmax::BNB_SETTINGS settings = rcmax::remove_settings<>();
	if (has_option("-nV")) {
		settings = rcmax::remove_settings<rcmax::BNB_Features::VAR_FIXING>(settings);
	}
	if (has_option("-nS")) {
		settings = rcmax::remove_settings<rcmax::BNB_Features::SUB_GRAD_MANY_ITER>(settings);
	}
	if (has_option("-nD")) {
		settings = rcmax::remove_settings<rcmax::BNB_Features::SUB_GRAD_DEFLECTED>(settings);
	}
	if (has_option("-nL")) {
		settings = rcmax::remove_settings<rcmax::BNB_Features::LARGE_LOCAL_SEARCH>(settings);
	}
	if (has_option("-nE")) {
		settings = rcmax::remove_settings<rcmax::BNB_Features::ERGODIC_SEQUENCE>(settings);
	}
	if (has_option("-v")) {
		settings = rcmax::remove_settings<rcmax::BNB_Features::NO_OUTPUT>(settings);
	}

	std::vector<int> x;
	auto start = std::chrono::system_clock::now();
	auto [nodes, iter, z, gap_root, lb_gap_root] = rcmax::solve_bnb(c, settings, time_limit, &x);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> time = end - start;
	auto solved = time.count() + 1e-6 < time_limit;
	if (!has_option("-o")) {
		if (write_header) {
			out << "m   n   nodes   iter     time      solved gap_root% z     x"<<std::endl;
		}
		out << std::setw(4) << std::left << m << std::setw(4) << std::left << n << std::setw(8) << std::left << nodes << std::setw(9) << std::left << iter  <<
			std::setw(10) << std::left << time.count() << std::setw(7) << std::left << solved << std::setw(10) << std::left << 100.0*gap_root << std::setw(6) << std::left << z;
		for (auto xv : x) {
			out << xv +1 << " ";
		}
		out << std::endl;
	}
	else {
		if (write_header) {
			out << "m,n,nodes,iter,time,solved,gaproot,lbgaproot,z,x" << std::endl;
		}
		out << m << "," << n << "," << nodes << "," << iter << "," << time.count() << "," << solved << "," << gap_root << "," << lb_gap_root << "," << z << ",";
		for (auto xv : x) {
			out << xv + 1 << ",";
		}
		out << std::endl;
	}
	return 0;
}