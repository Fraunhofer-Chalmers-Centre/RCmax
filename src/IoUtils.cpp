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

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

namespace rcmax {
	namespace io {
		void write_rcmax(std::string name, const std::vector<std::vector<int>>& c)
		{
			const auto m = c.size();
			const auto n = c.front().size();
			std::ofstream file;
			file.open(name, std::ofstream::out);
			for (size_t i = 0; i < m; ++i) {
				file << c[i][0];
				for (size_t j = 1; j < n; ++j) {
					file << "," << c[i][j];
				}
				file << std::endl;
			}
			file.close();
		}
		std::vector<std::vector<int>> read_rcmax(std::string name)
		{
			std::vector<std::vector<int>> c;
			std::ifstream infile(name);
			while (infile)
			{
				std::string line;
				if (!getline(infile, line)) break;
				c.emplace_back();
				std::istringstream ss(line);
				while (ss)
				{
					std::string s;
					if (!getline(ss, s, ',')) break;
					c.back().push_back(std::stoi(s));
				}
				if (c.back().empty()) {
					c.pop_back();
				}
			}

			return c;
		}

		std::default_random_engine generator;
		std::vector < std::vector<int>> generate_rcmax(int m, int n, int type, unsigned int seed)
		{
			generator.seed(seed);
			return generate_rcmax(m, n, type);

		}
		std::vector < std::vector<int>> generate_rcmax(int m, int n, int type)
		{
			std::vector<std::vector<int>> c(m, std::vector<int>(n, 0));
			if (type == 0) {
				std::uniform_int_distribution<int> distribution(10, 100);
				for (auto& ci : c) {
					for (auto& cij : ci) {
						cij = distribution(generator);
					}
				}
			}
			else if (type == 1) {
				std::uniform_int_distribution<int> p(1, 100);
				std::uniform_int_distribution<int> d(1, 20);
				for (int j = 0; j < n; ++j) {
					const auto cj = p(generator);
					for (int i = 0; i < m; ++i) {
						c[i][j] = cj + d(generator);
					}
				}
			}
			else {
				std::uniform_int_distribution<int> p(1, 100);
				std::uniform_int_distribution<int> d(1, 20);
				for (int i = 0; i < m; ++i) {
					const auto ci = p(generator);
					for (int j = 0; j < n; ++j) {
						c[i][j] = ci + d(generator);
					}
				}
			}
			return c;
		}
	} // namespace io
} // namespace rcmax