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

#pragma once
#include <string>
#include <vector>

namespace rcmax {
	namespace io {
		/**
		* generates an instance of RCmax, using a specific random seed
		* @param m the number of machines
		* @param n the number of tasks
		* @param type the type of random distribution, 
		*			0 => uniform p_ij ~ (10,100)
		*			1 => correlated jobs p_ij = a_j + b_ij, a_j ~ (1,100) b_ij ~ (1,20)
		*			2 => correlated machine p_ij = a_i + b_ij, a_i ~ (1,100) b_ij ~ (1,20)
		* @param seed the random seed to use
		* @return the processing times
		*/
		std::vector<std::vector<int>> generate_rcmax(int m, int n, int type, unsigned int seed);
		
		/**
		* generates an instance of RCmax, continuing from previous seed
		* @param m the number of machines
		* @param n the number of tasks
		* @param type the type of random distribution, 
		*			0 => uniform p_ij ~ (10,100)
		*			1 => correlated jobs p_ij = a_j + b_ij, a_j ~ (1,100) b_ij ~ (1,20)
		*			2 => correlated machine p_ij = a_i + b_ij, a_i ~ (1,100) b_ij ~ (1,20)
		* @return the processing times
		*/
		std::vector<std::vector<int>> generate_rcmax(int m, int n, int type = 0);

		/**
		* Write the RCmax instance to a file
		* @param name the name of the file
		* @param p the processing times
		*/
		void write_rcmax(std::string name, const std::vector<std::vector<int>>& p);
		/**
		* Read an RCmax instance from a file
		* @param name the name of the file
		* @return the processing times
		*/
		std::vector<std::vector<int>> read_rcmax(std::string name);
	} // namespace  io
}; // namespace  rcmax