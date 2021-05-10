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
#include <bitset>
#include <tuple>
#include <vector>

namespace rcmax {
namespace knapsack {


	/**
	* Solve the subproblem retrieved by relaxing that each task should be done once
	* This version is to not compute the corresponding solution but only its value
	*/
	class SubproblemSolverNoSol
	{
	public:
		/**
		* constructor
		* @param c the processing times (knapsack weights)
		* @param i_ub the current upper bound
		*/
		SubproblemSolverNoSol(const std::vector<std::vector<int>>& c, const int i_ub);

		/**
		* update weights
		* @param c new weights (processing times)
		*/
		void update_weights(const std::vector<std::vector<int>>& c);
		
		/**
		* set a new upper bound
		* @param i_ub the new upper bound
		*/
		void set_ub(const int i_ub);
		
		/**
		* solve the Lagragian subproblem using dynamic programming
		* @param u the current Lagragian multipliers
		* @param lb the current lower bound
		* @return the computed lower bound
		*/
		double solve(const std::vector<double>& u, int lb);

	private:
		const size_t m, n;
		int ub;
		std::vector<std::vector<int>> weight;
		std::vector<std::vector<double>> value;
	};

	/**
	* Solve the subproblem retrieved by relaxing that each task should be done once
	*/
	class SubproblemSolver
	{
	public:
		/**
		* constructor
		* @param c the processing times (knapsack weights)
		* @param i_ub the current upper bound
		*/
		SubproblemSolver(const std::vector<std::vector<int>>& c, const int i_ub);

		/**
		* update weights
		* @param c new weights (processing times)
		*/
		void update_weights(const std::vector<std::vector<int>>& c);
		
		/**
		* set a new upper bound
		* @param i_ub the new upper bound
		*/
		void set_ub(const int i_ub);
	
		/**
		* solve the Lagragian subproblem using dynamic programming
		* @param x the subproblem solution
		* @param c_red_p pointer to the reduced costs, nullptr => do not compute
		* @param u the current Lagragian multipliers
		* @param lb the current lower bound
		* @return the computed lower bound
		*/
		double solve(std::vector<std::vector<bool>>& x, std::vector<std::vector<double>>* c_red_p, const std::vector<double>& u, int lb);
		/**
		* reduced costs for using variable x_ij and x_is
		* @param i machine
		* @param j task 1
		* @param s task 2
		* @param u Lagrangian multiplier
		* @param lb node lower bound
		* @return the reduced cost
		*/
		double reduced_cost(size_t i, size_t j, size_t s, const std::vector<double>& u, int lb) const;
		
		/**
		* get the weight matrix (processing times)
		* @return weights
		*/
		const auto& get_weights() const { return weight; };

	private:
		const size_t m, n;
		int ub;
		std::vector<double> last_opt_diff;
		std::vector<int> fixed_w;
		std::vector<size_t> fixed_j;
		std::vector<size_t> last_j;
		std::vector<std::vector<int>> weight;
		std::vector<std::vector<std::vector<double>>> value;
	};
}// namespace  knapsack
}; // namespace  rcmax