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
namespace sgupdates {
	/**
	 MDS_update deflects the search direction according to the Modified Deflected Subgradient method (MSD) (Rachid Belgacem and Abdessamad Amir, 2018)
	* @param dir the previous subgradient search direction to be updated
	* @param sub_grad the current subgradient
	*/
	void MDS_update(std::vector<double>& dir, const std::vector<double>& sub_grad);
} // namespace sgupdates
namespace LR1 {
	/**
	* The Euclidian projection of the vector u onto the unit simplex Σ u = 1, u ≥ 0, specialized from Nonlinear Programming: Theory and Algorithms by (Mokhtar S. Bazaraa et al., 2006)
	* @param u the vector to be projected
	*/
	void project_onto_unit_simplex(std::vector<double>& u);

	/**
	* Find the l = max(1e-3, argmin_l max_i ((1-l) zc_i  + l zn_i)) and set zc := (1-l) zc  + l zn, this to heuristically improve an upper bound corresponding to the master problem of the Dantzig-Wolfe decomposition
	*  @param z_conv is a convex combination of previous z vectors
	*  @param z_new is the current z vector, i.e., z_i = Σ_j p_ij x_ij
	*/
	void opt_conv_comb(std::vector<double>& z_conv, const std::vector<double>& z_new);
	/**
	*	Solves the LR1 relaxation from (Velde et al., 1991) naming comes from (Martello et al., 1997), i.e., relax the makespan constraints.
	*	Note that this is equivalent to the LP bound and it provides a feasible LP dual solution for computing residual costs
	*/
	class  LR1Solver
	{
	public:
		/**
		* Constructor
		* @param im is the number of machines
		* @param in is the number of tasks
		*/
		LR1Solver(size_t im, size_t in);

		/**
		* primal_solution, assumes solve has been called
		* @return the best primal (integer) solution found
		*/
		const auto& primal_solution() const { return x_sol; };
		/**
		* dual_solution, assumes solve has been called
		* @return the Lagragian multiplier vector
		*/
		const auto& dual_solution() const { return u; };
		/**
		* lp_task_duals, assumes solve has been called
		* @return the LP dual variables corresponding to the task constraints, i.e., Σ_i x_ij = 1
		*/
		const auto& lp_task_duals() const { return v; };
		/**
		* reduced_costs, assumes solve has been called
		*  @return the LP reduced costs, i.e, c[i][j]*u[i] - v[j]
		*/
		const auto& reduced_costs() const { return cred; };
		/**
		* maximize the Lagrangian dual function using the subgradient method
		* @param ci is the processing times, somtimes denoted p_ij
		* @return the best lower bound, which coincides with the value of LP relaxation, if the method converges (which is true for all problem in our test-bed)
		*/
		double solve(const std::vector<std::vector<int>>& ci);
	private:
		size_t m, n;
		// use transpose for better cache locality
		std::vector<std::vector<int>> cT;
		std::vector<std::vector<double>> cred;
		std::vector<double> u, v;
		std::vector<int> x_sol;
		bool root{ true };
	};
} //namespace LR1
}// namespace  rcmax