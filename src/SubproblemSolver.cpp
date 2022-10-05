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

#include "SubproblemSolver.h"

#include <algorithm>
#include <numeric>
namespace rcmax {
	namespace {
		/**
		* ceil to the nearest integer, the point of this function is to ceil in a conservative way, so that round_up(int + eps) = int, to avoid problems with finite precision
		*	@param v the fractional number
		*	@return ceil(v)
		*/
		int round_up(const double v)
		{
			return static_cast<int>(v + .999999);
		}
	}

namespace knapsack {
	SubproblemSolverNoSol::SubproblemSolverNoSol(const std::vector<std::vector<int>>& c, const int i_ub) :
		m(c.size()),
		n(c.front().size()),
		ub(i_ub),
		weight(c),
		value(m, std::vector<double>(static_cast<size_t>(i_ub), 0.0))
	{}

	void SubproblemSolverNoSol::update_weights(const std::vector<std::vector<int>>& c)
	{
		weight = c;
	}
	void SubproblemSolverNoSol::set_ub(const int i_ub)
	{
		if (i_ub < ub) {
			ub = i_ub;
			for (auto& vv : value) {
				vv.resize(static_cast<size_t>(ub));
			}
		}
	}
	double SubproblemSolverNoSol::solve(const std::vector<double>& u, int lb)
	{
		std::vector<double> knap_sum(static_cast<size_t>(ub) - lb, 0.0);
		std::vector<size_t> fixed_j(n, m);
		std::vector<int> fixed_w(m, 0);
		std::vector<double> fixed_v(m, 0.0);

		for (size_t j = 0; j < n; ++j) {
			for (size_t i = 0; i < m; ++i) {
				if (weight[i][j] < ub) {
					if (fixed_j[j] < m) {
						fixed_j[j] = m;
						break;
					}
					fixed_j[j] = i;
				}
			}
			if (auto i = fixed_j[j]; i < m) {
				fixed_w[i] += weight[i][j];
				fixed_v[i] += u[j];
			}
		}

		if (*std::max_element(fixed_w.begin(), fixed_w.end()) > lb) {
			lb = *std::max_element(fixed_w.begin(), fixed_w.end());
		}
		if (lb >= ub) {
			return static_cast<double>(lb);
		}
		for (size_t i = 0; i < m; ++i) {
			auto& vv = value[i];
			auto ub_i = ub - fixed_w[i];
			std::fill(vv.begin(), vv.end(), 0.0);
			for (size_t j = 0; j < n; ++j) {
				if (weight[i][j] >= ub_i || i == fixed_j[j]) {
					continue;
				}
				const auto v = u[j];
				std::transform(vv.begin() + weight[i][j], vv.begin() + ub_i, vv.begin(), vv.begin(), [v](double v1, double v2) { return std::max(v1 + v, v2); });
			}
			for (int w = lb; w < ub; ++w) {
				knap_sum[static_cast<size_t>(w) - lb] += vv[ub - static_cast<size_t>(w) - 1] + fixed_v[i];
			}
		}
		const auto u_sum = std::reduce(u.begin(), u.end());

		auto opt_lb = static_cast<double>(ub);
		auto opt_w = ub;
		for (int w = lb; w < ub; ++w) {
			const auto z_lb = u_sum + w - knap_sum[static_cast<size_t>(w) - lb];
			const auto z_lb_i = round_up(z_lb);
			if (z_lb_i <= w && z_lb < opt_lb) {
				opt_lb = z_lb;
				opt_w = w;
			}
		}
		return opt_lb;
	}

	SubproblemSolver::SubproblemSolver(const std::vector<std::vector<int>>& c, const int i_ub) :
		m(c.size()),
		n(c.front().size()),
		ub(i_ub),
		fixed_w(c.size()),
		fixed_j(c.front().size()),
		last_j(c.front().size()),
		weight(c),
		value(m, std::vector<std::vector<double>>(n, std::vector<double>(static_cast<size_t>(i_ub), 0.0)))
	{}

	void SubproblemSolver::update_weights(const std::vector<std::vector<int>>& c)
	{
		weight = c;
		for (size_t i = 0; i < m; ++i) {
			std::fill(value[i][0].begin(), value[i][0].end(), 0.0);
		}
	}

	void SubproblemSolver::set_ub(const int i_ub)
	{
		if (i_ub < ub) {
			ub = i_ub;
			for (auto& vm : value) {
				for (auto& vmn : vm) {
					vmn.resize(static_cast<size_t>(ub));
				}
			}
		}
	}

	/**
	* solve the Lagrangian subproblem using dynamic programming
	* @param x the subproblem solution
	* @param c_red_p pointer to the reduced costs, nullptr => do not compute
	* @param u the current Lagrangian multipliers
	* @param lb the current lower bound
	* @return the computed lower bound
	*/
	double SubproblemSolver::solve(std::vector<std::vector<bool>>& x, std::vector<std::vector<double>>* c_red_p, const std::vector<double>& u, int lb)
	{
		std::vector<double> knap_sum(static_cast<size_t>(ub) - lb, 0.0);
		last_opt_diff = knap_sum;
		last_j.assign(m, 0);
		fixed_j.assign(n, m);
		fixed_w.assign(m, 0);
		std::vector<double> fixed_v(m, 0.0);
		for (size_t j = 0; j < n; ++j) {
			for (size_t i = 0; i < m; ++i) {
				x[i][j] = false;
				if (weight[i][j] < ub) {
					if (fixed_j[j] < m) {
						fixed_j[j] = m;
						break;
					}
					fixed_j[j] = i;
				}
			}
			if (auto i = fixed_j[j]; i < m) {
				fixed_w[i] += weight[i][j];
				fixed_v[i] += u[j];
				x[i][j] = true;
			}
		}
	
		if (*std::max_element(fixed_w.begin(), fixed_w.end()) > lb) {
			lb = *std::max_element(fixed_w.begin(), fixed_w.end());
		}
		if (lb >= ub) {
			return static_cast<double>(lb);
		}
		for (size_t i = 0; i < m; ++i) {	
			auto& vm = value[i];
			auto ub_i = ub - fixed_w[i];
			for (int w = weight[i][0]; w < ub_i; ++w) {
				vm[0][w] = u[0];
			}
			if (i == fixed_j[0]) {
				for (int w = weight[i][0]; w < ub_i; ++w) {
					vm[0][w] = 0.0;
				}
			}
			for (size_t j = 1; j < n; ++j) {
				if (weight[i][j] >= ub_i || i == fixed_j[j]) {
					continue;
				}
				std::copy(vm[last_j[i]].begin(), vm[last_j[i]].begin() + weight[i][j], vm[j].begin()); const auto v = u[j];
				std::transform(vm[last_j[i]].begin() + weight[i][j], vm[last_j[i]].begin() + ub_i, vm[last_j[i]].begin(), vm[j].begin() + weight[i][j],
					[&v](const double v1, const double v2) {return std::max(v1, v2 + v); });
				last_j[i] = j;
			}
			for (int w = lb; w < ub; ++w) {
				knap_sum[static_cast<size_t>(w) - lb] += vm[last_j[i]][static_cast<size_t>(w) - fixed_w[i]] + fixed_v[i];
			}
		}
		const auto u_sum = std::reduce(u.begin(), u.end());

		auto opt_lb = static_cast<double>(ub);
		auto opt_w = ub;
		for (int w = lb; w < ub; ++w) {
			const auto z_lb = u_sum + w - knap_sum[static_cast<size_t>(w) - lb];
			const auto z_lb_i = round_up(z_lb);
			if (z_lb_i <= w && z_lb < opt_lb) {
				opt_lb = z_lb;
				opt_w = w;
			}
		}
		if (opt_lb == ub) {
			return opt_lb;
		}
		for (int w = lb; w < ub; ++w) {
			last_opt_diff[static_cast<size_t>(w) - lb] = std::max(0.0, u_sum + w - knap_sum[static_cast<size_t>(w) - lb] - opt_lb);
		}

		for (size_t i = 0; i < m; ++i) {
			const auto& vm = value[i];
			auto wc = opt_w - fixed_w[i];
			auto ub_i = ub - fixed_w[i];
			size_t j = n - 1;
			for (; j > 0 && (weight[i][j] >= ub_i || i == fixed_j[j]); --j) {
				x[i][j] = i == fixed_j[j];				
			}
			for (size_t jp = j; jp > 0;) {
				--jp;
				if (weight[i][jp] >= ub_i || i == fixed_j[jp]) {
					x[i][jp] = i == fixed_j[jp];
					continue;
				}
				x[i][j] = vm[j][wc] != vm[jp][wc];
				if (x[i][j]) {
					wc -= weight[i][j];
				}
				j = jp;
			}
			x[i][j] = vm[j][wc] > 0.0 || i == fixed_j[j];
		}
		// early cheep reduction
		if (c_red_p == nullptr) {
			for (size_t i = 0; i < m; ++i) {
				auto& vm = value[i][last_j[i]];
				for (size_t j = 0; j < n; ++j) {
					if (weight[i][j] >= ub || fixed_j[j] < m || x[i][j]) {
						continue;
					}
					bool viol = true;
					for (auto w = std::max(lb - fixed_w[i], weight[i][j]); w < ub - fixed_w[i]; ++w) {
						if (viol) {
							auto c_red_z = vm[w] - vm[static_cast<size_t>(w) - weight[i][j]] - u[j] + last_opt_diff[static_cast<size_t>(w) - lb + fixed_w[i]];
							viol = round_up(c_red_z + opt_lb) >= ub;
						}
					}
					if (viol) {
						weight[i][j] = ub + 1;
					}
				}
			}
		}

		if (c_red_p != nullptr) {
			auto& c_red = *c_red_p;
			//	Bidirectional improvement
			std::vector<std::vector<double>> vvs(m), W(m);
			for (size_t i = 0; i < m; ++i) {
				auto ub_i = ub - fixed_w[i];
				vvs[i].assign(ub_i, 0.0);
				W[i].assign(ub_i, 0.0);
			}
			for (int j = static_cast<int>(n) - 1; j >= 0; --j) {
				// Compute W
				for (size_t i = 0; i < m; ++i) {
					auto ub_i = ub - fixed_w[i];
					if (i == fixed_j[j] || weight[i][j] >= ub_i) {
						continue;
					}
					const auto& vv = vvs[i];
					int jn = j - 1;
					while (jn >= 0 && (weight[i][jn] >= ub_i || i == fixed_j[jn]))
					{
						--jn;
					}
					if (jn < 0) {
						for (int w = std::max(lb - fixed_w[i] - weight[i][j], 0); w < ub_i; ++w) {
							W[i][w] = vv[ub_i - w - 1];
						}
					}
					else {
						const auto& vvn = value[i][jn];
						for (int w = std::max(lb - fixed_w[i] - weight[i][j], 0); w < ub_i - weight[i][j]; ++w) {
							W[i][w] = std::transform_reduce(vv.begin() + ub_i - w - 1, vv.end(), vvn.begin(), 0.0,
								[](double v1, double v2) {return std::max(v1, v2); },
								std::plus());
						}
						for (int w = lb - fixed_w[i]; w < ub_i; ++w) {
							W[i][w] = std::transform_reduce(vv.begin() + ub_i - w - 1, vv.end(), vvn.begin(), 0.0,
								[](double v1, double v2) {return std::max(v1, v2); },
								std::plus());
						}
					}
				}

				for (size_t i = 0; i < m; ++i) {
					auto ub_i = ub - fixed_w[i];

					c_red[i][j] = std::numeric_limits<double>::infinity();
					if (i == fixed_j[j]) {
						c_red[i][j] = 0.0;
						continue;
					}
					if (weight[i][j] >= ub_i) {
						continue;
					}

					auto opt_lb_ij = static_cast<double>(ub);
					for (int w = std::max(lb, fixed_w[i] + weight[i][j]); w < ub; ++w) {
						auto lb_ij_w = u_sum + w - u[j];
						for (size_t ii = 0; ii < m; ++ii) {
							if (ii != i) {
								if (weight[ii][j] >= ub - fixed_w[ii]) {
									lb_ij_w -= value[ii][last_j[ii]][w - fixed_w[ii]] + fixed_v[ii];
								}
								else {
									lb_ij_w -= W[ii][w - fixed_w[ii]] + fixed_v[ii];
								}

							}
						}
						lb_ij_w -= W[i][w - fixed_w[i] - weight[i][j]] + fixed_v[i];
						const auto lb_ij_w_i = round_up(lb_ij_w);
						if (lb_ij_w_i <= w && lb_ij_w < opt_lb_ij) {
							opt_lb_ij = lb_ij_w;
						}
					}
					opt_lb_ij = std::max(opt_lb_ij, static_cast<double>(fixed_w[i]) + weight[i][j] - 1e-5);
					c_red[i][j] = opt_lb_ij - opt_lb;

					// add item 
					auto& vv = vvs[i];
					const auto v = u[j];
					std::transform(vv.begin() + weight[i][j], vv.begin() + ub_i, vv.begin(), vv.begin(), [v](double v1, double v2) { return std::max(v1 + v, v2); });
				}
			}
		}
		return opt_lb;
	}

	double SubproblemSolver::reduced_cost(size_t i, size_t j, size_t s, const std::vector<double>& u, int lb) const {
	
		auto& vm = value[i][last_j[i]];
		auto c_red = std::numeric_limits<double>::infinity();
		for (auto w = std::max(lb - fixed_w[i], weight[i][j] + weight[i][s]); w < ub - fixed_w[i]; ++w) {
			c_red = std::min(c_red, std::max(vm[w] - vm[static_cast<size_t>(w) - weight[i][j] - weight[i][s]] - u[j] - u[s] + last_opt_diff[static_cast<size_t>(w) - lb + fixed_w[i]], 0.0));
		}
		return c_red;
	}

} // namespace knapsack
} // namespace rcmax