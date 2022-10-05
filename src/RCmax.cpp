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

#include "RCmax.h"

#include "LR1solver.h"
#include "RCmaxHeuristics.h"
#include "SubproblemSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
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

namespace { 
/**
* remove nodes that do not comply with the reduced costs and bounds, see paper for details
* @param c the processing times, sets, c[i][j] = ub+1 to denote a removed variable
* @param c_red the reduced costs, does not need to be additive
* @param c_red2 function to compute the reduced cost of two variables on the same machine
* @param lb is the current lower bound corresponding to the reduced  costs
* @param ub is the current upper bound 
* @return the number of variables removed
*/
template<class F>
int node_reduction(std::vector<std::vector<int>>& c, const std::vector<std::vector<double>>& c_red, const F& c_red2, const double lb, const int ub)
{
	int n_removed = 0;
	const auto m = c.size();
	const auto n = c.front().size();
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (c[i][j] < ub && round_up(lb + c_red[i][j]) >= ub) {
				c[i][j] = ub + 1;
				++n_removed;
			}
			for (size_t s = 0; s < j && c[i][j] <= ub && c_red[i][j] > 0.0; ++s) {
				if (c[i][s] >= ub || c_red[i][s] <= 0.0  || round_up(lb + c_red2(i, j, s)) < ub) {
					continue;
				}
				// both can not be one, see if s dominates j
				bool dominated = c[i][s] <= c[i][j];
				for (size_t r = 0; r < m && dominated; ++r) {
					dominated = dominated && (r == i || c[r][s] >= ub || c[r][s] >= c[r][j]);
				}
				if (dominated) {
					c[i][j] = ub + 1;
					++n_removed;
					continue;
				}
				// both can not be one, see if j dominates s
				dominated = c[i][s] >= c[i][j];
				for (size_t r = 0; r < m && dominated; ++r) {
					dominated = dominated && (r == i || c[r][j] >= ub ||  c[r][s] <= c[r][j]);
				}
				if (dominated) {
					c[i][s] = ub + 1;
					++n_removed;
				}
			}
		}
	}
	std::vector<size_t> fixed_j(n, m + 1);
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (c[i][j] < ub) {
				if (fixed_j[j] < m) {
					fixed_j[j] = m;
				}
				if (fixed_j[j] == m + 1) {
					fixed_j[j] = i;
				}
			}
		}
	}
	std::vector<int> z(m, 0);
	for (size_t j = 0; j < n; ++j) {
		if (auto i = fixed_j[j]; i < m) {
			z[i] += c[i][j];
		}
	}

	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (c[i][j] <= ub && z[i] + c[i][j] >= ub && i != fixed_j[j]) {
				c[i][j] = ub + 1;
				++n_removed;
			}
		}
	}
	return n_removed;
}

/**
* compute the step length in the dual space
* @param delta the scale factor
* @param lb the node lower bound
* @param ub the upper bound
* @param sg_norm2 the square norm of the subgradient vector 
*/
double step_length(const double delta, const double lb, const double ub, const double sg_norm2)
{
	return delta * (ub - lb) / sg_norm2;
}
/**
* perform the subgradient maximization
* @param sub_solver the knapsack solver for the R||Cmax subproblem
* @param sub_x the subproblem solution
* @param u the Lagrangian multiplier vector
* @param z_lb_frac the lower bound (not rounded upwards)
* @param z_ub the upper bound
* @param x the best primal solution found
* @param settings the feature to use and not use
* @param root if it is the root node
* @param x_ub pointer to the best known primal solution (output)
* @param pert_ofs an offset factor that pertubes the heuristic to settle ties in a deterministic (cyclic) fashion
* @return the number of iterations
*/
int sub_gradmax(knapsack::SubproblemSolver& sub_solver, std::vector<std::vector<bool>>& sub_x, std::vector<double>& u, double& z_lb_frac, int& z_ub, std::vector<int>& x, const BNB_SETTINGS& settings, const bool root, std::vector<int>* x_ub, double& pert_ofs)
{

	const auto& c = sub_solver.get_weights();
	const auto m = c.size();
	const auto n = c.front().size();
	const auto delta_min = 1e-3;
	auto iter_before_div = root ? 40 : 12;
	size_t n_iter = root ? 1000 : 100;
	if (!settings[idx(BNB_Features::SUB_GRAD_MANY_ITER)]) {
		n_iter = root ? 300 : 10;
		iter_before_div = root ? 30 : 3;
	}

	auto delta = root ? 0.5 : 1e-2;
	auto iter_since_last_div = 0;
	auto z_lb = round_up(z_lb_frac);

	std::vector<double> u_dir(n, 0.0), sg(n, 0.0);
	std::vector<std::vector<double>> x_seq(m, std::vector<double>(n, 0.0));
	double ergodic_weight = 0.0;
	auto best_u = u;
	size_t iter = 0;

	for (; iter < n_iter; ++iter_since_last_div) {
		++iter;
		auto z_lr_lb = sub_solver.solve(sub_x, nullptr, u, z_lb);
		
		if (z_lr_lb > z_lb_frac + 1e-6) {			
			z_lb_frac = z_lr_lb;
			best_u = u;
			iter_since_last_div = 0;
		}
		if (iter_since_last_div > iter_before_div) {
			iter_since_last_div = 0;
			delta *= .5;
		}

		z_lb = std::max(z_lb, round_up(z_lr_lb));

		if (z_lb >= z_ub) {
			break;
		}

		for (size_t j = 0; j < n; ++j) {
			sg[j] = 1.0;
			for (size_t i = 0; i < m; ++i) {
				if (sub_x[i][j]) {
					sg[j] -= 1.0;
				}
			}
		}

		const auto sg_norm2 = std::transform_reduce(sg.begin(), sg.end(), 0.0, std::plus<>(), [] (const double v) { return v * v; });
		
		if (iter == 1 || !settings[idx(BNB_Features::SUB_GRAD_DEFLECTED)]) {
			u_dir = sg;
		}
		else {
			sgupdates::MDS_update(u_dir, sg);
		}
		const auto tk = step_length(delta, z_lr_lb, z_ub, sg_norm2);
		std::transform(u.begin(), u.end(), u_dir.begin(), u.begin(), [tk] (const double u, const double d) {
			return std::max(u + tk * d, 0.0); });

		if (delta <= 1e-2 && settings[idx(BNB_Features::ERGODIC_SEQUENCE)]) {
			ergodic_weight += tk;
			for (size_t i = 0; i < m; ++i) {
				for (size_t j = 0; j < n; ++j) {
					x_seq[i][j] *= (ergodic_weight - tk) / ergodic_weight;
					if (sub_x[i][j]) {
						x_seq[i][j] += tk / ergodic_weight;
					}
				}
			}
		}

		if (root) {
			rcmax_heuristic(x, sub_x, c);
			if (iter % 10 == 0) {
				auto z_lr_ub = local_search_heuristic_cycle3<true>(x, c, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
				if (x_ub && z_lr_ub < z_ub) {
					*x_ub = x;
				}
				z_ub = std::min(z_ub, z_lr_ub);
			}
			else {
				auto z_lr_ub = local_search_heuristic_cycle3<false>(x, c, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
				if (x_ub && z_lr_ub < z_ub) {
					*x_ub = x;
				}
				z_ub = std::min(z_ub, z_lr_ub);
			}
		}
		sub_solver.set_ub(z_ub);

		if (z_lb >= z_ub || sg_norm2 <= .5 || delta_min > delta) {
			break;
		}
	}
	if (z_lb < z_ub) {
		if (!root) {
			rcmax_heuristic(x, sub_x, c);
			auto z_lr_ub = local_search_heuristic_cycle3<true>(x, c, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
			if (x_ub && z_lr_ub < z_ub) {
				*x_ub = x;
			}
			z_ub = std::min(z_ub, z_lr_ub);
		}
		if (settings[idx(BNB_Features::ERGODIC_SEQUENCE)]) {
			x = rcmax_heuristic_ergo(x_seq);
			auto z_lr_ub = local_search_heuristic_cycle3<true>(x, c, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
			if (x_ub && z_lr_ub < z_ub) {
				*x_ub = x;
			}
			z_ub = std::min(z_ub, z_lr_ub);
		}
	}
	u = best_u;
	return static_cast<int>(iter);
}
} // namespace
using NodeDef = std::pair< std::vector<std::vector<int>>, std::vector<double>>;
using QueueItem = std::pair<double, NodeDef>;
std::tuple<int, int, int, double, double> solve_bnb(const std::vector<std::vector<int>>& c, const BNB_SETTINGS settings, const double time_limit, std::vector<int>* x_ub)
{
	const bool verbose = !settings[idx(BNB_Features::NO_OUTPUT)];
	const auto m = c.size();
	const auto n = c.front().size();
	int print_cnt = 0;
	auto pert_ofs = 0.0;
	auto start = std::chrono::system_clock::now();
	LR1::LR1Solver lr1solver(m, n);
	auto z_lr1 = lr1solver.solve(c);

	auto x = rcmax_heuristic(c);
	auto z_ub = local_search_heuristic_cycle3<true>(x, c, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
	if (x_ub) {
		*x_ub = x;
	}

	x = lr1solver.primal_solution();
	auto v = local_search_heuristic_cycle3<true>(x, c, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
	if (x_ub && v < z_ub) {
		*x_ub = x;
	}
	z_ub = std::min(v, z_ub);

	auto z_lb_frac = z_lr1;

	knapsack::SubproblemSolver sub_solver(c, z_ub);
	knapsack::SubproblemSolverNoSol sub_solve_no_x(c, z_ub);

	std::vector<double> u = lr1solver.lp_task_duals();

	std::vector<std::vector<bool>> sub_x(m, std::vector<bool>(n));
	auto c_red = lr1solver.reduced_costs();

	std::priority_queue<QueueItem, std::vector<QueueItem>, std::less<QueueItem>> tree;
	tree.emplace(QueueItem{ z_lb_frac, { c, u } });
	bool root = true;
	auto z_root_gap = (static_cast<double>(z_ub) - round_up(z_lb_frac)) / z_ub;
	auto z_root_lb = round_up(z_lb_frac);
	auto n_nodes = 0;
	auto n_iter = 0;
	while (!tree.empty())
	{
		if (std::chrono::duration<double> time = std::chrono::system_clock::now() - start; time_limit < time.count()) {
			break;
		}
		++n_nodes;
		auto& node = tree.top();
		auto z_lb_node = node.first;
		auto c_node = node.second.first;
		u = node.second.second;
		if (root) {
			node_reduction(c_node, c_red, [&c_red](size_t i, size_t j, size_t s) {return c_red[i][j] + c_red[i][s]; }, z_lr1, z_ub);
			z_lr1 = lr1solver.solve(c_node); // exploit dual degeneracy
			c_red = lr1solver.reduced_costs();
			node_reduction(c_node, c_red, [&c_red](size_t i, size_t j, size_t s) {return c_red[i][j] + c_red[i][s]; }, z_lr1, z_ub);
		}
		sub_solver.update_weights(c_node);
		tree.pop();

		auto z_lb = round_up(z_lb_node);
		if (z_lb >= z_ub || z_lb_node >= z_ub) {
			continue;
		}
		if (root) {
			n_iter += sub_gradmax(sub_solver, sub_x, u, z_lb_node, z_ub, x, settings, root, x_ub, pert_ofs);
			node_reduction(c_node, c_red, [&c_red](size_t i, size_t j, size_t s) {return c_red[i][j] + c_red[i][s]; }, z_lr1, z_ub);
			z_root_lb = round_up(z_lb_node);
		}
		else {
			n_iter += sub_gradmax(sub_solver, sub_x, u, z_lb_node, z_ub, x, settings, root, x_ub, pert_ofs);
		}
		root = false;
		z_lb = round_up(z_lb_node);
		if (z_lb >= z_ub || z_lb_node >= z_ub) {
			continue;
		}
		c_node = sub_solver.get_weights();
		if (verbose && 10 * (++print_cnt) > std::min(n_nodes, 1000) + 9) {
			print_cnt = 0;
			std::chrono::duration<double> time = std::chrono::system_clock::now() - start;
			std::stringstream timestr;
			timestr << time.count() << "s";
			std::cout << "time: " << std::setw(10) << std::left << timestr.str() << "open nodes = " << std::setw(8) << std::left << tree.size() << "closed nodes = " << std::setw(8) << std::left << n_nodes << "node bounds: " << std::setw(10) << std::right << z_lb_node << " < " << std::setw(10) << std::left << z_ub << std::endl;
		}
		if (settings[idx(BNB_Features::VAR_FIXING)]) {
			auto z_lr_lb = sub_solver.solve(sub_x, &c_red, u, z_lb);
			if (round_up(z_lr_lb) >= z_ub) {
				continue;
			}
			auto c_red2 = [&sub_solver, &u, &z_lb](size_t i, size_t j, size_t s) { return sub_solver.reduced_cost(i, j, s, u, z_lb); };
			node_reduction(c_node, c_red, c_red2, z_lr_lb, z_ub);
			x = rcmax_heuristic(c_red);
			auto z_lr_ub = local_search_heuristic_cycle3<true>(x, c_node, settings[idx(BNB_Features::LARGE_LOCAL_SEARCH)], pert_ofs);
			if (z_lr_ub < z_ub && x_ub) {
				*x_ub = x;
			}
			z_ub = std::min(z_ub, z_lr_ub);
			if (z_lb >= z_ub || z_lb_node >= z_ub) {
				continue;
			}
		}
		else {
			z_lr1 = lr1solver.solve(c_node);
			if (round_up(z_lr1) >= z_ub) {
				continue;
			}
			c_red = lr1solver.reduced_costs();
			node_reduction(c_node, c_red, [&c_red](size_t i, size_t j, size_t s) {return c_red[i][j] + c_red[i][s]; }, z_lr1, z_ub);
		}

		std::vector<size_t> fixed_j(n, m + 1);
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				if (c_node[i][j] < z_ub) {
					if (fixed_j[j] < m) {
						fixed_j[j] = m;
					}
					else if (fixed_j[j] == m + 1) {
						fixed_j[j] = i;
					}
				}
			}
		}
		// check phantom criterion
		for (auto [j, n_to_check] = std::pair(size_t(0), n); n_to_check > 0; ++j, --n_to_check) {
			if (j == n) {
				j = 0;
			}
			for (size_t i = 0; i < m && fixed_j[j] >= m; ++i) {
				if (c_node[i][j] >= z_ub) {
					continue;
				}
				for (size_t s = 0; s < n; ++s) {
					if (size_t r = fixed_j[s]; r < m && r != i && c[r][s] >= c[r][j] && c[i][j] >= c[i][s]) {
						c_node[i][j] = z_ub + 1;
						n_to_check = n;
						for (size_t i2 = 0; i2 < m; ++i2) {
							if (c_node[i2][j] < z_ub) {
								if (fixed_j[j] < m) {
									fixed_j[j] = m;
									break;
								}
								fixed_j[j] = i2;
							}
						}
						break;
					}
				}
			}
		}
		// perform branching  
		auto jb = n;
		std::vector<std::pair<double, int>> cib;
		auto z_best_score = 0.0;
		sub_solve_no_x.set_ub(z_ub);
		for (size_t j = 0; j < n; ++j) {
			if (fixed_j[j] < m) {
				continue;
			}
			std::vector<std::pair<double, int>> ci;
			for (size_t i = 0; i < m; ++i) {
				if (c_node[i][j] < z_ub) {
					ci.emplace_back(c_red[i][j], static_cast<int>(i));
				}
			}
			if (ci.empty()) {
				jb = n;
				break;
			}
			else if (ci.size() == 1) {
				fixed_j[j] = ci.front().second;
				continue;
			}

			std::sort(ci.begin(), ci.end());
			const auto z_lr_score = m - ci.size() + 20.0 * ci.front().first;

			if (jb == n || (z_lr_score > z_best_score)) {
				jb = j;
				cib = ci;
				z_best_score = z_lr_score;
			}
		}
		if (n != jb) {
			auto c_left = c_node, c_right = c_node;
			bool to_left = true;
			for (auto [ctmp, i] : cib) {
				if (to_left) {
					to_left = false;
					c_right[i][jb] = z_ub + 1;
				}
				else {
					to_left = true;
					c_left[i][jb] = z_ub + 1;
				}
			}
			sub_solve_no_x.update_weights(c_left);
			auto z_lr_left = sub_solve_no_x.solve(u, z_lb);
			if (round_up(z_lr_left) < z_ub) {
				tree.emplace(std::pair<double, NodeDef>{ z_lr_left, { c_left, u } });
			}
			sub_solve_no_x.update_weights(c_right);
			auto z_lr_right = sub_solve_no_x.solve(u, z_lb);

			if (round_up(z_lr_right) < z_ub) {
				tree.emplace(std::pair<double, NodeDef>{ z_lr_right, { c_right, u } });
			}
		}
		if (n_nodes <= 1) {
			z_root_gap = (static_cast<double>(z_ub) - round_up(z_lb_node)) / z_ub;
		}
	}
	if (n_nodes <= 1) {
		z_root_gap = 0.0;
	}

	z_root_gap = std::max(z_root_gap, 0.0);
	auto z_lb_gap = (static_cast<double>(z_ub) - round_up(z_root_lb)) / z_ub;
	return std::tuple(n_nodes, n_iter, z_ub, z_root_gap, std::max(z_lb_gap, 0.0));
}
} // namespace rcmax