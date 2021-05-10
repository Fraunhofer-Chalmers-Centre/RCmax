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
// The features corresponding to the sensitivity analysis in the paper
enum class BNB_Features
{
	SUB_GRAD_DEFLECTED,
	SUB_GRAD_MANY_ITER,
	VAR_FIXING,
	LARGE_LOCAL_SEARCH,
	ERGODIC_SEQUENCE,
	NO_OUTPUT,
	Count
};
// determines the index of a feature
constexpr size_t idx(const BNB_Features feat) { return static_cast<size_t>(feat); }
using BNB_SETTINGS = std::bitset<idx(BNB_Features::Count)>;
// The default is to use all features, this function can be used to remove features
template<BNB_Features... donts>
constexpr BNB_SETTINGS remove_settings(BNB_SETTINGS settings = BNB_SETTINGS("111111"))
{
	(settings.set(idx(donts), false),...);
	return settings;
}
/**
 * solves the R||Cmax problem using branch-and-bound and the Lagrangian relaxation of machine-task-assignment, see paper for details and main.cpp for an example call
 *
 * @param p Matrix of processing times, p[i][j] is the processing time of task j by machine i.
 * @param settings bitset of settings to use, use e.g., remove_settings<BNB_Features::LARGE_LOCAL_SEARCH>() to use all but the large local search feature
 * @param time_limit is the maximal computation (wall)time in seconds, if the time limit is reached optimality cannot be guaranteed
 * @param x optional solution ouput, use nullptr to discard the solution 
 * @return (n_nodes, n_iter, z_ub, z_root_gap, z_root_lb_gap):
 * 			n_nodes - the number of branch-and-bound nodes explored
 *          n_iter - the total number of subgradient iterations
 *          z_ub - the best obtained solution value
 *          z_root_gap = (ub_root - lb_root) / ub_root
 *          z_root_lb_gap = (ub - lb_root) / ub  
 */
std::tuple<int, int, int, double, double> solve_bnb(const std::vector<std::vector<int>>& p, const BNB_SETTINGS settings = BNB_SETTINGS("111111"), const double time_limit = 120.0, std::vector<int>* x = nullptr);

}; // namespace  rcmax