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

#include "LR1solver.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>

namespace rcmax {
namespace sgupdates{	
	/** 
	* ADS_param compute the deflection parameter of the Average Direction Strategy (ADS) (Sherali and Ulular, 1989)
	* @param dir the previous subgradient search direction
	* @param sub_grad the current subgradient
	* @return the ADS paramter
	*/
	double ADS_param(std::vector<double>& dir, const std::vector<double>& sub_grad)
	{
		const auto dir_norm2 = std::transform_reduce(dir.begin(), dir.end(), 0.0, std::plus<>(), [](const double v) { return v * v; });
		if (dir_norm2 < 1e-6)
			return 0.0;
		const auto sg_norm2 = std::transform_reduce(sub_grad.begin(), sub_grad.end(), 0.0, std::plus<>(), [](const double v) { return v * v; });
		return std::sqrt(sg_norm2) / std::sqrt(dir_norm2);
	}
	/** 	
	* MGT_param compute the deflection parameter of the Modified Gradient Technique (MGT) (Camerini et al., 1975)
	* @param dir the previous subgradient search direction
	* @param sub_grad the current subgradient
	* @param eta scale parameter in (0,2], where 0 is to use the subgradient as is, 1 is to project down in to the halfspace of the last halfspace, and 2 is to reflect in the halfspace
	* @return the ADS paramter
	*/
	double MGT_param(std::vector<double>& dir, const std::vector<double>& sub_grad, const double eta = 1.5)
	{
		const auto dir_norm2 = std::transform_reduce(dir.begin(), dir.end(), 0.0, std::plus<>(), [](const double v) { return v * v; });
		const double sg_dir_dot = std::transform_reduce(dir.begin(), dir.end(), sub_grad.begin(), 0.0);
		return std::max(0.0, -eta * sg_dir_dot / dir_norm2);
	}
	/**
	* MDS_update deflects the seachdirection according to the Modified Deflected Subgradient method (MSD) (Rachid Belgacem and Abdessamad Amir, 2018)
	* @param dir the previous subgradient search direction to be updated
	* @param sub_grad the current subgradient
	*/
	void MDS_update(std::vector<double>& dir, const std::vector<double>& sub_grad)
	{
		const double sg_dir_dot = std::transform_reduce(dir.begin(), dir.end(), sub_grad.begin(), 0.0);
		const auto dir_norm2 = std::transform_reduce(dir.begin(), dir.end(), 0.0, std::plus<>(), [](const double v) { return v * v; });
		const auto sg_norm2 = std::transform_reduce(sub_grad.begin(), sub_grad.end(), 0.0, std::plus<>(), [](const double v) { return v * v; });
		const auto alpha = std::max(0.0, -sg_dir_dot / std::sqrt(dir_norm2 * sg_norm2));
		const auto eta = 1.0 / (2.0 - alpha) - 1e-6;
		const auto deflect_param = (1.0 - alpha) * MGT_param(dir, sub_grad, eta) + alpha * ADS_param(dir, sub_grad);
		std::transform(dir.begin(), dir.end(), sub_grad.begin(), dir.begin(), [deflect_param](const double dir, const double sg) {return sg + deflect_param * dir; });
	}
} // namespace sgupdates

namespace LR1 {
	/**
	* The Euclidian projection of the vector u onto the unit simplex Σ u = 1, u ≥ 0, specialized from Nonlinear Programming: Theory and Algorithms by (Mokhtar S. Bazaraa et al., 2006)
	* @param u the vector to be projected
	*/
	void project_onto_unit_simplex(std::vector<double>& u) {
		auto diff = (1.0 - std::accumulate(u.begin(), u.end(), 0.0)) / u.size();
		for (auto& v : u) {
			v += diff;
		}
		bool all_ok = true;
		size_t n_free = u.size();
		for (auto& v : u) {
			if (v <= 0.0) {
				v = 0.0;
				all_ok = false;
				--n_free;
			}
		}
		while (!all_ok)
		{
			all_ok = true;
			diff = (1.0 - std::accumulate(u.begin(), u.end(), 0.0)) / n_free;
			for (auto& v : u) {
				if (v != 0.0) {
					v += diff;
					if (v < 0.0) {
						v = 0.0;
						all_ok = false;
						--n_free;
					}
				}
			}
		}
	}
	/**
	* Find the l = max(1e-3, argmin max_i ((1-l) zc_i  + l zn_i)) and set zc := (1-l) zc  + l zn, this to heuristically improve an upper bound corresponding to the master problem of the Dantzig-Wolfe decomposition
	*  @param z_conv is a convex combination of previous z vectors
	*  @param z_new is the current z vector, i.e., z_i = Σ_j p_ij x_ij
	*/
	void opt_conv_comb(std::vector<double>& z_conv, const std::vector<double>& z_new) {
		auto m = z_conv.size();
		auto a = 0.0, b = -1.0; // current left line a + l b
		auto c = -1.0, d = 1.0; // current right line c + l d
		auto zmax = 0.0;
		auto l = 0.0;
		for (size_t i = 0, iter = 0; iter < m; ++iter, ++i) {
			if (i == m) {
				i = 0;
			}
			// Check if new line cuts of current sol, e + l * f
			auto e = z_conv[i];
			auto f = z_new[i] - z_conv[i];
			if (e + l * f < zmax + 1e-8) {
				continue;
			}
			// switch basis
			if (f >= 0.0) {
				c = e;
				d = f;
			}
			else {
				a = e;
				b = f;
			}
			// update solution // a + l b = c + l d
			l = (a - c) / (d - b);
			if (l < 0.0) {
				l = 0.0;
				a = c;
			}
			else if (l > 1.0) {
				l = 1.0;
				c = a - d + b;
			}
			zmax = a + b * l;
			iter = 0;
		}
		l = std::max(l, 2e-3); // ensure a small step to escape non-differentiable parts
		std::transform(z_conv.begin(), z_conv.end(), z_new.begin(), z_conv.begin(), [l](const double zc, const double zn) {
			return (1 - l) * zc + l * zn;
			});
	}

	/**
	* Constructor
	* @param im is the number of machines
	* @param in is the number of tasks
	*/
	LR1Solver::LR1Solver(size_t im, size_t in)
		: m(im)
		, n(in)
		, cT(in, std::vector<int>(im))
		, cred(im, std::vector<double>(in, 0.0))
		, u(im, 1.0 / im)
		, v(in, 0.0)
		, x_sol(in) {
	}

	/**
	* maximize the Lagrangian dual function using the subgradient method
	* @param ci is the processing times, somtimes denoted p_ij
	* @return the best lower bound, which coincides with the value of LP relaxation, if the method converges (which is true for all problems in our test-bed)
	*/
	double LR1Solver::solve(const std::vector<std::vector<int>>& ci) {
		// use transpose for better cache locality
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				cT[j][i] = ci[i][j];
			}
		}
		std::vector<double> z(m, 0.0), z_conv(m);
		auto iter_before_div = 4 * (m + 15);
		size_t n_iter = 100 * (m + 15);
		auto delta = 0.5;
		if (!root) {
			iter_before_div = 2 * (m + 15);
			n_iter = 25 * (m + 15);
			delta = 1e-2;
		}
		size_t iter_since_last_div = 0;
		auto lb = 0.0;
		auto ub = std::numeric_limits<double>::infinity();
		auto ub_frac = ub;
		auto x = x_sol;
		auto best_u = u;
		auto u_dir = u;
		for (size_t iter = 0; iter < n_iter; ++iter_since_last_div, ++iter) {
			project_onto_unit_simplex(u);
			std::fill(z.begin(), z.end(), 0.0);
			auto lbi = 0.0;
			for (size_t j = 0; j < n; ++j) {
				auto minv = std::numeric_limits<double>::infinity();
				for (size_t i = 0; i < m; ++i) {
					if (auto vj = u[i] * cT[j][i]; vj < minv) {
						x[j] = static_cast<int>(i);
						minv = vj;
					}
				}
				z[x[j]] += cT[j][x[j]];
				lbi += minv;
			}
			auto ubi = *std::max_element(z.begin(), z.end());
			if (ubi < ub) {
				ub = ubi;
				x_sol = x;
				if (ub <= lb + 0.99) {
					break;
				}
			}
			if (lbi > lb + 1e-8 * lb) {
				lb = lbi;
				best_u = u;
				iter_since_last_div = 0;
				if (ub <= lb + 0.99 || ub_frac <= lb + lb * 1e-3) {
					break;
				}
			}
			if (iter_since_last_div > iter_before_div) {
				iter_since_last_div = 0;
				delta *= .5;
			}
			if (iter == 0) {
				z_conv = z;
			}
			else {
				opt_conv_comb(z_conv, z);
			}
			ub_frac = std::min(ub_frac, *std::max_element(z_conv.begin(), z_conv.end()));
			if (ub_frac <= lb + lb * 1e-3) {
				break;
			}
			auto mean_z = std::accumulate(z.begin(), z.end(), 0.0) / m;
			auto& sg = z;
			auto sg_norm2 = 0.0;
			for (auto& e : sg) {
				e -= mean_z;
				sg_norm2 += e * e;
			}
			if (sg_norm2 < 1e-5) {
				break;
			}
			if (iter == 0) {
				u_dir = sg;
			}
			else {
				sgupdates::MDS_update(u_dir, sg);
			}
			const auto tk = delta * (ub_frac - lbi) / sg_norm2;
			std::transform(u.begin(), u.end(), sg.begin(), u.begin(), [tk](const double u, const double d) {
				return u + tk * d; });
		}
		// compute reduced costs
		u = best_u;
		for (size_t j = 0; j < n; ++j) {
			auto minv = std::transform_reduce(cT[j].begin(), cT[j].end(), u.begin(), std::numeric_limits<double>::infinity(),
				[](double v1, double v2) {return std::min(v1, v2); },
				std::multiplies<>());
			v[j] = minv;
			for (size_t i = 0; i < m; ++i) {
				cred[i][j] = cT[j][i] * u[i] - minv;
			}
		}
		return lb;
	}
}//namespace LR1
}//namespace rcmax