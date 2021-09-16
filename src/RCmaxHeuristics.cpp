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

#include "RCmaxHeuristics.h"

#include <algorithm>
#include <cmath>
#include <numeric>
namespace rcmax {

std::vector<int> rcmax_heuristic(const std::vector<std::vector<int>>& c)
{
	const auto m = c.size();
	const auto n = c.front().size();
	std::vector<int> z(m, 0);
	std::vector<int> x(n, 0);
	for (size_t j = 0; j < n; ++j) {
		size_t i_min = 0;
		for (size_t i = 1; i < m; ++i) {
			if (z[i] + c[i][j] < z[i_min] + c[i_min][j]) {
				i_min = i;
			}
		}
		z[i_min] += c[i_min][j];
		x[j] = static_cast<int>(i_min);
	}
	return x;
}
std::vector<int> rcmax_heuristic(const std::vector<std::vector<double>>& c_red)
{
	const auto n = c_red.front().size();
	std::vector<int> x(n, 0);
	for (size_t j = 0; j < n; ++j) {
		x[j] = static_cast<int>(std::distance(c_red.begin(),
											  std::min_element(c_red.begin(), c_red.end(),
															   [j] (const auto& l, const auto& r) {return l[j] < r[j]; })));
	}
	return x;
}

std::vector<int> rcmax_heuristic_ergo(const std::vector<std::vector<double>>& x_seq)
{
	const auto n = x_seq.front().size();
	std::vector<int> x(n, 0);
	for (size_t j = 0; j < n; ++j) {
		x[j] = static_cast<int>(std::distance(x_seq.begin(),
											  std::max_element(x_seq.begin(), x_seq.end(),
															   [j] (const auto& l, const auto& r) {return l[j] < r[j]; })));
	}
	return x;
}

void rcmax_heuristic(std::vector<int>& x, const std::vector<std::vector<bool>>& x_sub, const std::vector<std::vector<int>>& c)
{
	const auto m = c.size();
	const auto n = c.front().size();
	std::vector<int> z(m, 0);
	std::fill(x.begin(), x.end(), static_cast<int>(m));
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (x_sub[i][j]) {
				++x[j];
			}
		}
	}
	for (size_t j = 0; j < n; ++j) {
		auto n_assigned = 0;
		for (size_t i = 0; i < m; ++i) {
			if (x_sub[i][j]) {
				x[j] = static_cast<int>(i);
				z[i] += c[i][j];
				++n_assigned;
			}
		}
		if (n_assigned > 1) {
			x[j] = static_cast<int>(m + 1);
		}
	}
	for (size_t j = 0; j < n; ++j) {
		if (x[j] < m)
			continue;
		size_t i_min = 0;
		if(x[j] > m){
			auto z_min = *std::max_element(z.begin(), z.end()) + 1;
			for (size_t i = 0; i < m; ++i) {
				if (x_sub[i][j] && z[i] < z_min) {
					i_min = i;
					z_min = z[i];
					z[i] -= c[i][j];
				}
			}
		}
		else {
			for (size_t i = 1; i < m; ++i) {
				if (z[i] + c[i][j] < z[i_min] + c[i_min][j]) {
					i_min = i;
				}
			}
		}
		z[i_min] += c[i_min][j];
		x[j] = static_cast<int>(i_min);
	}
	return;
}

template<class T>
T local_search_heuristic(std::vector<int>& x, const std::vector<std::vector<T>>& c)
{
	const auto m = c.size();
	const auto n = c.front().size();
	std::vector<T> z(m, 0);
	for (size_t j = 0; j < n; ++j) {
		z[x[j]] += c[x[j]][j];
	}
	auto zm = *std::max_element(z.begin(), z.end());
	auto zi = zm;

	//The neighbour hood is to reassign a task
	for (size_t j = 0, n_checked = 0; n_checked < n; ++j, ++n_checked) {
		if (j == n)
			j = 0;

		auto best_i = x[j];
		if (z[best_i] < zi)
			continue;
		for (size_t i = 0; i < m; ++i) {
			if (z[i] + c[i][j] < zm) {
				zm = z[i] + c[i][j];
				best_i = static_cast<int>(i);
			}
		}

		if (zm < zi) {
			z[x[j]] -= c[x[j]][j];
			x[j] = best_i;
			z[x[j]] += c[x[j]][j];
			zm = zi = *std::max_element(z.begin(), z.end());
			n_checked = 0;
		}
	}
	return zi;
}

/* Neighbourhood is a 2-cycle.
i.e.
	For any j1 : i_max is machine of io_sol[j1]
	For any j2 : i_2 (machine of k12) of io_sol[j2]
*/
template<class T>
T local_search_heuristic_swap(std::vector<int>& x, const std::vector<std::vector<T>>& c)
{	//The neighbour hood is to swap two task, sending one from the highest contributor and reciving one in return
	local_search_heuristic(x, c);
	const auto m = c.size();
	const auto n = c.front().size();
	std::vector<T> z(m, 0);
	for (size_t j = 0; j < n; ++j) {
		z[x[j]] += c[x[j]][j];
	}
	auto zi = *std::max_element(z.begin(), z.end());

	std::vector<std::vector<size_t>> tasks(m);
	for (size_t j = 0; j < n; ++j) {
		tasks[x[j]].push_back(j);
	}

	for (size_t j1 = 0, n_checked = 0; n_checked < n; ++j1, ++n_checked) {
		if (j1 == n)
			j1 = 0;
		const auto i1 = x[j1];
		if (z[i1] < zi)
			continue;
		const auto c_from_i1 = c[i1][j1];
		auto improvement = 0;
		auto best_j2 = n;
		auto best_i2 = m;
		for (size_t i2 = 0; i2 < m; ++i2) {
			if (i2 == i1)
				continue;
			const auto c_slack_i2 = zi - z[i2];
			const auto c_to_i2 = c[i2][j1];
			if (c_to_i2 < c_slack_i2)
				continue;
			for (auto j2 : tasks[i2]) {
				if (c_from_i1 > improvement + c[i1][j2] && c_slack_i2 + c[i2][j2] > improvement + c_to_i2) {
					improvement = std::min(c_from_i1 - c[i1][j2], c_slack_i2 + c[i2][j2] - c_to_i2);
					best_j2 = j2;
					best_i2 = i2;
				}
			}
		}
		if (improvement > 0) {
			n_checked = 0;
			x[j1] = static_cast<int>(best_i2);
			x[best_j2] = i1;
			local_search_heuristic(x, c);
			for (auto& v : z) {
				v = 0;
			}
			for (size_t j = 0; j < n; ++j) {
				z[x[j]] += c[x[j]][j];
			}
			zi = *std::max_element(z.begin(), z.end());

			for (auto& vec : tasks) {
				vec.clear();
			}
			for (size_t j = 0; j < n; ++j) {
				tasks[x[j]].push_back(j);
			}
		}
	}
	return zi;
}


/* Neighbourhood is a send two from imax and recive one in return
i.e.
For any j1 : i_max is machine of x[j1]
For any j2 : i_max is machine of x[j2], j2!=j1
For any j3 : i_2 of x[j2]
*/
template<class T>
T local_search_heuristic_swap21(std::vector<int>& x, const std::vector<std::vector<T>>& c)
{
	local_search_heuristic_swap(x, c);
	const auto m = c.size();
	const auto n = c.front().size();
	std::vector<T> z(m, 0);
	for (size_t j = 0; j < n; ++j) {
		z[x[j]] += c[x[j]][j];
	}
	auto zm = *std::max_element(z.begin(), z.end());

	std::vector<std::vector<size_t>> tasks(m);
	for (size_t j = 0; j < n; ++j) {
		tasks[x[j]].push_back(j);
	}

	//send j1 and j2 from max_i recieve j3 from i
	for (size_t j1 = 0, n_checked = 0; n_checked < n; ++j1, ++n_checked) {
		if (j1 == n)
			j1 = 0;
		const auto i1 = x[j1];
		if (z[i1] < zm)
			continue;
		auto improvement = 0;
		auto best_j2 = n, best_j3 = n;
		auto best_i2 = m;

		for (auto j2 : tasks[i1]) {
			if (j1 == j2)
				continue;
			const auto c_from_i1 = c[i1][j1] + c[i1][j2];
			for (size_t i2 = 0; i2 < m; ++i2) {
				if (i2 == i1)
					continue;
				const auto c_slack_i2 = zm - z[i2];
				const auto c_to_i2 = c[i2][j1] + c[i2][j2];
				if (c[i2][j1] < c_slack_i2 || c[i2][j2] < c_slack_i2)
					continue;
				for (auto j3 : tasks[i2]) {
					if (c_from_i1 > improvement + c[i1][j3] && c_slack_i2 + c[i2][j3] > improvement + c_to_i2) {
						improvement = std::min(c_from_i1 - c[i1][j3], c_slack_i2 + c[i2][j3] - c_to_i2);
						best_j2 = j2; best_j3 = j3;
						best_i2 = i2;
					}
				}
			}
		}
		if (improvement > 0) {
			n_checked = 0;
			x[j1] = static_cast<int>(best_i2);
			x[best_j2] = static_cast<int>(best_i2);
			x[best_j3] = i1;

			local_search_heuristic_swap(x, c);
			for (auto& v : z) {
				v = 0;
			}
			for (size_t j = 0; j < n; ++j) {
				z[x[j]] += c[x[j]][j];
			}
			zm = *std::max_element(z.begin(), z.end());
			for (auto& vec : tasks) {
				vec.clear();
			}
			for (size_t j = 0; j < n; ++j) {
				tasks[x[j]].push_back(j);
			}
		}
	}
	return zm;
}
/*Neighbourhood defined by:
	For any j1 : i_max is machine of x[j1]
	For any k12 in A(j1,!i_max)
	For any j2 : i_2 (machine of k12) of io_sol[j2]
	For any k23 : A(j2,!i_2 !i_max)
	For any j3 : i_3 (machine of k23) of io_sol[j3]
	For any k31 : A(j3,i_max)
*/
template<bool DO_4>
int local_search_heuristic_cycle3(std::vector<int>& x, const std::vector<std::vector<int>>& c_in, const bool do_large, double& pert_ofs)
{
	if (!do_large) {
		return local_search_heuristic_swap(x, c_in);
	}
	local_search_heuristic_swap(x, c_in);
	const auto m = c_in.size();
	const auto n = c_in.front().size();
	std::vector<std::vector<long long>> c(m, std::vector<long long>(n));
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			c[i][j] = c_in[i][j] * static_cast<long long>(n) * 100 + static_cast<long long>(5 * (1.0 + std::sin(static_cast<double>(i) * n + j + pert_ofs)));
		}
	}
	pert_ofs += 37.0;

	local_search_heuristic_swap21(x, c);
	std::vector<long long> z(m, 0);
	for (size_t j = 0; j < n; ++j) {
		z[x[j]] += c[x[j]][j];
	}
	auto zm = *std::max_element(z.begin(), z.end());

	std::vector<std::vector<size_t>> tasks(m);
	for (size_t j = 0; j < n; ++j) {
		tasks[x[j]].push_back(j);
	}

	for (size_t j1 = 0, n_checked = 0; n_checked < n; ++j1, ++n_checked) {
		if (j1 == n)
			j1 = 0;
		const auto i1 = static_cast<size_t>(x[j1]);
		if (z[i1] < zm)
			continue;
		const auto c_from_i1 = c[i1][j1];
		auto improvement = 0;
		auto best_i1 = m, best_i2 = m, best_i3 = m, best_i4 = m;
		auto best_j2 = n, best_j3 = n, best_j4 = n;

		for (size_t i2 = 0; i2 < m; ++i2) {
			if (i2 == i1)
				continue;
			const auto c_slack_i2 = zm - z[i2];
			const auto c_to_i2 = c[i2][j1];
			if (c_to_i2 < c_slack_i2)
				continue;

			for (auto j2 : tasks[i2]) {
				const auto c_from_i2 = c[i2][j2];
				if (c_to_i2 + improvement >= c_from_i2 + c_slack_i2)
					continue;

				for (size_t i3 = 0; i3 < m; ++i3) {
					if (i3 == i2 || i3 == i1)
						continue;
					const auto c_slack_i3 = zm - z[i3];
					const auto c_to_i3 = c[i3][j2];
					//Check if nothing needs to be sent back to i1 (not complete the cycle) //not always reducable to send operations
					if (c_to_i3 + improvement < c_slack_i3) {
						best_j2 = j2;	best_j3 = n; best_j4 = n;
						best_i1 = i2, best_i2 = i3; best_i3 = m; best_i4 = m;
						improvement = std::min({ c_from_i1, c_from_i2 - c_to_i2 + c_slack_i2, c_slack_i3 - c_to_i3 });
					}
					for (auto j3 : tasks[i3]) {
						const auto c_from_i3 = c[i3][j3];
						if (c_to_i3 + improvement >= c_from_i3 + c_slack_i3)
							continue;
						{
							// send back to i1
							const auto c_to_i1 = c[i1][j3];
							if (c_from_i1 > improvement + c_to_i1) {
								best_j2 = j2;	best_j3 = j3; best_j4 = n;
								best_i1 = i2, best_i2 = i3;	best_i3 = i1;	best_i4 = m;
								improvement = std::min({ c_from_i1 - c_to_i1, c_from_i2 - c_to_i2 + c_slack_i2, c_slack_i3 + c_from_i3 - c_to_i3 });
							}
							// send back to i2
							if (c_from_i2 + c_slack_i2 > improvement + c_to_i2 + c[i2][j3]) {
								best_j2 = j2;	best_j3 = j3; best_j4 = n;
								best_i1 = i2, best_i2 = i3;	best_i3 = i2;	best_i4 = m;
								improvement = std::min({ c_from_i1, c_from_i2 - c_to_i2 - c[i2][j3] + c_slack_i2, c_slack_i3 + c_from_i3 - c_to_i3 });
							}
						}
						if constexpr (DO_4) {
							for (size_t i4 = 0; i4 < m; ++i4) {
								if (i4 == i3 || i4 == i2 || i4 == i1)
									continue;
								const auto c_slack_i4 = zm - z[i4];
								const auto c_to_i4 = c[i4][j3];
								//Check if nothing needs to be sent back to i1 (not complete the cycle) //not always reducable to send operations
								if (c_to_i4 + improvement < c_slack_i4) {
									best_j2 = j2;	best_j3 = j3; best_j4 = n;
									best_i1 = i2, best_i2 = i3;	best_i3 = i4;	best_i4 = m;
									improvement = std::min({ c_from_i1, c_from_i2 - c_to_i2 + c_slack_i2, c_from_i3 + c_slack_i3 - c_to_i3, c_slack_i4 - c_to_i4 });
								}
								for (auto j4 : tasks[i4]) {
									const auto c_from_i4 = c[i4][j4];
									if (c_to_i4 + improvement >= c_from_i4 + c_slack_i4)
										continue;
									// send back to i1
									const auto c_to_i1 = c[i1][j4];
									if (c_from_i1 > improvement + c_to_i1) {
										best_j2 = j2;	best_j3 = j3; best_j4 = j4;
										best_i1 = i2, best_i2 = i3;	best_i3 = i4;	best_i4 = i1;
										improvement = std::min({ c_from_i1 - c_to_i1, c_from_i2 - c_to_i2 + c_slack_i2, c_slack_i3 + c_from_i3 - c_to_i3, c_slack_i4 + c_from_i4 - c_to_i4 });
									}
									// send back to i2
									if (c_from_i2 + c_slack_i2 > improvement + c_to_i2 + c[i2][j4]) {
										best_j2 = j2;	best_j3 = j3; best_j4 = j4;
										best_i1 = i2, best_i2 = i3;	best_i3 = i4;	best_i4 = i2;
										improvement = std::min({ c_from_i1, c_from_i2 - c_to_i2 - c[i2][j4] + c_slack_i2, c_slack_i3 + c_from_i3 - c_to_i3, c_slack_i4 + c_from_i4 - c_to_i4 });
									}
									// send back to i3
									if (c_from_i3 + c_slack_i3 > improvement + c_to_i3 + c[i3][j4]) {
										best_j2 = j2;	best_j3 = j3; best_j4 = j4;
										best_i1 = i2, best_i2 = i3;	best_i3 = i4;	best_i4 = i3;
										improvement = std::min({ c_from_i1, c_from_i2 - c_to_i2 + c_slack_i2, c_slack_i3 + c_from_i3 - c_to_i3 - c[i3][j4], c_slack_i4 + c_from_i4 - c_to_i4 });
									}
								}
							}
						}
					}
				}
			}
		}
		if (improvement > 0) {
			n_checked = 0;
			x[j1] = static_cast<int>(best_i1);
			x[best_j2] = static_cast<int>(best_i2);
			if (best_j3 < n)
				x[best_j3] = static_cast<int>(best_i3);
			if (best_j4 < n)
				x[best_j4] = static_cast<int>(best_i4);
			local_search_heuristic_swap21(x, c);
			for (auto& v : z) {
				v = 0;
			}
			for (size_t j = 0; j < n; ++j) {
				z[x[j]] += c[x[j]][j];
			}
			zm = *std::max_element(z.begin(), z.end());
			for (auto& vec : tasks) {
				vec.clear();
			}
			for (size_t j = 0; j < n; ++j) {
				tasks[x[j]].push_back(j);
			}
		}
	}

	for (auto& v : z) {
		v = 0;
	}
	for (size_t j = 0; j < n; ++j) {
		z[x[j]] += c_in[x[j]][j];
	}
	zm = *std::max_element(z.begin(), z.end());
	return static_cast<int>(zm);
}
template int local_search_heuristic_cycle3<false>(std::vector<int>& , const std::vector<std::vector<int>>& , const bool,double& );
template int local_search_heuristic_cycle3<true>(std::vector<int>&, const std::vector<std::vector<int>>&, const bool,double&);
} // namespace rcmax