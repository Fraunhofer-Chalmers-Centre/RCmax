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
#include <vector>

namespace rcmax {
	/**
	* Greedy heuristic that assigns one task at the time to the first free machine
	* @param p the processing times p_ij
	* @return the solution vector x[j] = i
	*/
	std::vector<int> rcmax_heuristic(const std::vector<std::vector<int>>& p);

	/**
	* Greedy heuristic that chooses the minimal reduced cost solution for each task
	* @param cred the reduced costs
	* @return the solution vector x[j] = i
	*/ 
	std::vector<int> rcmax_heuristic(const std::vector<std::vector<double>>& cred);

	/**
	* Feasibility heuristic choosing the most common assignment for each task (highest frequency)
	* @param x_seq the frequency vector (fractional solution)
	* @return the solution vector x[j] = i
	*/
	std::vector<int> rcmax_heuristic_ergo(const std::vector<std::vector<double>>& x_seq);

	/**
	* Feasibility heuristic that makes an assignment based on the subproblem solution x_sub greedily by choosing the best of the available machines or the best of all if none is available
	* @param x the solution vector x[j] = i
	* @param x_sub the subproblem solution
	* @param p the processing times p_ij
	*/
	void rcmax_heuristic(std::vector<int>& x, const std::vector<std::vector<bool>>& x_sub, const std::vector<std::vector<int>>& p);

	/**
	* Local search heuristic, the neighbourhood is to reassign a task
	* @param x  the solution vector x[j] = i
	* @param p the processing times p_ij
	* @return the makespan
	*/
	template<class T>
	T local_search_heuristic(std::vector<int>& x, const std::vector<std::vector<T>>& p);
	
	/**
	* Local search heuristic, the neighbourhood is a 2-cycle, i.e.,
	*	For any j1 : i_max is machine of io_sol[j1]
	*	For any j2 : i_2 (machine of k12) of io_sol[j2]
	* @param x  the solution vector x[j] = i
	* @param p the processing times p_ij
	* @return the makespan
	*/
	template<class T>
	T local_search_heuristic_swap(std::vector<int>& x, const std::vector<std::vector<T>>& p);

	/**
	* Local search heuristic, the neighbourhood is to send two from imax and recive one in return, i.e.,
	*	For any j1 : i_max is machine of x[j1]
	*	For any j2 : i_max is machine of x[j2], j2!=j1
	*	For any j3 : i_2 of x[j2]
	* @param x  the solution vector x[j] = i
	* @param p the processing times p_ij
	* @return the makespan
	*/
	template<class T>
	T local_search_heuristic_swap21(std::vector<int>& x, const std::vector<std::vector<T>>& p);

	/**
	* Local search heuristic, the neighbourhood is defined by:
	*	For any j1 : i_max is machine of x[j1]
	*	For any k12 in A(j1,!i_max)
	*	For any j2 : i_2 (machine of k12) of io_sol[j2]
	*	For any k23 : A(j2,!i_2 !i_max)
	*	For any j3 : i_3 (machine of k23) of io_sol[j3]
	*	For any k31 : A(j3,i_max)
	* @param DO_4 extends to do a 4 cycle neighbourhood
	* @param x  the solution vector x[j] = i
	* @param p the processing times p_ij
	* @param do_large if false, swap21 is used instead
	* @param pert_off an offset that is used to settle ties
	* @return the makespan
	*/
	template<bool DO_4>
	int local_search_heuristic_cycle3(std::vector<int>& x, const std::vector<std::vector<int>>& p, const bool do_large, double& pert_off);
}; // namespace  rcmax