#ifndef _TAUCS_SOLVER_H_
#define _TAUCS_SOLVER_H_

/********************************************************************
	Version:	 1.5

	created:	 9:11:2009   9:04
	filename: 	 taucs_solver.h
	author:		 Liangliang Nan
	contact:     liangliang.nan@gmail.com
	purpose:	 Wrapper of TAUCS for easy use 

*********************************************************************

	How to use: see the examples

*********************************************************************

	Change log:
	------------------------------------------------
	Jul  9, 2014 - better encapsulation (all number types are double)

	Dec 29, 2012 - add API for an array of rhs using the same matrix
			     - delete the multifile for non-symmetry solver

*********************************************************************/

#include <vector>
#include <string>



class TaucsMatrix;

class TaucsSolver
{
public:
	static std::string title() { return "[TaucsSolver]: "; }

	// solve for "A*x=b"
	// A: the symmetry coefficient matrix, 
	// b: the right side column vector
	// x: the result
	static bool solve_symmetry(
		const TaucsMatrix& A, 
		const std::vector<double>& b, 
		std::vector<double>& x
		);
	
	
	// solve for "A*x=b"
	// A: the non-symmetry coefficient matrix, 
	// b: the right side column vector
	// x: the result
	static bool solve_non_symmetry(
		const TaucsMatrix& A, 
		const std::vector<double>& b, 
		std::vector<double>& x
		);


	// solve for "A*x=b" in least square sence
	// A: the coefficient m * n matrix (m >= n)
	// b: the right side column vector
	// x: the result
	static bool solve_linear_least_square(
		const TaucsMatrix& A, 
		const std::vector<double>& b, 
		std::vector<double>& x
		);

	//////////////////////////////////////////////////////////////////////////

	// API for an array of rhs using the same matrix

	// solve for "A*x=b"
	// A: the symmetry coefficient matrix, 
	// B: the array of right side column vector
	// X: the array of result vectors
	static bool solve_symmetry(
		const TaucsMatrix& A, 
		const std::vector<std::vector<double>>& B, 
		std::vector<std::vector<double>>& X
		);

	// solve for "A*x=b"
	// A: the non-symmetry coefficient matrix, 
	// B: the array of right side column vector
	// X: the array of result vectors
	static bool solve_non_symmetry(
		const TaucsMatrix& A, 
		const std::vector<std::vector<double>>& B, 
		std::vector<std::vector<double>>& X
		);

	// solve for "A*x=b" in least square sence
	// A: the coefficient m * n matrix (m >= n)
	// B: the array of right side column vector
	// X: the array of result vectors
	static bool solve_linear_least_square(
		const TaucsMatrix& A, 
		const std::vector<std::vector<double>>& B, 
		std::vector<std::vector<double>>& X
		);
};


#endif