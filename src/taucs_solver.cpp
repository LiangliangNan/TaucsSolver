#include "taucs_solver.h"
#include "taucs_matrix.h"
#include "taucs_util.h"
#include <iostream>


#define  TAUCS_CORE_DOUBLE
extern "C" {
#include <taucs.h>
}


bool TaucsSolver::solve_symmetry(const TaucsMatrix& matrix, 
								 const std::vector<double>& rhs, 
								 std::vector<double>& result)
{
	int num_row = matrix.row_dimension();
	int num_col = matrix.column_dimension();

	if (num_row != num_col) {
		std::cout << title() << "num_row != num_col" << std::endl;
		return false;
	}
	
	if (num_row != rhs.size()) {
		std::cout << title() << "num_row != rhs.size()" << std::endl;
		return false;
	}

	//////////////////////////////////////////////////////////////////////////

	// A
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)matrix.get_taucs_matrix();

	// X
	result.resize(num_col);

	//////////////////////////////////////////////////////////////////////////

	// solve
	void* F = NULL;
	char* factor[] = {"taucs.factor.LLT=true", NULL};
	char* solve[]  = {"taucs.factor=false", NULL};
	void* opt_arg[] = { NULL };

	/* this should work, factor, solve, free */
	int slove_rc = taucs_linsolve(A, &F, 1, &(result[0]), (void*)&(rhs[0]), factor, opt_arg);
	if (slove_rc != TAUCS_SUCCESS)
		std::cout << title() << "solve failed!" << std::endl;

	int free_rc = taucs_linsolve(NULL, &F, 0, NULL, NULL, factor, opt_arg);
	if (free_rc != TAUCS_SUCCESS)
		std::cout << title() << "free failed!" << std::endl;

	// clean
	//taucs_ccs_free(A); // A will be free by TaucsMatrix

	return (slove_rc == TAUCS_SUCCESS) && (free_rc == TAUCS_SUCCESS);
}


bool TaucsSolver::solve_non_symmetry(const TaucsMatrix& matrix, 
								   const std::vector<double>& rhs, 
								   std::vector<double>& result)
{
	int num_row = matrix.row_dimension();
	int num_col = matrix.column_dimension();

	if (num_row != num_col) {
		std::cout << title() << "num_row != num_col" << std::endl;
		return false;
	}
	
	if (num_row != rhs.size()) {
		std::cout << title() << "num_row != rhs.size()" << std::endl;
		return false;
	}

	//////////////////////////////////////////////////////////////////////////

	// A
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)matrix.get_taucs_matrix();

	// X
	result.resize(num_col);

	//////////////////////////////////////////////////////////////////////////

	int*    perm = NULL;
	int*    invperm = NULL;

	// ordering
	taucs_ccs_order(A, &perm, &invperm,	"colamd");
	if ( perm == NULL || invperm == NULL) {
		std::cout << title() << "ordering failed" << std::endl;
		return false;
	}

	taucs_io_handle* LU = taucs_io_create_multifile("taucs.L");
	if (LU == NULL) {
		std::cout << title() << "can not create multifile" << std::endl;
		return false;
	}

	// factorization
	int memory_mb = int(taucs_available_memory_size() / 1048576.0);
	int rc = taucs_ooc_factor_lu(A, perm, LU, memory_mb * 1048576.0);
	if (rc != TAUCS_SUCCESS) {
		std::cout << title() << "factorization failed" << std::endl;
		return false;
	}

	// solve
	rc = taucs_ooc_solve_lu(LU, &(result[0]), (void*)&(rhs[0]));
	if (rc != TAUCS_SUCCESS) {
		std::cout << title() << "solving failed" << std::endl;
		return false;
	}

	// delete the temporal multifile
	rc = taucs_io_delete(LU);
	if (rc != TAUCS_SUCCESS) {
		std::cout << title() << "delete multifile file failed" << std::endl;
		return false;
	}

	// clean
	//taucs_ccs_free(A); // A will be free by TaucsMatrix

	return (rc == TAUCS_SUCCESS);
}


bool TaucsSolver::solve_linear_least_square(const TaucsMatrix& matrix, 
											const std::vector<double>& rhs, 
											std::vector<double>& result)
{
	int num_row = matrix.row_dimension();
	int num_col = matrix.column_dimension();
	
	if (num_row < num_col) {
		std::cout << title() << "num_row < num_col" << std::endl;
		return false;
	}

	if (num_row != rhs.size()) {
		std::cout << title() << "num_row != rhs.size()" << std::endl;
		return false;
	}

	//////////////////////////////////////////////////////////////////////////

	// A
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)matrix.get_taucs_matrix();
	taucs_ccs_matrix* At = TaucsUtil::MatrixTranspose(A);
	taucs_ccs_matrix* AtA = TaucsUtil::Mul2NonSymmMatSymmResult(At, A);

	std::vector<double> AtB(num_col);
	TaucsUtil::MulNonSymmMatrixVector(At, &(rhs[0]), &(AtB[0]));

	// X
	result.resize(num_col);

	//////////////////////////////////////////////////////////////////////////

	// solve
	void* F = NULL;
	char* factor[] = {"taucs.factor.LLT=true", NULL};
	char* solve[]  = {"taucs.factor=false", NULL};
	void* opt_arg[] = { NULL };

	int slove_rc = taucs_linsolve(AtA, &F, 1, &(result[0]), &(AtB[0]), factor, opt_arg);
	if (slove_rc != TAUCS_SUCCESS)
		std::cout << title() << "solve failed" << std::endl;

	int free_rc = taucs_linsolve(NULL, &F, 0, NULL, NULL, factor, opt_arg);
	if (free_rc != TAUCS_SUCCESS) 
		std::cout << title() << "free failed" << std::endl;

	// clean
	//taucs_ccs_free(A); // A will be free by TaucsMatrix
	taucs_ccs_free(At);
	taucs_ccs_free(AtA);

	return (slove_rc == TAUCS_SUCCESS) && (free_rc == TAUCS_SUCCESS);
}


//////////////////////////////////////////////////////////////////////////
// api for array rhs

bool TaucsSolver::solve_symmetry(const TaucsMatrix& matrix, 
								 const std::vector<std::vector<double>>& rhs, 
								 std::vector<std::vector<double>>& result)
{
	int num_row = matrix.row_dimension();
	int num_col = matrix.column_dimension();

	if (num_row != num_col) {
		std::cout << title() << "num_row != num_col" << std::endl;
		return false;
	}

	for (unsigned int i=0; i<rhs.size(); ++i) {
		if (num_row != rhs[i].size()) {
			std::cout << title() << "num_row != rhs.size()" << std::endl;
			return false;
		}	
	}

	//////////////////////////////////////////////////////////////////////////

	// A
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)matrix.get_taucs_matrix();

	// X
	result.resize(rhs.size());
	for (unsigned int i=0; i<result.size(); ++i)
		result[i].resize(num_col);

	//////////////////////////////////////////////////////////////////////////
	// first factor, then solve, then free 

	void* F = NULL;
	char* factor[] = {"taucs.factor.LLT=true", NULL};
	char* solve[]  = {"taucs.factor=false", NULL};
	void* opt_arg[] = { NULL };

	// factor
	int factor_rc = taucs_linsolve(A, &F, 0, NULL, NULL, factor, opt_arg);
	if (factor_rc != TAUCS_SUCCESS) 
		std::cout << title() << "factorization failed" << std::endl;
	
	// solve
	int solve_rc = TAUCS_ERROR;
	for (unsigned int i=0; i<rhs.size(); ++i) {
		solve_rc = taucs_linsolve(A, &F, 1, &(result[i][0]), (void*)&(rhs[i][0]), solve, opt_arg);
		if (solve_rc != TAUCS_SUCCESS) {
			std::cout << title() << "solve for the " << i << "th vector failed" << std::endl;
			break;
		}
	}

	int free_rc = taucs_linsolve(NULL, &F, 0, NULL, NULL, factor, opt_arg);
	if (free_rc != TAUCS_SUCCESS)
		std::cout << title() << "free failed!" << std::endl;

	// clean
	//taucs_ccs_free(A); // A will be free by TaucsMatrix

	return (solve_rc == TAUCS_SUCCESS) && (free_rc == TAUCS_SUCCESS);
}


bool TaucsSolver::solve_non_symmetry(const TaucsMatrix& matrix, 
								   const std::vector<std::vector<double>>& rhs, 
								   std::vector<std::vector<double>>& result)
{
	int num_row = matrix.row_dimension();
	int num_col = matrix.column_dimension();

	if (num_row != num_col) {
		std::cout << title() << "num_row != num_col" << std::endl;
		return false;
	}

	for (unsigned int i=0; i<rhs.size(); ++i) {
		if (num_row != rhs[i].size()) {
			std::cout << title() << "num_row != rhs.size()" << std::endl;
			return false;
		}	
	}

	//////////////////////////////////////////////////////////////////////////

	// A
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)matrix.get_taucs_matrix();

	// X
	result.resize(rhs.size());
	for (unsigned int i=0; i<result.size(); ++i)
		result[i].resize(num_col);

	//////////////////////////////////////////////////////////////////////////

	int*    perm = NULL;
	int*    invperm = NULL;

	// ordering
	bool order_rc = true;
	taucs_ccs_order(A, &perm, &invperm,	"colamd");
	if ( perm == NULL || invperm == NULL) {
		std::cout << title() << "ordering failed" << std::endl;
		order_rc = false;
	}

	bool multfile_rc = true;
	taucs_io_handle* LU = taucs_io_create_multifile("taucs.L");
	if (LU == NULL) {
		std::cout << title() << "can not create multifile" << std::endl;
		multfile_rc = false;
	}

	// factorization
	int memory_mb = int(taucs_available_memory_size() / 1048576.0);
	int factor_rc = taucs_ooc_factor_lu(A, perm, LU, memory_mb * 1048576.0);
	if (factor_rc != TAUCS_SUCCESS)
		std::cout << title() << "factorization failed" << std::endl;

	// solve
	int solve_rc = TAUCS_ERROR;
	for (unsigned int i=0; i<rhs.size(); ++i) {
		solve_rc = taucs_ooc_solve_lu(LU, &(result[i][0]), (void*)&(rhs[i][0]));
		if (solve_rc != TAUCS_SUCCESS) {
			std::cout << title() << "solve for the " << i << "th vector failed" << std::endl;
			break;
		}
	}

	// delete the temporal multifile
	int delete_rc = taucs_io_delete(LU);
	if (delete_rc != TAUCS_SUCCESS) 
		std::cout << title() << "delete multifile file failed" << std::endl;

	// clean
	//taucs_ccs_free(A); // A will be free by TaucsMatrix

	return order_rc && multfile_rc && (factor_rc == TAUCS_SUCCESS) && (solve_rc == TAUCS_SUCCESS) && (delete_rc == TAUCS_SUCCESS);
}


bool TaucsSolver::solve_linear_least_square(const TaucsMatrix& matrix, 
											const std::vector<std::vector<double>>& rhs, 
											std::vector<std::vector<double>>& result)
{
	int num_row = matrix.row_dimension();
	int num_col = matrix.column_dimension();

	if (num_row < num_col) {
		std::cout << title() << "num_row < num_col" << std::endl;
		return false;
	}

	for (unsigned int i=0; i<rhs.size(); ++i) {
		if (num_row != rhs[i].size()) {
			std::cout << title() << "num_row != rhs.size()" << std::endl;
			return false;
		}	
	}

	//////////////////////////////////////////////////////////////////////////

	// AtA
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)matrix.get_taucs_matrix();
	taucs_ccs_matrix* At = TaucsUtil::MatrixTranspose(A);
	taucs_ccs_matrix* AtA = TaucsUtil::Mul2NonSymmMatSymmResult(At, A);

	// AtB
	std::vector<std::vector<double>> AtB(rhs.size());
	for (unsigned int i=0; i<rhs.size(); ++i) {
		AtB[i].resize(num_col);
		TaucsUtil::MulNonSymmMatrixVector(At, &(rhs[i][0]), &(AtB[i][0]));
	}

	// X
	result.resize(rhs.size());
	for (unsigned int i=0; i<rhs.size(); ++i) {
		result[i].resize(num_col);
	}

	//////////////////////////////////////////////////////////////////////////
	// first factor, then solve, then free 

	void* F = NULL;
	char* factor[] = {"taucs.factor.LLT=true", NULL};
	char* solve[]  = {"taucs.factor=false", NULL};
	void* opt_arg[] = { NULL };
	
	// factor
	int factor_rc = taucs_linsolve(AtA, &F, 0, NULL, NULL, factor, opt_arg);
	if (factor_rc != TAUCS_SUCCESS) 
		std::cout << title() << "factorization failed" << std::endl;

	// solve
	int solve_rc = TAUCS_ERROR;
	for (unsigned int i=0; i<AtB.size(); ++i) {
		solve_rc = taucs_linsolve(AtA, &F, 1, &(result[i][0]), &(AtB[i][0]), solve, opt_arg);
		if (solve_rc != TAUCS_SUCCESS) {
			std::cout << title() << "solve for the " << i << "th vector failed" << std::endl;
			break;
		}
	}

	// free
	int free_rc = taucs_linsolve(NULL, &F, 0, NULL, NULL, factor, opt_arg);
	if (free_rc != TAUCS_SUCCESS)
		std::cout << title() << "free failed" << std::endl;

	// clean
	//taucs_ccs_free(A); // A will be free by TaucsMatrix
	taucs_ccs_free(At);
	taucs_ccs_free(AtA);

	return (factor_rc == TAUCS_SUCCESS) && (solve_rc == TAUCS_SUCCESS) && (free_rc == TAUCS_SUCCESS);
}