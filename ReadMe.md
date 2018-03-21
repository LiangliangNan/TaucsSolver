**TaucsSolver** is a wrapper class for [TAUCS](http://www.tau.ac.il/~stoledo/taucs/), a powerfull linear system solver. 
See [here](http://www.tau.ac.il/~stoledo/taucs/) for the license and the availability note.

The taucs_addon.h file is modified from SOMECODE (I forget where it comes from. Please remind me if you know).

---

### Available functions
 * bool solve_symmetry();
 * bool solve_non_symmetry();
 * bool solve_linear_least_square();
 
Note: Corresponding APIs are also included for solving a bunch of rhs for the same coefficient matrix.

---

### How to use ? 
Quite easy! See the examples in "example/test.cpp" :-)

---

Please feel free to contact me if any question or comment.

Liangliang Nan

liangliang.nan@gmail.com

Sept. 2009

