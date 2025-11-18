# FMSQP.m — A MATLAB Feasible-Mode SQP Solver

This repository provides **FMSQP.m**, a standalone MATLAB implementation of a **F**ilter-based **M**ulti-start **S**equential **Q**uadratic **P**rogramming (SQP) method for solving constrained nonlinear optimization problems.  
The solver is lightweight, readable, and designed for research prototyping.

---

## 1. Overview

`FMSQP.m` solves constrained optimization problems of the form:

$$
\begin{aligned}
\min_x\quad & f(x) \\
\text{s.t.}\quad & c_i(x)=0,\quad i=1,\dots,m_e, \\
& h_j(x)\le 0,\quad j=1,\dots,m_i.
\end{aligned}
$$

**Important:** the iterates produced by the algorithm are **not required** to be feasible at every step. When the original problem is infeasible, the solver returns a point that minimizes a suitable measure of constraint violation while keeping the objective value small (i.e., it attempts to find a point that best satisfies the constraints and has a low objective).

The implementation is fully written in MATLAB without external dependencies, making it easy to read and extend.

---

## 2. Key Features

- Pure MATLAB implementation in a **single file**: `FMSQP.m`
- Handles:
  - Nonlinear objectives
  - Equality constraints
  - Inequality constraints
- Uses gradient information
- Includes example problem files (`Ex3.*`, `Ex5.*`) and a main driver (`MainG.m`) for running a benchmark test problem library.

---

## 3. Function Interface

Call the solver using:

### Updated Function Interface (correct call format)

The actual solver signature is:

```matlab
function [x, fx, output] = FMSQP(X0, opts)
````
### Required User-Supplied Function Files

The solver relies on **four user-defined MATLAB function files** that define the objective and constraint system:

| File       | Description |
|------------|-------------|
| `funf.m`   | Returns the **objective function value** \( f(x) \). |
| `gradf.m`  | Returns the **gradient of the objective** \( \nabla f(x) \). |
| `func.m`   | Returns the **constraint function values** \( c(x) \le 0 \). |
| `gradc.m`  | Returns the **Jacobian (gradient matrix)** of the constraints \( \nabla c(x) \). |

**Inputs**

* `X0` — initial point (column vector).
* `opts` — options structure (optional). Recognized fields:

  * `opts.varbose` — display verbosity (default `0`).
  * `opts.epsilon` — stopping tolerance (default `1e-5`).
  * `opts.nmax` — maximum iterations (default `500`).

**Outputs**

* `x` — final iterate returned by the solver.
* `fx` — objective value at `x`.
* `output` — struct with diagnostic information.

### Output Structure (`output`)

The `output` struct returned by `FMSQP` contains the following fields:

| Field     | Description |
|-----------|-------------|
| `m`       | Total number of constraints (size of `c(x)`). |
| `vx`      | Final constraint violation at `x` (vector of `c(x)` values). |
| `nit`     | Number of iterations performed. |
| `nf`      | Number of objective function evaluations. |
| `ng`      | Number of gradient evaluations (objective and/or constraints). |
| `time`    | Total elapsed computation time in seconds (`toc`). |

**Example: inspecting output**

```matlab
fprintf('Number of iterations: %d\n', output.nit);
fprintf('Number of function evaluations: %d\n', output.nf);
fprintf('Number of gradient evaluations: %d\n', output.ng);
fprintf('Final constraint violation (max): %g\n', max(output.vx));
fprintf('Elapsed time: %g seconds\n', output.time);
````

---

## 4. Minimal Example

```matlab
% Objective function
fun = @(x) deal( x(1)^2 + x(2)^2, [2*x(1); 2*x(2)] );

% Initial points (each column is one starting point)
x0 = [0, 0; 
      1, 0; 
      0, 1; 
      1, 1]';

% Run solver
[x, lambda, info] = FMSQP(X0);
```

---

## 5. Algorithm Summary

`FMSQP.m` implements a **multi-start SQP framework combined with a filter technique**:

1. Generate multiple initial points and run SQP iterations from each.  
2. Evaluate objective, gradient, and constraint functions at the current iterate.  
3. Linearize inequality constraints and build a QP subproblem.  
4. Solve the QP to obtain a search direction.  
5. Perform step-size selection (line search).  
6. Update primal variables and multipliers.  
7. Use a **filter** to discard iterates that are dominated in terms of objective and constraint violation, helping the algorithm:  
   - Avoid infeasible stationary points  
   - Escape local minima that are not globally optimal  
8. Continue until tolerances for KKT residuals or constraint violation are satisfied or maximum iterations are reached.  

This combination of **multi-start strategy + filter method** allows the solver to explore multiple basins of attraction and improve the chance of finding a high-quality feasible solution or a best-effort point when the problem is infeasible.


---


## 6. Repository Structure

```
FMSQP/
├── FMSQP.m        % Main solver
├── MainG.m        % Benchmark example
├── Ex*          % Example problems (Ex3.10, Ex5.1, etc.)
└── README.md
```

---

## 7. Notes

* Intended for **small or medium-scale** optimization problems
* For large-scale applications, consider:

  * sparse matrices
  * structure-exploiting BFGS updates
  * iterative solvers

---

## 8. License

This project is licensed under the MIT License. You are free to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of this software, provided that the following conditions are met:

1. The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
2. The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement.

For the full license text, see the [LICENSE](LICENSE) file.

---

## 10. Contact

Author: **Wenhao Fu**
E-mail: wenhfu@usts.edu.cn
For questions or suggestions, please open an issue on GitHub.
