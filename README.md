# FMSQP.m — A MATLAB Feasible-Mode SQP Solver

This repository provides **FMSQP.m**, a standalone MATLAB implementation of a **F**ilter-based **M**ulti-start **S**equential **Q**uadratic **P**rogramming (SQP) method for solving constrained nonlinear optimization problems.  
The solver is lightweight, readable, and designed for research prototyping.

---

## 1. Overview

`FMSQP.m` solves constrained optimization problems of the form:

\[
\begin{aligned}
\min_x\quad & f(x) \\
\text{s.t.}\quad & c_i(x)=0,\quad i=1,\dots,m_e, \\
& h_j(x)\le 0,\quad j=1,\dots,m_i.
\end{aligned}
\]

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

```matlab
[x, lambda, info] = FMSQP(fun, x0, Aeq, beq, Aineq, bineq, lb, ub, options)
````

### **Inputs**

| Argument         | Description                                     |
| ---------------- | ----------------------------------------------- |
| `fun`            | Function handle returning `[f, g]` at input `x` |
| `x0`             | Initial point                                   |
| `Aeq`, `beq`     | Equality constraints (`Aeq * x = beq`)          |
| `Aineq`, `bineq` | Inequality constraints (`Aineq * x ≤ bineq`)    |
| `lb`, `ub`       | Variable bound constraints                      |
| `options`        | Structure defining algorithm parameters         |

### **Function handle format**

Your function must follow:

```matlab
function [f, g] = fun(x)
```

where

* `f` = objective value
* `g` = gradient column vector

### **Outputs**

| Output   | Description                                                         |
| -------- | ------------------------------------------------------------------- |
| `x`      | Final solution vector                                               |
| `lambda` | Estimated Lagrange multipliers                                      |
| `info`   | Struct containing stopping status, iteration count, residuals, etc. |

---

## 4. Minimal Example

```matlab
% Objective function
fun = @(x) deal( x(1)^2 + x(2)^2, [2*x(1); 2*x(2)] );

% Initial point
x0 = [1; 1];

% No constraints
Aeq = []; beq = [];
Aineq = []; bineq = [];
lb = []; ub = [];

% Options
options.MaxIter = 100;
options.TolX = 1e-6;
options.TolFun = 1e-6;
options.Display = 'on';

% Run solver
[x, lambda, info] = FMSQP(fun, x0, Aeq, beq, Aineq, bineq, lb, ub, options);
```

---

## 5. Algorithm Summary

`FMSQP.m` follows a classical feasible-mode SQP framework:

1. Evaluate objective, gradient, and constraints
2. Linearize equality/inequality constraints
3. Form the QP subproblem:
   [
   \min_d \ \nabla f(x_k)^\top d + \tfrac12 d^\top B_k d
   ]
4. Solve for a search direction
5. Perform step-size selection (line search)
6. Update primal variables and multipliers
7. Check KKT residuals and feasibility
8. Stop when tolerances are satisfied

This formulation provides a clean platform for exploring new SQP techniques or verifying theoretical properties.

---

## 6. Options

Common fields in the `options` struct include:

```
options.MaxIter     % maximum number of iterations
options.TolFun      % optimality tolerance
options.TolX        % step-length tolerance
options.Display     % 'on' or 'off'
```

If omitted, default values inside `FMSQP.m` are applied.

---

## 7. Repository Structure

```
FMSQP/
├── FMSQP.m        % Main solver
├── MainG.m        % Driver example
├── Data_G.mat     % Sample data
├── Data_RS.mat    % Sample data
├── Ex*.m          % Example problems (Ex3.10, Ex5.1, etc.)
└── README.md
```

---

## 8. Notes

* Intended for **small or medium-scale** optimization problems
* For large-scale applications, consider:

  * sparse matrices
  * structure-exploiting BFGS updates
  * iterative solvers
* Algorithm performance depends on initial feasibility
* Recommended for research codebases and classroom demonstrations

---

## 9. License

Specify your desired license (e.g., MIT).
If omitted, all rights are reserved by default.

---

## 10. Contact

Author: **Wenhao Fu**
For questions or suggestions, please open an issue on GitHub.
