# FMSQP.m — A MATLAB Filter-based Multi‑start SQP Solver for Inequality‑Only Problems

This repository contains **FMSQP.m**, a MATLAB implementation of a **Sequential Quadratic Programming (SQP)** solver tailored for nonlinear optimization with **inequality constraints only**. The method supports **multiple starting points** and uses a **filter technique** to improve global convergence, especially to escape infeasible stationary points or poor local minima.

---

## 1. Overview

`FMSQP.m` solves optimization problems of the form:

$$
\begin{aligned}
\min_x\quad & f(x) \\
\text{s.t.}\quad & c_i(x) \le 0,\ i = 1,\dots,m.
\end{aligned}
$$

Key properties:

- It does **not require intermediate iterates to be feasible**.  
- If the problem is infeasible, the algorithm aims to find a point that **minimizes violation** of the inequality constraints while also trying to reduce the objective.  
- By running from **multiple initial points**, it increases the chance of finding a good compromise or a feasible solution.  
- A **filter strategy** is employed: iterates that are dominated in terms of constraint violation and objective value are rejected, helping to avoid convergence to bad infeasible points or poor local minima.

---

## 2. Key Features

- **Single-file implementation** in MATLAB: `FMSQP.m`  
- Supports:
  - Nonlinear objective  
  - Nonlinear inequality constraints  
- **Multi-start**: you can supply multiple initial guesses to explore different basins.  
- **Filter acceptance mechanism**: no penalty function needed — filter handles trade-off between feasibility and optimality.  
- Suitable for research, algorithm development, and numerical experiments.

---

## 3. Function Interface

Call format:

```matlab
[x, fx, output] = FMSQP(X0, opts);
````

* `X0` — initial points, an \(n \times J\) matrix where each column is a different starting point (multi-start).
* `opts` — options struct with fields:

  * `opts.varbose` (verbosity, default `0`)
  * `opts.epsilon` (tolerance, default `1e-5`)
  * `opts.nmax` (maximum number of iterations, default `500`)

Output:

* `x` — final point (from one of the runs) that passed the filter or best compromise.
* `fx` — objective value at `x`.
* `output` — struct with diagnostic info: number of iterations, number of function / gradient evaluations, final constraint violations, and time.

---

## 4. Example Usage

```matlab
% Define multiple initial points (each column is one point)
X0 = [0, 1, 0, 1;
      0, 0, 1, 1];

opts.varbose = 1;
opts.epsilon = 1e-6;
opts.nmax = 300;

[x, fx, output] = FMSQP(X0, opts);

fprintf('Found x = [%g, %g]\n', x(1), x(2));
fprintf('f(x) = %g\n', fx);
fprintf('Constraint violation: %g\n', output.vx);
```

---

## 5. Algorithmic Overview

1. From **each initial point**, run an SQP iteration.
2. At every iteration:

   * Evaluate objective ( f(x) ) and its gradient.
   * Evaluate constraint violations ( c(x) ) and compute their gradients if needed.
   * Linearize the constraints around the current point.
   * Build and solve a feasible quadratic programming subproblem.
   * Use a **line search** to get a trial step.
   * Use a **filter** to decide whether to accept the trial point:

     * A point is accepted if it is not dominated in both objective value and constraint violation by any point in the filter.
     * Otherwise, it is rejected, which helps the method avoid bad or infeasible regions.
3. Continue until convergence criteria (optimality and violation tolerances) are met or maximum iterations reached.
4. Among all runs (from different starts), pick the best accepted point according to the filter.

This **multi-start + filter** strategy provides good robustness: the method can escape infeasible stationary points and poor local minima.

---

## 6. Required User-Defined Functions

You should provide these four MATLAB files in the same directory / on path:

* `funf.m` — computes \( f(x) \)
* `gradf.m` — computes \( \nabla f(x) \)
* `func.m` — computes constraint vector \( c(x) \le 0 \)
* `gradc.m` — computes Jacobian \( \nabla c(x) \) (each column is gradient of one constraint)

---

## 7. Output Structure

The `output` struct returned by `FMSQP` includes:

* `output.m` — number of constraints
* `output.vx` — final constraint violation values
* `output.nit` — number of iterations
* `output.nf` — number of function evaluations
* `output.ng` — number of gradient evaluations
* `output.time` — elapsed time (in seconds)

---

## 8. License

This project is licensed under the **MIT License**.
See the [LICENSE](LICENSE) file for full terms.

---

## 9. Contact / Contribution

* **Author**: Wenhao Fu
* **E-mail**: wenhfu@usts.edu.cn
* Contributions welcome! Feel free to open issues or pull requests to improve the multi-start logic, filter strategy, or performance.

