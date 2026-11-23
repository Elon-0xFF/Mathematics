# Fixed-Point Iteration: From Existence to Linearization and Chaos

> This document systematically organizes the core theory around **fixed points** in first-order iteration $x_{n+1}=f(x_n)$:
> how to decide whether fixed points exist, how to construct a transformation that linearizes the iteration,
> and what role fixed points play in chaotic systems.
> An advanced reading list (all in English) is included at the end for further study.

---

## Contents

1. [Fundamentals of Fixed Points](#1-fundamentals-of-fixed-points)
2. [Existence Criteria for Fixed Points](#2-existence-criteria-for-fixed-points)
3. [Why Fixed-Point Methods Simplify Iteration](#3-why-fixed-point-methods-simplify-iteration)
4. [Iterations That Can Be Linearized](#4-iterations-that-can-be-linearized)
5. [Systematically Finding a Linearizing Transformation g(x)](#5-systematically-finding-a-linearizing-transformation-gx)
6. [Fixed Points in Chaotic Systems](#6-fixed-points-in-chaotic-systems)
7. [Failure Modes and Limitations](#7-failure-modes-and-limitations)
8. [Quick Reference Table](#8-quick-reference-table)
9. [Suggested Further Reading (in English)](#9-suggested-further-reading-in-english)

---

## 1. Fundamentals of Fixed Points

**Definition**: A point $p$ is called a **fixed point** of $f$ if $f(p)=p$.

For the iteration
$$x_{n+1}=f(x_n),$$
fixed points are the candidate values where the iteration may "stay put" or "converge to".

### The Full General Term in the Linear Case

For $f(x)=ax+b$, write the $n$-th iterate as $f^{(n)}(x)$:

- **When $a=1$**: $f^{(n)}(x)=x+nb$ (arithmetic progression)
- **When $a\neq 1$**: the fixed point is $p=\dfrac{b}{1-a}$, and
$$\boxed{\,f^{(n)}(x)=a^n x+\frac{b(1-a^n)}{1-a}=a^n\bigl(x-p\bigr)+p\,}$$

> **Memory aid**: fixed point $p=\dfrac{b}{1-a}$; the iteration is just $a^n(x-p)+p$.

---

## 2. Existence Criteria for Fixed Points

### 2.1 Continuous Self-Maps (Intermediate Value Theorem)

If $f:[a,b]\to[a,b]$ is continuous, then there exists $p\in[a,b]$ with $f(p)=p$.

Proof idea: set $h(x)=f(x)-x$; then $h(a)\ge0$ and $h(b)\le0$, so by the Intermediate Value Theorem $h$ has a zero.

More generally, if there exist $x_0,x_1$ with $f(x_0)\ge x_0$ and $f(x_1)\le x_1$, then a continuous $f$ must have a fixed point between them.

### 2.2 Contraction Mapping (Banach Fixed-Point Theorem)

On a complete metric space $(X,d)$, if there exists $0\le L<1$ such that
$$d(f(x),f(y))\le L\,d(x,y),\quad\forall x,y\in X,$$
then $f$ has a **unique** fixed point, and the iteration converges to it from *any* starting value.

This is the theoretical cornerstone of numerical methods such as Newton's method, Jacobi iteration, and Gauss–Seidel iteration.

### 2.3 Brouwer Fixed-Point Theorem (Higher-Dimensional Generalization)

If $f:B^n\to B^n$ is continuous ($B^n$ the closed $n$-ball), then a fixed point exists.
This generalizes the Intermediate Value Theorem to higher dimensions, but it **does not guarantee uniqueness** and does not construct the fixed point.

### 2.4 Practical Decision Flow

1. Is the domain compact, and does $f$ map the domain into itself? → Use the IVT / Brouwer.
2. Is $f$ a contraction ($|f'|<1$ holding uniformly)? → Use Banach: get uniqueness + convergence.
3. Otherwise, solve $f(x)-x=0$ numerically (bisection, Newton's method).

---

## 3. Why Fixed-Point Methods Simplify Iteration

The core idea: **subtracting the fixed point often turns the error into a geometric progression.**

Let the error be $e_n=x_n-p$, where $f(p)=p$:

- **Linear** $f(x)=ax+b$:
$$e_{n+1}=a x_n+b-p=a(x_n-p)+\underbrace{(ap+b-p)}_{=0}=a\,e_n$$
The error is exactly geometric, and the general term follows immediately.

- **Nonlinear but differentiable at $p$**: the first-order approximation
$$e_{n+1}\approx f'(p)\,e_n$$
is precisely **local linearization**. Writing $\lambda=f'(p)$:
  - $|\lambda|<1$: the error shrinks and the iteration converges (**attracting fixed point**);
  - $|\lambda|>1$: the error grows and the iteration diverges (**repelling fixed point**);
  - $|\lambda|=1$: higher-order analysis is needed.

---

## 4. Iterations That Can Be Linearized

### 4.1 Linear Recurrence $x_{n+1}=ax_n+b$
Construct $p=b/(1-a)$ directly; the error is geometric — see Section 1.

### 4.2 Linear-Fractional Recurrence $x_{n+1}=\dfrac{ax_n+b}{cx_n+d}$

The fixed points satisfy the quadratic equation $cp^2+(d-a)p-b=0$. If there are two distinct fixed points $p_1,p_2$, set
$$y_n=\frac{x_n-p_1}{x_n-p_2},$$
then $y_{n+1}=k\,y_n$ (geometric), where $k=\dfrac{a-cp_1}{a-cp_2}$.
Binet's formula for the Fibonacci numbers is essentially this idea.

### 4.3 Nonlinear Recurrences Reducible to Linear Ones

| Original iteration | Substitution | Becomes |
|---|---|---|
| $x_{n+1}=a x_n^b$ | $y=\ln x$ | $y_{n+1}=b\,y_n+\ln a$ (linear) |
| $x_{n+1}=k\sqrt{x_n}$ | $y=\ln x$ | $y_{n+1}=\tfrac12 y_n+\ln k$ |
| $x_{n+1}=\dfrac{x_n}{1+x_n}$ | $y=1/x$ | $y_{n+1}=y_n+1$ (arithmetic) |
| $x_{n+1}=2x_n(1-x_n)$ | $x_n=\tfrac12(1-y_n)$ | $y_{n+1}=y_n^2$, then take logarithms to linearize |
| $x_{n+1}=4x_n(1-x_n)$ | $x=\sin^2\theta$ | $\theta_{n+1}=2\theta_n$ (angle doubling) |

### 4.4 Newton's Method

$x_{n+1}=x_n-\dfrac{f(x_n)}{f'(x_n)}$ is a fixed-point iteration with $g(x)=x-f(x)/f'(x)$,
whose fixed points are exactly the roots of $f(x)=0$. Near a simple root, $g'(x^*)=0$, hence **local quadratic convergence**.

---

## 5. Systematically Finding a Linearizing Transformation g(x)

The core idea: **conjugacy**. Find an invertible transformation $\varphi$ such that
$$\varphi\circ f\circ\varphi^{-1}=L,$$
where $L$ is a linear function. Then the original iteration becomes $y_{n+1}=L(y_n)$, and the general term follows immediately.

### 5.1 Local: Schröder's Equation

Look for $\varphi$ satisfying
$$\varphi(f(x))=\lambda\,\varphi(x),\qquad \lambda=f'(p).$$
Once found, $\varphi(x_n)=\lambda^n\varphi(x_0)$, so
$$x_n=\varphi^{-1}\bigl(\lambda^n\varphi(x_0)\bigr).$$

**Constructing $\varphi$ — the Kœnigs power-series method**:
Assume $\varphi$ is analytic near $p$, $\varphi(x)=\sum_{k\ge1}c_k(x-p)^k$; substitute into the functional equation and compare coefficients term by term to determine the $c_k$ recursively.

### 5.2 Global: Common Substitution Patterns

1. **Linear-fractional** → cross-ratio transformation $\varphi(x)=\dfrac{x-p_1}{x-p_2}$.
2. Involving **powers / exponentials** → take logarithms.
3. Involving a **reciprocal-fraction structure** → take reciprocals.
4. Nested **square roots** → square both sides, or use a trigonometric substitution.
5. **Symmetric structure** $f(1-x)=1-f(x)$ → try $\varphi(x)=x-\tfrac12$ or $\ln\dfrac{x}{1-x}$.

### 5.3 Systematic Trial Steps

1. Find the fixed point $p$;
2. Compute $\lambda=f'(p)$;
3. If $\lambda\neq0,1$, try the Schröder power series;
4. Recognize the algebraic structure (fractional / power / trigonometric) and apply the corresponding substitution;
5. Use symmetry to guess $\varphi$;
6. Verify that $\varphi(f(x))=L(\varphi(x))$ closes into a linear form.

> **There is no universal algorithm**, but for the standard types of problems the transformation can always be found by "observing symmetries + trial and error".

---

## 6. Fixed Points in Chaotic Systems

### 6.1 What Chaos Is

**Apparently random, long-term unpredictable** behavior arising in deterministic nonlinear systems, generated by the internal dynamics rather than by external noise.

### 6.2 Devaney's Definition of Chaos

For continuous $f:X\to X$, $f$ is called chaotic if it simultaneously satisfies:

1. **Sensitive dependence on initial conditions**: there exists $\delta>0$ such that in every neighborhood of every point there is a point whose orbit eventually separates from it by more than $\delta$;
2. **Topological transitivity**: for any two nonempty open sets $U,V$, there exists $n$ with $f^{(n)}(U)\cap V\neq\varnothing$ (the system cannot be decomposed into two invariant subsets);
3. **Density of periodic points**: the periodic points are dense in $X$.

### 6.3 Classic Example: the Logistic Map $x_{n+1}=rx_n(1-x_n)$

| Parameter $r$ | Behavior |
|---|---|
| $[0,1]$ | converges to 0 |
| $(1,3)$ | converges to the positive fixed point $1-1/r$ |
| $(3,\,3.5699\ldots)$ | period-doubling bifurcations (2→4→8→…) |
| $\approx 3.5699$ | Feigenbaum point: onset of chaos |
| $=4$ | fully chaotic, $\lambda=\ln2>0$, solution $x_n=\sin^2(2^n\pi\theta)$ |

### 6.4 How Chaos Relates to Fixed Points

- Chaotic systems typically have **infinitely many periodic orbits** (fixed points are period-1 orbits), but they are **all unstable** (repelling, $\lvert f'\rvert>1$);
- Example: at $r=4$ both fixed points $x=0$ and $x=3/4$ are repelling;
- Chaotic orbits do not converge to any fixed point; instead they wander ergodically across the interval.

### 6.5 Tools for Detecting Chaos

- **Lyapunov exponent** $\displaystyle\lambda=\lim_{n\to\infty}\frac1n\sum_{i=0}^{n-1}\ln|f'(x_i)|$: $\lambda>0$ indicates chaos;
- Continuous power spectrum (in contrast to the discrete spectrum of periodic motion);
- Positive topological entropy;
- Numerical experiments: the rate at which nearby trajectories separate.

---

## 7. Failure Modes and Limitations

| Situation | Problem |
|---|---|
| No fixed point (e.g. $f(x)=x+1$) | the iteration diverges; there is nothing to build on |
| Multiple fixed points | different initial values converge to different fixed points; basins of attraction must be analyzed |
| Periodic orbit but no fixed point (e.g. $f(x)=-x$) | the fixed-point method yields no general term |
| $\lvert\lambda\rvert=1$ | linearization fails; higher-order terms are needed |
| Chaotic regime | all fixed points are unstable; no simple closed-form general term |

---

## 8. Quick Reference Table

| Question | Key tool |
|---|---|
| Does a fixed point exist? | Intermediate Value Theorem / Banach contraction / Brouwer |
| Unique? Convergent? | Contraction property $\lvert f'\rvert<1$ |
| Rate of convergence | $\lambda=f'(p)$: if $0$, quadratic (Newton); otherwise linear |
| How to linearize | Schröder's equation + power series; cross-ratio for fractional maps; substitute according to algebraic structure |
| Is it chaotic? | Lyapunov exponent $>0$; Devaney's three conditions |

---

## 9. Suggested Further Reading (in English)

The books below are ordered by increasing difficulty and cover fixed-point theory, discrete dynamical systems, and chaos:

1. **Steven H. Strogatz — *Nonlinear Dynamics and Chaos***
   The best introductory book on nonlinear dynamics and chaos: intuitive and richly illustrated. Read it first to build intuition; it covers fixed points, bifurcations, Lyapunov exponents, and chaos.

2. **Robert L. Devaney — *An Introduction to Chaotic Dynamical Systems***
   A rigorous chaos textbook starting from one-dimensional maps; Devaney's definition of chaos comes out of this line of work. Symbolic dynamics and the Sharkovskii theorem are explained clearly.

3. **John Milnor — *Dynamics in One Complex Variable***
   The classic on complex dynamics (Julia sets, Fatou sets, fixed points and local linearization); requires some background in complex analysis.

4. **Walter Rudin — *Principles of Mathematical Analysis* (Ch. 9, contraction mappings)**
   The standard rigorous source for the Banach fixed-point theorem and the contraction mapping principle.

5. **V. I. Arnold — *Ordinary Differential Equations***
   Treats fixed points (equilibria), linearization, stability, and the Poincaré map from the viewpoint of differential equations; extremely strong geometric intuition.

6. **R. Clark Robinson — *An Introduction to Dynamical Systems: Continuous and Discrete***
   Covers continuous and discrete dynamical systems together — fixed points, stability, bifurcations, and chaos treated in a unified way; friendly exercises.

7. **Earl A. Coddington & Norman Levinson — *Theory of Ordinary Differential Equations***
   A classic ODE theory reference; a rigorous foundation for existence–uniqueness and stability theory.

8. **Predrag Cvitanović et al. — *Chaos: Classical and Quantum***
   A free online book (chaosbook.org); starts from periodic-orbit theory and goes deep into the quantitative analysis of chaos.

> **Suggested reading order**: Strogatz (intuition) → Devaney (rigor for discrete chaos) → Robinson (unified continuous + discrete) → Arnold / Milnor (advanced geometry and complex dynamics).

---

*This document was compiled from a series of earlier discussions on fixed-point iteration, for systematic review and further exploration.*
