---
layout: post
title: "On sphere packing: Review"
author: "Nguyen Doan Dai"
categories: math
tags: math
---

$$
\newcommand{\vol}{\text{Vol}}
$$

(This is a summary of my study during the Spring semester of my first year at Ulm)

### Introduction

Let $$r > 0$$, denote $$B(p, r)$$ the ball of radius $$r > 0$$ centered at $$p$$. onsider the space $$\mathbb{R}^n$$, a sphere packing configuration is a set of points $$\mathcal{P} \subsetneq \mathbb{R}^n$$ such that for each point $$x, y \in \mathcal{P}$$, we have that $\|x - y\| \geq 2r$, meaning that if we assign to each point $$x \in \mathcal{P}$$ a ball centered at $$x$$ and of radius $$r$$, then the collection of balls $$\Omega = \{B(x, r) \mid x \in \mathcal{P}\}$$ satisfies that any two balls are either disjoint or intersecting at the boundary. We define the density of the sphere packing $$\mathcal{P}$$ as

$$
    \Delta(\mathcal{P}) = \sup_{p \in \mathbb{R}^n} \limsup_{r' \to \infty} \frac{\vol(\Omega \cap B(p, r'))}{\vol(B(p, r))}
$$


It can easily be shown that a packing attaining the supremum always exists, and moreover, Groemer {% cite Groemer1963 --file on_sphere_packing %} shows that there exists optimal packings for which the limit convergence is uniform with respect to $$p$$. Thus, it is reasonable to use the upper density as a characterisation of the global density of a packing.

In particular, we can define the density bound $$\Delta$$ as 

$$
    \Delta = \sup_{\mathcal{P}} \Delta(\mathcal{P})
$$


A special case is when all points in $$\mathcal{P}$$ lie in an organised manner, such as on a lattice $$\Lambda$$, which gives rise to the notion of a _lattice packing_. Combining some finite number of translation of a lattice packing whilst ensuring that no two points are of distance less than $$2r$$ from each other, we obtain a periodic packing. It is trivial that the density of a periodic packing is at least that of an optimal packing; it is nonetheless surprising that one can approach the optimal density by a periodic packing {% cite Cohn2003 --file on_sphere_packing %}.

Thus, we may concern ourselves with the analogous question of optimal density of a lattice packing. For a lattice $$\Lambda$$, its packing density is given by

$$
    \Delta(\Lambda) = \frac{\pi^{\frac{n}{2}}}{\Gamma\left(\frac{n}{2} + 1\right)} \left(\frac{r}{2}\right)^n \frac{1}{|\Lambda|}.
$$


One may ignore the two factors, which corresponds to the volume of $$B(0, r)$$ and the scaling factor, which now leaves us the centre density. For $$r = 1$$, this is precisely $\frac{1}{\|\Lambda\|}$.

When $$r = 1$$, is it clear that the shortest element of $$\Lambda$$ has length $$2$$. If it is less than $$2$$, some two balls intersect; if it is larger than $$2$$, one can reduce it to obtain a tighter packing. The question of minimising $\|\Lambda\|$ is equivalent, after a rescaling, to the problem of finding the shortest element of a unit-covolume lattice. In particular, if we consider $$\Lambda$$ of unit covolume, i.e. $\vol(\mathbb{R}^n / \Lambda) = \frac{1}{\|\Lambda\|} = 1$, denote $$\lambda_1(\Lambda)$$ the length of the shortest nonzero element of $$\Lambda$$. Then we define $$\gamma_n = \max_{\Lambda} \gamma_1(\Lambda)^2$$ the $$n^{\text{th}}$$ Hermite constant.

The link between Hermite constant and the optimal (centre) density of a lattice packing is made clear by the following relation.

$$
    \gamma_n^n = \frac{4}{|\Lambda|^2}
$$


### Linear programming bound

The case $$n = 2$$ was first solved by Fejes Toth, and in fact we have an alternative elementary solution {%cite Leppmeier2019 --file on_sphere_packing%}. Following Toth's approach, Hales announced a solution for $$n = 3$$ in 1998 {%cite Hales2005%}, formally checked in 2005 {%cite HALES2017 --file on_sphere_packing %}. This approach is generally deemed complicated and not generalisable to higher dimension.

In the early 2000s, Cohn and Elkies {% cite Cohn2003 --file on_sphere_packing %} introduced a new approach. The essential idea is that the term $\frac{1}{\|\Lambda\|}$ appears in the Poisson formula for higher dimension. The proof is simple enough that I will reproduce below, yet insightful enough to be included before we discuss possible progress.

**Theorem 1.** (Theorem 3.1, {% cite Cohn2003 --file on_sphere_packing %}) Suppose $$f : \mathbb{R}^n \to \mathbb{R}$$ is a continuous and integrable function, is not identically zero, and satisfies the following conditions:
1. $$f(x) \leq 0$$ for $\|x\| \geq 1$, and
2. $$\widehat{f}(t) \geq 0$$ for all $$t$$.

Then, the centre density of $$n$$-dimensional sphere packings is bounded from above by

$$
    \frac{f(0)}{2^n \widehat{f}}.
$$

_Proof_. Consider a periodic packing given by $$N$$ translations of a lattice $$\Lambda$$ by vectors $$v_1, v_2, \cdots, v_N$$. If we choose the scaling so that no points are closer than $$1$$, then we can pack with spheres of radius $$\frac{1}{2}$$, and the centre density is given by $$\frac{N}{2^n \Lambda}$$.

Consider such a function $$f$$, Poisson formula gives

$$
    \sum_{x \in \Lambda} f(x + v) = \frac{1}{|\Lambda|} \sum_{t \in \Lambda*} \widehat{f}(t) \left| \sum_{1 \leq j \leq N} e^{2\pi i \langle v_j, t\rangle} \right|.
$$

All terms are non-negative, so the sum is bounded by the term corresponding to $$t = 0$$, namely $\frac{N^2\widehat{f}(0)}{\|\Lambda\|}$. On the other hand, for the left hand side, $$x + v_j - v_k$$ is the difference between two centers in the packing, thus $\|x + v_j - v_k\| \geq 1$ if and only if $$j \neq k$$. Whenever it happens, the terms is then non-positive, so the sum may be bounded from above by $$N f(0)$$, corresponding to all the terms when $$j = k$$. Thus

$$
    N f(0) \geq \frac{N^2 \widehat{f}(0)}{|\Lambda|},
$$

and rearranging the inequality gives our result. <span style="float:right;">$$\square$$</span>

Hunting for good functions $$f$$ can then be seen as a linear programming problem (with infinitely many constraints) and it turned out to be quite difficult. One can see that the theorem remains even if we "rotate" $$f$$, thus by taking average over the rotation, we may assume that $$f$$ is radial. In their paper, using the same technique for the case of isodual lattices (which contains the lattices giving the best density for dimension 8 and 24 at the time, namely $$E_8$$ and Leech), Cohn and Elkies uses $g_k(x) = L_k^\alpha(2\pi\|x\|^2) e^{-\pi\|x\|^2}$ where $$L^a_k$$ denotes Laguerre polynomial with respect to the measure $$e^{-x} x^\alpha dx$$ on $$\mathbb{R}_{\geq 0}$$. These functions form a basis for the radial eigenfunctions of the Fourier transform with eigenvalues $$(-1)^k$$, and the problem now boils down to finding a good linear combination of these functions subjected to the constraints in the theorem above (or, in this case, some analogs thereof).

With this approach, they got an uppoer bound very close to the density given by $$E_8$$ and Leech lattice, suggesting that with a suitable choice of $$f$$, we might solve the packing problem for dimension $$8$$ and $$24$$. That "suitable choice of $$f$$" turns out to be much more difficult than a numerical problem. As Cohn and Elkies noted, such a function $$f$$ must have its zeros _and_ zeros of its Fourier transform located at $$\sqrt{2n}$$ for integers $$n \geq 1$$ (note that since $$f$$ - and thus $$\widehat{f}$$ - is radial, we can consider $$f : \mathbb{R}_{\geq 0} \to \mathbb{R}$$), and it was unknown at the time how such a task could be done.

In 2016, Viazovska gave such a construction using theory of modular forms for dimension $$8$$, and further collaboration lead to analogous solution for dimension 24. Work carried out after this breakthrough finally culminated into what is known as Fourier-Viazovksa-Radchenko interpolation {% cite Radchenko2019 --file on_sphere_packing %}, which allows to prove that such a construction is in fact unique.

### Semidefinite programming and three-point bound

The next question is then if this approach can be used for other dimension and/or if there are variants of this approach that allow us to do so.

For the first question, it is in fact conjectured that the answer is no, i.e. the Cohn-Elkies linear programming bound is conjectured not to be sharped in any other dimension larger than $$2$$, and it was shown to be the case for $$n \in \{3,4, 5, 6, 12, 16\}$$. To do so, we need two ingredients: a method to bound from below the optimal bound that can be achieved by linear programming method, say, by $$A$$, and a method to bound from above the optimal packing density, say, by $$B$$. If $$A > B$$, then we know that the linear programming approach is not sharp.

For the first part, as the name implies, this linear programming problem has a dual as a maximisation problem, whose optimal value gives a lower bound on that of the primal. If this program were finite-dimensional, the two optimal values would coincide, but as it is infinite-dimensional, such strong duality is more subtle. It was conjectured by Cohn {% cite Cohn2002 --file on_sphere_packing %}, and proven by Cohn, de Laat, and Salmon {% cite Cohn2022 --file on_sphere_packing %}. If the primal is already hard, the dual is even harder, as now we are not looking for optimal functions, but optimal measures, where such a measure should be supported on a discrete set of radii, and thus be singular. There is no known family of measures analogous to the family of Gaussian used by Cohn and Elkies above.

Nevertheless, Cohn and Triantafillou {% cite Cohn2021 --file on_sphere_packing %} optimised over the space of modular forms to find feasible points for the dual program, thus established lower bounds for dimension $$n = 12$$ and $$n = 16$$. They proposed that this approach should be generalisable to other dimensions; the reason why their result is limited to dimension divisible by 4 is their usage of only modular forms of even weight to generate possible measures (specifically, the modular forms for $$\Gamma_0(N)$$, which has a nice basis of modular forms with rational coefficients in $$q$$-expansion). By using $$\Gamma_1(N)$$ and forms of odd weights, in theory, one can carry out the computations for even dimensions. This is precisely what was done by de Courcy-Ireland, Dostert, and Viazovska {% cite DeCourcy-Ireland2024 --file on_sphere_packing %}: using modular forms of weight $$3$$ for $$\Gamma_0(48)$$, they achieved a good enough bound for $$n = 6$$.

In a different direction, Li {% cite Li2022 --file on_sphere_packing %} used discrete reduction to obtain finite-dimensional linear programming programs, which a priori give weaker bounds than solving the dual of Cohn-Elkies linear programming problem, but such bounds turned to be good enough for $$n \in \{3, 4, 5\}$$.

For the second part, recall that in the proof of the theorem by Cohn and Elkies, we used the argument that if $$j \neq k$$, $\|x + v_j - v_k\| \geq 1$, thus the corresponding term $$f(x + v_j - v_k)$$ is non-negative, and we bounded the sum by $$Nf(0)$$. In essence, we are comparing centers of two spheres, and this is generalisable to $$n$$ centers, giving rise to what is known as $$n$$-point bound. Consider the case $$n = 3$$, as studied by Cohn, de Laat, and Salmon, an analog of Theorem 1 is the following.

**Theorem 2.** (Theorem 1.1, {% cite Cohn2022 --file on_sphere_packing %}) Let $$0 < r < R$$, and suppose $$f_2 : \mathbb{R}^n \to \mathbb{R}$$ and $$f_3 : \mathbb{R}^n \times \mathbb{R}^n \to \mathbb{R}$$ continuous and integrable such that $$\widehat{f_2}(0) = 1$$, $$\widehat{f_2} \geq 0$$ and $$f_3$$ is positive definite as a kernel. If

$$
    f_2(x) + f_3(0, 0) + f_3(0, x) + f_3(x, x) \geq 0 \text{ for } r \leq x \leq R,  
$$


$$
    f_2(x) \leq 0 \text{ for } | x | \geq R,
$$

and

$$
    f_3(x, y) \leq 0 \text{ whenever } r \leq |x|, |y| \leq R \text{ and } |x - y | \geq r,
$$

then the centre sphere packing density in $$\mathbb{R}^n$$ is at most $$\frac{f_2(0) + f_3(0, 0)}{2^n}$$.

Using this theorem, and considering essentially linear combinations of some polynomials times the exponential, they arrived at upper bounds which are lower than the upper bound of the optimal value given by the linear programming approach, which suffices to prove that the linear programming approach is not optimal for the dimensions mentioned above.

The proof of the theorem goes similarly as above, though it is not clear how to generalise this theorem to bigger number of points. I find that it is easier to see the theorem for lattice packing.

**Theorem 3.** (Theorem 1.4, {% cite Cohn2022 --file on_sphere_packing %}) Let $$r > 0$$, and define

$$
    S_n = \{ (x, y) \in \mathbb{R}^n \times \mathbb{R}^n : \|x\|, \|y\|, \|x + y\|, \|x-y\| \in \{0\} \cup [r, \infty)\}.
$$

Suppose $$f : \mathbb{R}^n \times \mathbb{R}^n \to \mathbb{R}$$ is a continuous, integrable function with $$\widehat{f} \geq 0$$, $$\widehat{f}(0, 0) = 1$$, and

$$
    f(x, y) \leq 0 \text{ for } (x, y) \in S \setminus \{(0, 0)\}.  
$$

Then the centre density of optimal lattice packing in $$\mathbb{R}^n$$ is at most $$\sqrt{f(0, 0)}$$.

The proof essentially goes as for Theorem 1, and moreover it is directly generalisable to the following, as noted by the authors.

**Theorem 4.** Let $$r > 0$$, and define

$$
    S_n = \left\{ x \in \left(\mathbb{R}^{n}\right)^d : \forall (w_i)_{i = 1}^d \ni gcd(w_1, \cdots, w_d) = 1, \left|\left| \sum_i w_i x_i\right|\right| \in \{0\} \cup [r, \infty)\right\}.
$$

Suppose $$f : \left(\mathbb{R}^n\right)^d \to \mathbb{R}$$ is a continuous, integrable function with $$\widehat{f} \geq 0$$, $$\widehat{f}(0, \cdots, 0) = 1$$, and

$$
    f(x) \leq 0 \text{ for } x \in S \setminus \{(0, \cdots, 0)\}.  
$$

Then the centre density of optimal lattice packing in $$\mathbb{R}^n$$ is at most $$\sqrt[d]{f(0, \cdots, 0)}$$.

### References

{% bibliography --file on_sphere_packing %}