---
layout: post
title: "On sphere packing: Discrete reductions"
author: "Nguyen Doan Dai"
categories: math
tags: math
---

### Li's discrete reductions

Let us focus on Li's proof for dimension $$n \in \{3, 4, 5\}$$ {% cite Li2022 --file on_sphere_packing %}. A function $$f : \mathbb{R}^n \to \mathbb{R}$$ is called an auxiliary function if it is smooth, and satisfies the following conditions.
- $f(x) = O\left(\frac{1}{(1 + \|x\|)^k}\right)$ for all $$k$$, 
- $$f(0) > 0$$, 
- $$\widehat{f}(0) \geq 1$$,
- $$f(x) < 0$$ for all $$\|x\| \geq r$$,
- $$\widehat{f} \geq 0$$.

We can analogously define the notion of auxiliary functions over $$\mathbb{Z}$$ and $$\mathbb{Z}_m = \mathbb{Z}/m\mathbb{Z}$$, with the Fourier transform now replaced by its counterparts on the corresponding groups. Here, in case of $$\mathbb{Z}_m^n$$, we use the coset representation $$(-m/2, m/2]$$ for $$\mathbb{Z}_m$$, and we normalise the Fourier transform by 

$$
    \widehat{f}(y) = m^{-\frac{n}{2}} \sum_{x \in \mathbb{Z}_m^n} f(x) e^{-2\pi \frac{i}{m} \langle x, y \rangle},
$$

so the analog of "$$\widehat{f}(0) \geq 1$$" for $$\mathbb{Z}_m^n$$ is $$\widehat{f}(0) \geq m^{-\frac{n}{2}}$$.

The crucial result in the paper by Li is the following theorem, whose proof is included for completeness. Recall that such an auxiliary function $$f$$ over $$\mathbb{R}^n$$ verifies the hypothesis of theorem by Cohn and Elkies, thus giving an upper bound of $$\left(\frac{r}{2}\right)^d$$.

**Theorem 1.** (Theorem 3.2, {% cite Li2022 --file on_sphere_packing %}) For any $$r > 0$$ and positive $$n$$ and $$m \geq 2r$$, if there exists an auxiliary function $$f$$ on $$\mathbb{R}^n$$, then there exists an auxiliary function $$g_m$$ on $$\mathbb{Z}_m^n$$ with the same $$r$$ and $$m^{-\frac{n}{2}} \frac{g_m(0)}{\widehat{g_m}(0)} \leq \frac{f(0)}{\widehat{f}(0)}$$.

_Proof_. Consider an auxiliary function $$f$$ on $$\mathbb{R}^d$$. First, we discretise $$f$$ by defining 

$$
    \widehat{g}(x) = \sum_{y \in \mathbb{Z}^n} \widehat{f}(x + y).
$$

Then $$\widehat{g}$$ is well-defined on the torus $$\mathbb{T}^n$$ as both $$f$$ and $$\widehat{f}$$ are rapidly decreasing. Then Fourier's inverse transform gives for all $$x \in \mathbb{Z}^d$$,

$$
    g(x) = \int_{\mathbb{T}^n} e^{2\pi i \langle x, y \rangle} \widehat{g}(y) dy = \int_{\mathbb{R}^n} e^{2\pi i \langle x, y \rangle} \widehat{f}(y) dy = f(x),
$$

or in other words, $g = f \Big\|\_\{\mathbb{Z}^n\}$. Now, define

$$
    g_m(x) = \sum_{y \in m\mathbb{Z}^n} g(x + y),
$$

we have that since $$g \in L^1(\mathbb{Z}^n)$$ and $$\sum_{x \in \mathbb{Z}_m^n} \|g_m(x)\| \leq \sum_{x \in \mathbb{Z}^n} \|g(x)\|$$, $$g_m \in L^1(\mathbb{Z}_m^n)$$.

Then the Fourier transform of $$g_m$$ is given by

$$
    \widehat{g_m}(y) = m^{-\frac{n}{2}} \sum_{x \in \mathbb{Z}_m^n} e^{-2\pi \frac{i}{m} \langle x, y \rangle} g_m(x) = m^{-\frac{n}{2}} \sum_{x \in \mathbb{Z}^n} e^{-2\pi i \langle x, y / m \rangle} g(x) = m^{-\frac{n}{2}} \widehat{g}(y/m)
$$

for all $$y \in \mathbb{Z}_m^n$$. Thus, 

$$
    \widehat{g_m}(0) = m^{-\frac{n}{2}} \widehat{g}(0) \geq m^{-\frac{n}{2}} \widehat{f}(0) \geq m^{-\frac{n}{2}}, 
$$

and $$\widehat{g_m} \geq 0$$. On the other hand, using the coset representation $$(-m/2, m/2]$$ for $$\mathbb{Z}_m$$ with $$\|x\|$$ for $$x \in \mathbb{Z}^n_m$$ as embedded in $$\mathbb{R}^n$$. As $$m \geq 2r$$, if $$\|x\| \geq r$$ in $$\mathbb{Z}_m^n$$, then $$\|x + y\| \geq r$$ in $$\mathbb{Z}^n$$, implying that for $$\|x\| \geq r$$ in $$\mathbb{Z}_m^n$$, we have

$$
    g_m(x) = \sum_{y \in m\mathbb{Z}_m^n} g(x + y) \geq 0.
$$


As $$\| y \| \geq m > r$$, for all $$y \in m\mathbb{Z}^n \setminus \{0\}$$, 

$$
    g_m(0) = \sum_{y \in m\mathbb{Z}^n} g(y) \leq g(0) = f(0),
$$

as desired. <span style="float:right;">$$\square$$</span>

Now finding the optimal $$g_m$$ is a finite-dimensional linear programming problem (albeit of size potentially very large).

### Why $$\mathbb{Z}^n$$ and $$\mathbb{Z}_m^n$$? Why not something else?

Intuitively, we are sampling the plane at the integral point, but the key constraint is $$f(x) < 0$$ for all $$\|x\| \geq r$$, which is a condition on the ball $$B(0, r)$$, and we are approximating the ball by hypercubes. At the limit - whatever that means at the moment - as hypercubes get smaller and smaller, the approximation gets finer and finer, and it does not matter. Nevertheless, Li's proof for $$n \in \{3, 4, 5\}$$ did not require $$m$$ to be too large, e.g. $$m = 30$$ for $$n = 4$$, or $$m = 24$$ for $$n = 5$$, and the size of the liner programming problem quickly becomes an obstacle as $$n$$ grows. It is thus interesting to improve the method, say, by a better approximation.

At the same time, one still wishes that there is some nice structure to the set of sampling points so that much, if not all, of the proof above gets carried over easily. When $$n = 2$$, one of such candidates is the honeycomb lattice: the intuition goes that an inscribed hexagon approximates a circle better than an inscribed square. But it is not immediate how to generalise that to higher dimension.

But I abandoned this path for now. Many reasons motivated me to consider regular polyhedra, e.g. the radial nature of $$f$$ motivates that of $$g$$ and eventually $$g_m$$ but even in three dimensional, only cube out of all the Platonic solids tiles the space, since the other solids have nonzero Dehn invariant. The same problem persists for higher dimensions. In theory, one could consider other polytopes, e.g. truncated octahedron, which also tiles the space, but the calculation quickly proved to be cumbersome whilst not hinting any immediate benefits.

The same reasons justified abandonning the consideration of more general lattices instead of $$\mathbb{Z}^n$$. The lattices being conjectured to be optimal have high degrees of symmetry, and it makes little sense (whilst still possible!) that one would require a highly asymetric structure to prove that linear programming bound is ineffective.

### Convergence.

One of the questions left open is if the bound given by $$g_m$$, called discrete Cohn-Elkies linear programming bound, converges to that given by Cohn-Elkies original linear programming problem as $$m$$ and $$r$$ tends to infinity. This question has two parts.
1. Given $$f$$, does $$g_m$$ give the bound given by $$f$$ as $$m$$ tends to infinity?
2. Does the bound given by an optimal $$g_m$$ converge to the bound given by an optimal $$f$$ as $$m$$ and $$r$$ tends to infinity?

The first question has positive answer, provided that $$m$$ and $$r$$ approach infinity in a particular way. In short, in his note {% cite Cohn2002 --file on_sphere_packing %}, Cohn wrote that "rapidly decreasing" could be replaced by "Schwartz", and it turns out that this is the case (Proposition 3.5, {% cite Cohn2022 --file on_sphere_packing %}). If $$f$$ is a Schwartz function, so is $$\widehat{f}$$, and in particular we know that both of them are rapidly decreasing.

Essentially we need to control two differences, $$f(0) - g_m(0)$$ and $$\widehat{g}(0) - \widehat{f}(0)$$, and to prove that $$m^{-\frac{n}{2}} \frac{g_m(0)}{\widehat{g_m}(0)}$$ converges to $$\frac{f(0)}{\widehat{f}(0)}$$, it suffices to prove that these two differences tends to $$0$$, which is trivial by faster-than-polynomial decay rate of both $$f$$ and $$\widehat{f}$$. For example, 


$$
    f(0) - g_m(0) = \sum_{y \in m\mathbb{Z}^n \setminus \{0\}} g(y) = \sum_{i = 1}^\infty \sum_{y \in m\mathbb{Z}^n, \|y\|_\infty = m\cdot i} f(y)
$$


$$
    = \sum_{i = 1}^\infty \sum_{y \in m\mathbb{Z}^n, \|y\|_\infty = m\cdot i} O\left(\frac{1}{(1 + \|y\|)^k}\right) = \sum_{i = 1}^\infty O\left(\frac{2^n m^{n-1}}{(1 + \|m i\|)^k}\right) = o(1),
$$

by choosing $$k$$ sufficiently large. The same goes for $$\widehat{g}(0) - \widehat{f}(0)$$.

The second question, however, is more complicated. Essentially, it is possible (and in fact to be expected a priori) that an optimal $$g_m^*$$ by solving the finite linear programming problem does not correspond to any admissible function $$f$$, so there might be a gap between the finite case and the infinite case.

My guess is that although there might not always be a function $$f^*$$ giving rise to $$g_m^*$$, one might always approach the bound given by $$g_m^*$$ with a function $$f$$, and moreover, perhaps one might always find a function $$f$$ giving rise to an auxiliary function $$g_m$$ arbitrarily close to $$g_m^*$$.

### References

{% bibliography --cited --file on_sphere_packing %}