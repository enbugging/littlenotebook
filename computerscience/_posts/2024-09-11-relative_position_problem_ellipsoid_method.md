---
layout: post
title: "Relative position problem: ellipsoid method"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\vol}[1]{\text{vol}\left(#1\right)}
\newcommand{\conv}[1]{\text{conv}\left(#1\right)}
\newcommand{\sgn}[1]{\text{sgn}\left(#1\right)}
$$

(This is mainly taken from Chapter 8 of my Bachelor thesis at École Polytechnique under supervision of Prof. Sarah J. Berkemer and Prof. Yann Ponty, titled ["Polynomial-time parametric optimisation"](https://arxiv.org/abs/2304.14962))

## 1. Introduction

We have seen the motivation and the statement of Relative position problem [(Nguyen, 2024)](relative_position_problem_motivation.html), as well as a (family of) method akin to simplex method in linear programming [(Nguyen, 2024)](relative_position_problem_simplex_method.html). As I have shown, simplex method with an appropriate choice of pivot rule can demonstrate a great performance. Nonetheless, we see that in the worst case, its complexity can be linear with respect with the number of vertices of polytope, which, in the context of RNA polytopes, can be exponential of the sequence's length. On the one hand, in theory, this means we have not solved Relative position problem in polynomial time with respect to the dimension. On the other hand, in practice, a RNA sequence can be as long as 200,000 nts, and even if one limits themselves to the effective range of thermodynamics method, it can still be up to 700 nts. This may prove to be a bottleneck.
	
In this post, inspired by work of Hornus {% cite Hornus2017 -f rna_polytopes %}, I will develop an algorithm based on the ideas from ellipsoid method in linear programming, which motivates our naming. As we will see, the proof of convergence for both methods are closely related. In particular, it allows us to solve Relative position problem in weakly polynomial time.

**Theorem 1.** Relative position problem is solvable in weakly polynomial time.

Similar to other posts, in this post, we consider the case where the polytope $P$ is non-degenerate. Otherwise, $P$ has empty interior, so $P = \partial P$, and if $x \in P$, we return $x \in \partial P$.

## 2. Notations and preliminary on pherical geometry

We keep notation used in my other posts on the topic [(Nguyen, 2024)](intro_to_RNA_polytopes.html). Regarding the notions introduced below, readers are referred to Ratcliffe's _Foundations of Hyperbolic Manifolds_ {% cite Ratcliffe2007 -f rna_polytopes %} for more information.

Let $k \leq d+1$, an intersection of $\mathbb{S}^d$ and an $k$-plane is called a great $k$-sphere, which for terminology consistency, we shall call a $k$-plane whenever the context is clear. The case $k = 2$ gives great circles dividing $\mathbb{S}^d$ into halfspheres. And, just as any line segment is the shortest path connecting the two endpoints, any arc of a great circle is a geodesic on $\mathbb{S}^d$, which gives the notion of distance between two points $x, y \in \mathbb{S}^d$ as the length of the arc connecting them, which is necessarily unique (even if such an arc is not). Geometrically, this is equal to the angle between the two vector $x$ and $y$, giving the formula $d_{\mathbb{S}^d}(x, y) = \cos^{-1} \langle x, y \rangle$.

Given this metric, given two points $x, y \in \mathbb{S}^d$ and $d \geq d\_{\mathbb{S}^d}(x, y)$, we can define a spherical ellipsoid $\mathcal{E}$ on $\mathbb{S}^d$ as the set of points $z$ such that 
\\[
    d\_{\mathbb{S}^d} (x, z) + d\_{\mathbb{S}^d} (z, y) \leq d.
\\]
Equivalently, $\mathcal{E}$ can be viewed as the intersection of $\mathbb{S}^d$ and an elliptic hypercone.

For any two distinct non-antipodal points $x, y \in \mathbb{S}^d$ which we shall call a proper pair of points, there exists a unique geodesic segment defined by the arc connecting them, which we denote as $[x, y]_{\mathbb{S}^d}$. The subscript shall be neglected whenever the context is apparent.

## 3. Geometric ideas

For $u \in \mathbb{S}^d$, we denote $C(u) = \mathbb{S}^d \cap \\{v \mid \langle u, v \rangle = 0\\}$,  $C_{\leq 0}(u) = \mathbb{S}^d \cap \\{v \mid \langle u, v \rangle \leq 0\\}$, and $C_+(u) = \mathbb{S}^d \cap \\{v \mid \langle u, v \rangle > 0\\}$. Then, we have the two following lemmas.

**Lemma 2.** (rephrasing from Hornus's lemma {% cite Hornus2017 -f rna_polytopes %}) Let $x, y \in \partial P$ such that $x \neq y$, $u = \frac{y - x}{\\| y - x \\|}$ and $S$ be the normal spherical polytope of $x$. Then, one has $S \subseteq C_{\leq 0}(u)$.

_Proof._ For any $v \in C_+(u)$, one has $\frac{1}{\\| y - x \\|} (\langle y, v \rangle - \langle x, v \rangle) = \left\langle \frac{y - x}{\\| y - x \\|}, v \right\rangle > 0$, so $\langle y, v \rangle > \langle x, v \rangle$ and $v \not\in S$. Therefore, $S \subseteq \mathbb{S}^d \setminus C_+(u) = C_{\leq 0} (u)$. <span style="float:right;">$\square$</span>

**Lemma 3.** Let $x \in P$, $p \in \mathbb{S}^d$, and $y = \sigma_P(p)$. Assume that $\langle x, p \rangle < \langle y, p \rangle$, and let $u = \frac{y - x}{\\| y - x \\|}$ which is well-defined since $x \neq y$, then $p \in C_+(u)$.

_Proof._ By definition, $\langle y, p \rangle = h_P(p) > \langle x, p \rangle$, whence follows $\frac{1}{\\|y-x\\|}\langle p, y-x\rangle = \langle p, u \rangle > 0$, or $p \in C_{> 0}(u)$. <span style="float:right;">$\square$</span>

In simple terms, what these two lemmata say is that if $x \in \partial P$, it admits a normal spherical polytope $S$. Then, for any vector $p \in \mathbb{S}^d$, if $\langle x, p \rangle < h_P(p)$, there exists a hypercircle separating $p$ and $S$, given by $C(u)$ where $u = \frac{\sigma_P(p) - x}{\\|sigma_P(p) - x\\|}$. Analogous to ellipsoid method in linear programming, this plays the role of a separating affine hyperplane.

Now to construct the ellipsoid, suppose that after some $i$ iterations, we have bound $S$ by a spherical polytope $S_i$. Analogous to ellipsoid method in linear programming, we then bound $S_i$ by a spherical ellipsoid $\mathcal{E}_i$, and consider the center $u_i$ of $\mathcal{E}_i$. 

Choose $p = u_i$, and $u = \frac{\sigma_P(p) - x}{\\|\sigma_P(p)-x\\|}$. If $x$ lies on the boundary of $P$, by Lemma 2, $S$ will be contained in the region $\mathcal{E}\_i \cap C\_{\leq 0} (u)$, which can be bounded by another spherical ellipsoid $\mathcal{E}\_{i+1}$. By Lemma 3,  $u\_i \in \mathcal{E}\_i \cap C\_{> 0} (u)$, so we can hope that similarly to ellipsoid in linear programming, we have $\frac{\vol{\mathcal{E}\_{i+1}}}{\vol{\mathcal{E}\_i}} < C_1$ for some constant $0 < C_1 < 1$. As we shall see, we can associate to each spherical ellipsoid $\mathcal{E} \subsetneq \mathbb{S}^d$ an ellipsoid $E \subsetneq \mathbb{R}^d$ such that $\vol{\mathcal{E}}$ tends to $0$ if and only if $\vol{E}$ tends to $0$, and we show that indeed if we associate $E_i$ to $\mathcal{E}\_i$ and $E\_{i+1}$ to $\mathcal{E}\_{i+1}$, then $\frac{\vol{E\_{i+1}}}{\vol{E\_i}} < C_2$ for some $0 < C_2 < 1$.

### 3.1. Spherical ellipsoid representation

In $\mathbb{R}^{d+1}$, an ellipsoid $E$ can be represented by a positive definite matrix $Q$ and a center $c$, as $E = \\{x\mid(x-c)^T Q^{-1} (x-c)\\}$, which allows closed-form formulae for representing $E_i$, computing its volume, and practical implementation. Unfortunately, spherical ellipsoids lack such a representation, which not only makes it difficult to implement the algorithm, but also to compute $\vol{\mathcal{E}_i}$ even when $d = 2$. To overcome this obstacle, we propose the following solution.

Let $\mathcal{C}(\mathcal{E}) = \\{\alpha x\mid x \in \mathcal{E}, \alpha \geq 0\\}$, then $\mathcal{C}(\mathcal{E})$ is a convex cone such that $\mathcal{E} = \mathcal{C}(\mathcal{E}) \cap \mathbb{S}^d$. Let $u$ be the center of $\mathcal{E}$, and consider the hyperplane $H$ tangent to $\mathbb{S}^d$ at $u$. Identify $H$ with $\mathbb{R}^d$ and $u$ with the origin, then we have that $\mathcal{C}(\mathcal{E}) \cap H$ is an ellipsoid $E$ in $\mathbb{R}^d$, centered at the origin. Moreover, by endowing $H$ with a suitable basis, we can even have $E$ admitting vectors in the canonical basis as its axes.

![](/emilesnotebook/assets/img/rna_polytopes/ellipsoid_representation.png#center)
<div align="center">Figure 1. Spherical ellipsoid and its representation. Solid line and dashed line represent $E_H$ and $\mathcal{E}$, respectively. The plane represents $H$.</div>

Note that there is an abuse of notation here. There is the non-degenerated ellipsoid in $\mathbb{R}^d$, and the degenerated ellipsoid in $\mathbb{R}^{d+1}$ which is contained in $H$. The letter $E$ is used to denote both objects, which causes no confusion due to our identification.

One can show that such a pair $(u, E)$ is uniquely defined, which gives a representation of $\mathcal{E}$. In what follows, whenever appropriate, to mean a spherical ellipsoid $\mathcal{E}$, we will use such a pair $(u, E)$ of a vector $u \in \mathbb{S}^d$ and an ellipsoid $E \subseteq \mathbb{R}^d$ given by a positive-definite matrix $Q$ as $E = \\{x \mid x^T Q^{-1} x \leq 1\\}$.

Finally, we need to show that $\vol{\mathcal{E}}$ tends to $0$ as $\vol{E}$ tends to $0$.

**Lemma 4.** With notations as above, as $\vol{E}$ tends to $0$, so does $\vol{\mathcal{E}}$.

_Proof._ The hyperplane $H$ divides $\mathbb{R}^{d+1}$ into two halves, one of which contains the origin. Denote this halfspace by $\mathcal{H}$, and denote $C = u \cap \mathcal{H}$, which is a cone with $E$ as its base. Denote $\mathcal{B} = B_d(0, 1)$ the unit ball in $\mathbb{R}^{d+1}$. Since $\mathcal{B} \subseteq \mathcal{H}$, we have that $\vol{u \cap \mathcal{B}} \leq \vol{C}$. Simple integration gives, 
\\[
    \frac{\vol{\mathcal{E}}}{d+1} =  \vol{u \cap \mathcal{B}} \leq \vol{C} = \frac{\vol{E}}{d+1},
\\]
as desired. <span style="float:right;">$\square$</span>

## 4. Method description

Similar to ellipsoid method in linear programming, our algorithm has two phases: finding an initial bounding spherical ellipsoid, and refining it.

### 4.1. Phase 1: Finding an initial bounding spherical ellipsoid

We construct a non-degenerate simplex $Q = \conv{\\{u_0, u_1, \cdots, u_{d+1}\\}}$ in $d+2$ iterations.

Then, we consider the points $\mathcal{u}\_i \in \mathbb{S}^d$ corresponding to the face formed by $\\{u\_j\mid j \neq i\\}$, and since any two facets are adjacent, in particular they cannot be parallel to each other, thus no two points $\mathcal{u}\_i$ and $\mathcal{u}\_j$ are antipodal. Thus, with these $d+2$ points, we have $d+2$ great hypercircles, denoted by $C(\mathcal{u\_0}), C(\mathcal{u\_1}), \cdots, C(\mathcal{u\_{d+1}})$, which divides $\mathbb{S}^d$ into spherical polytopes.

Now if the algorithm does not terminate by now, then it must be the case that for each $i = 0, 1, \cdots, d+1$, the normal spherical polytope $S \in \mathcal{S}(P)$ of $x$, if non-empty, must lie in one of the halfspheres defined by $C(\mathcal{u}\_i)$. Therefore, amongst the spherical polytopes of dimension $d$ generated by $C(\mathcal{u\_0}), C(\mathcal{u\_1}), \cdots, C(\mathcal{u\_{d+1}})$, exactly one contains $S$.

On the other hand, any normal spherical polytope $S' \in \mathcal{S}(Q)$ can be contained in a ball (in $\mathbb{S}^d$) whose radius is strictly less than $\frac{\pi}{2}$. Indeed, suppose the contrary, then there exists some normal spherical polytope $S' \in \mathcal{S}(Q)$ and some $p \in \mathbb{S}^d$ such that both $p$ and $-p$ are in $S'$. The polytope $S'$ corresponds to some point $x' \in \partial Q$, and by construction, one has that
\\[
    Q \subseteq \\{y\mid\langle y, p \rangle \leq \langle x', p \rangle\\} \cap \\{y\mid\langle y, p \rangle \geq \langle x', p \rangle \\} = \\{y\mid\langle y, p \rangle = \langle x', p \rangle\\},
\\]
contradicting the non-degeneracy of $Q$.

Therefore, suppose $S$ is contained in some $S' \in \mathcal{S}(Q)$, which in turns lies in the ball $B_{\mathbb{S}^d}(u, \theta)$ where $1 < \theta < \frac{\pi}{2}$, we can identify with $(u_1, E_1)$ where $E_1 = B(1, r) \subseteq \mathcal{R}^d$ is a ball of radius $r = \tan^{-1} \theta$ necessarily finite.

As for how to construct $u_1$ and $E_1$ algorithmically, to each great hypercircle $C$ that contains a facet of $S'$, we denote by $H_C$ the hyperplane containing $C$. The collection of such hyperplanes $H_C$ forms a cone $u$. Intersecting $u$ with a hyperplane $H$ tangent to $\mathbb{S}^d$ such that $u \cap H$ is finite, then it is easy to find a suitable $r$ for $E_1$.

### 4.2. Phase 2: Refinding spherical ellipsoid

At the $i^{\text{th}}$ iteration of this phase, we have our spherical ellipsoid $\mathcal{E}_i$ represented by the pair $(u_i, E_i)$. We test the vector $u_i$ and see if $\langle x, u_i \rangle = h_P(u_i)$. If this is true, then $x \in \partial P$; else, we introduce the separating great hypercircle given by $C(w)$ where $w = \frac{\sigma_P(u_i) - x}{\\|\sigma_P(u_i) - x\\|}$. Also consider the hyperplane $K$ such that $C(w) = K \cap \mathbb{S}^d$. Let $\ell = K \cap H_i$, we have two cases.

- Either $\ell$ does not intersect $E\_i$, at which point we know that $C(w)$ does not intersect $\mathcal{E}\_i$. Given that the center $u\_i$ lies in $C\_{> 0}$, it is thus necessarily the case that $S \subset \mathcal{E}\_i \subset C\_{> 0}(u\_i)$, but by Lemma 2, $S \subset C\_{\leq 0}(u\_i)$, so $S$ is empty, and $x \in \mathring{P}$.
- Or $\ell$ intersects $E_i$. We then construct an ellipsoid of minimal volume $E' \subset H$ containing the intersection of $E$ and the halfplane of $H$ defined by $\ell$ and not containing $u_i$.

Note that $E'$ does _not_ give the representation of $\mathcal{E}\_{i+1}$. Whilst one does have that $u_{i+1} = \frac{c'}{\\|c'\\|}$ where $c'$ is the center of $E'$, $H$, which contains $E'$, is tangent to $\mathbb{S}^d$ at $u\_i$ and not $u\_{i+1}$. To correct this, we then consider the cone $C' = \\{\alpha v\mid\alpha \geq 0, v \in E'\\}$, $H_{i+1}$ the affine hyperplane tangent to $\mathbb{S}^d$ at $u\_{i+1}$, and compute $E\_{i+1} = H\_{i+1} \cap C'$. Note that let $\mathcal{E}\_{i+1} = C' \cap \mathbb{S}^d$, then it must contain the intersection of $\mathcal{E}\_i$ and the halfsphere defined by $C(w)$ not containing $u\_i$.

Finally, if $\vol{E_{i+1}} \leq \varepsilon$ for some fixed $\varepsilon > 0$, then we terminate and report that $x \in \mathring{P}$, as $S$ is contained in a spherical polytope too small that it must be empty.

## 5. Complexity

With Lemma 4, it suffices to exhibit some constant $0 < C < 1$ depending only on $d$ such that $\frac{\vol{E_{i+1}}}{\vol{E_i}} \leq C$. This is the content of the following theorem.

**Theorem 5.** With notations defined as above and assume that either $d \geq 3$ or $d = 1$, one has that $\frac{\vol{E_{i+1}}}{\vol{E_i}} \leq e^{-\frac{1}{2(d+1)}}$.

_Proof._ A similar result from the analysis of ellipsoid method in linear programming states that $\frac{\vol{E'}}{\vol{E_i}} \leq e^{-\frac{1}{2(d+1)}}$ , so what is left to show is that $\vol{E_{i+1}} \leq \vol{E_i}$, or equivalently, the correction does not increase the volume of the ellipsoid. First, we decompose the correction into two steps:
1. We "rotate" the hyperplane $H$ into another hyperplane $H'$ passing through $c'$ and admitting $c'$ as its normal vector, i.e. $H'$ is parallel with $H_{i+1}$. Let $E'' = C' \cap H'$ be the ellipsoid $E'$ but "rotated" to hyperplane $H'$.
2. We move the hyperplane $H'$ toward the origin until it touches $\mathbb{S}^d$, i.e. when it coincides with $H_{i+1}$.

Let $e_1, e_2, \cdots e_d$ be the canonical basis of $\mathbb{R}^d$. By possibly a rotation, assume without loss of generality that $u_i = (0, \cdots, 0, 1)^T$, so 
$$H = \left\{\begin{pmatrix} z \\\ 0 \end{pmatrix} \Big| z \in \mathbb{R}^d\right\}$$
. Write $$c' = \begin{pmatrix} f \\ 1 \end{pmatrix}$$ where $f = (f_1, f_2, \cdots, f_d)^T \in \mathbb{R}^d$. Moreover, let $P$, $P'$, and $P_{i+1}$ be the positive-definite matrices corresponding to $E'$ in $H$, $E''$ in $H'$, and $E_{i+1}$ in $H_{i+1}$, respectively. With possibly a rotation, we can assume without loss of generality that $P$ admits the canonical basis of $\mathbb{R}^d$ as its eigenvectors. Let $\lambda_1, \lambda_2, \cdots, \lambda_d$ be the eigenvalues of $P$, and for $j = 1, 2, \cdots, d$, let $$y_j = \begin{pmatrix} \sgn{f_j} \lambda_j e_j + f \\\ 1 \end{pmatrix}$$, where $$\sgn{x} = \begin{cases} 1 & \text{ if } x \geq 0 \\\ -1 & \text{ if } x < 0 \end{cases}$$.
    
We describe the projection from $H$ to $H'$ where the point $$y = \begin{pmatrix} z \\ 1 \end{pmatrix} \in H$$ is mapped to the point $\mu(y) y \in H'$, where $z \in \mathbb{R}^d$. The condition $\mu(y) y \in H'$ is equivalently to 
\\[
    f^T f + 1 = \\| c' \\|^2 = \left\langle \mu(y) y, c' \right\rangle = \mu(y) (z^T f + 1) \Rightarrow \mu(y) = \frac{f^T f + 1}{z^T f + 1}.
\\]
In particular, $\mu(y_j) = \mu_j = \frac{f^T f + 1}{\lambda_j \sgn{f_j} f_j + f^T f + 1} = \frac{f^T f + 1}{\lambda_j \|f_j\| + f^T f + 1}$.

One key insight, whose proof shall be omitted for the sake of brevity, is under this projection, axes of $E'$ are mapped to those of $E''$. Whilst we cannot at the moment explicitly endow $H'$ with a basis, one can still compute the eigenvalues of $P'$ by $\\| \mu(y_j) y_j - c' \\|$. And finally by multiplying with a scaling factor $\frac{1}{\\|c'\\\|}$ , we determine the eigenvalues $\delta_j$ of $P_{i+1}$. Let $F = f^T f + 1$, calculation shows

$$
\begin{split}
    \| \mu(y_j) y_j - c' \|^2
    & = \left\| \mu_j \begin{pmatrix}
        \sgn{f_j} \lambda_j e_j + f\\
        1
    \end{pmatrix} - \begin{pmatrix}
        f \\
        1
    \end{pmatrix}\right\|^2 \\
        & = 
    \left\| \begin{pmatrix}
        \sgn{f_j} \mu_j \lambda_j e_j + (\mu_j - 1) f\\
        \mu_j - 1
    \end{pmatrix}\right\|^2 \\
    & = (\mu_j - 1)^2 + (\sgn{f_j} \mu_j \lambda_j - (\mu_j - 1) f_j)^2 + (\mu_j - 1)^2 \sum_{j \neq i} f^2_j \\
    & = (\mu_j - 1)^2 + \mu^2_j \lambda^2_j - 2\sgn{f_j}\mu_j (\mu_j - 1) \lambda_j f_j \\
        & + (\mu_j - 1)^2 f^2_j + (\mu_j - 1)^2 \sum_{k \neq j} f^2_k \\
    & = \mu^2_j \lambda^2_j - 2\mu_j (\mu_j - 1) \lambda_j |f_j| + (\mu_j - 1)^2 F \\
    & = \left(\frac{F}{\lambda_j |f_j| + F}\right)^2 \lambda^2_j - 2 \frac{F}{\lambda_j |f_j| + F} \left(\frac{F}{\lambda_j |f_j| + F} - 1\right) \lambda_j |f_j| \\
    & + \left(\frac{F}{\lambda_j |f_j| + F} - 1\right)^2 F \\
    & = \left(\frac{F}{\lambda_j |f_j| + F}\right)^2 \lambda^2_j + 2 \frac{F\lambda^2_j f^2_j}{(\lambda_j |f_j| + F)^2} + \frac{\lambda^2_j f^2_j}{(\lambda_j |f_j| + F)^2} F \\
    & = \left(\frac{F}{\lambda_j |f_j| + F}\right)^2 \lambda^2_j + 3 \frac{F\lambda^2_j f^2_j}{(\lambda_j |f_j| + F)^2}\\
\Rightarrow \delta^2_j
    & = \frac{1}{F} \| \mu(y_j) y_j - c' \|^2 = \frac{3f^2_j\lambda^2_j}{\lambda_j |f_j| + F} + \frac{F\lambda^2_j}{(\lambda_j |f_j| + F)^2} = \frac{3f^2_j + F}{(\lambda_j |f_j| + F)^2} \lambda^2_j\\
\end{split}
$$
By Cauchy-Schwartz inequality, one obtains

$$
\begin{split}
    \prod_{j = 1}^{d} \frac{3f^2_j + F}{(\lambda_j |f_j| + F)^2}
    & \leq \prod_{j = 1}^{d} \frac{3f^2_j + F}{F^2} = \prod_{j = 1}^{d} \left(\frac{1}{F} + \frac{3f^2_j}{F^2}\right) \\
    & \leq \left[\frac{1}{d}\sum_{j=1}^{d} \left(\frac{1}{F} + \frac{3f^2_j}{F^2} \right)\right]^d = \left[\frac{1}{F} + \frac{3}{dF}- \frac{3}{dF^2}\right]^d \leq 1
\end{split}
$$
where the last inequality holds for $d \geq 3$, as
\\[
    1 - \frac{1}{F} - \frac{3}{dF} + \frac{3}{dF^2} = \left(1 - \frac{1}{F}\right)\left(1 - \frac{3}{dF}\right) \geq 0 \Rightarrow \frac{1}{F} + \frac{3}{dF} - \frac{3}{dF^2} \leq 1.
\\]
As for the case $d = 1$, it is clear that one has $\lambda_1 = 2|f_1|$, so in fact we have a stronger inequality,
\\[
    \frac{3f^2\_1 + F}{(\lambda\_1 \|f\_1\| + F)^2} = \frac{4f^2\_1 + 1}{(3f^2\_1 + 1)^2} \leq 1.
\\]
Finally, we conclude by remarking that
\\[
    \vol{E\_{i+1}} = \vol{B\_d(0, 1)}\left(\prod\_{j = 1}^d \delta_j\right)^{\frac{1}{2}} \leq \vol{B\_d(0, 1)}\left(\prod\_{j = 1}^d \lambda\_j\right)^{\frac{1}{2}} = \vol{E'}, 
\\]
where $B_d(0, 1)$ denotes the unit ball in $\mathbb{R}^d$, as desired. <span style="float:right;">$\square$</span>

We have attempted to prove for $d = 2$, but the calculations proved to be cumbersome. However, we remark that although a priori we do not know if the convergence rate holds for $d = 2$, this poses no difficulty, as one can easily lift a polytope from $\mathbb{R}^3$ to $\mathbb{R}^4$ by considering  the following polytope $Q$.

Let $z = (z\_1, z\_2, z\_3)^T \in P$, and consider two points $z\_+ = (z\_1, z\_2, z\_3, 1)^T$ and $z\_- = (z\_1, z\_2, z\_3, 1)^T$. Let $P' = \\{(y\_1, y\_2, y\_3, 0)^T \mid (y\_1, y\_2, y\_3)^T \in P\\}$. Finally, let $Q = \conv{\\{z\_+, z\_-, P'\\}}$, and $x' = (x^T, 0)^T = (x\_1, x\_2, x\_3, 0)^T$, then we run the algorithm with $Q$ and $x'$. One can see that $x \in \partial P$ if and only if $x' \in \partial Q$, so the algorithm runs correctly.

Thus for a given $\varepsilon$, for $d \neq 2$, the algorithm will terminate after at most $d+2 + 2(d+1)\ln \frac{\vol{E\_1}}{\epsilon}$ iterations, and in case $d = 2$, it will terminate after at most $d+3 + 2(d+2)\ln \frac{\vol{E\_1}}{\epsilon}$ iterations. Each iteration calls the supporting function $h_P$ and the extremal function $\sigma_P$ exactly once, and all the other linear algebra operations are performed in polynomial time with respect to $d$.

We can also write formally
$$
Q = \conv{\left\{\begin{pmatrix}
		z \\ 1
	\end{pmatrix}, \begin{pmatrix}
		z \\ -1
	\end{pmatrix}\right\} \cup \left\{\begin{pmatrix}
		y \\
		0 
	\end{pmatrix} \Big| y \in P\right\}}
$$
Finally, let $$x' = \begin{pmatrix} x \\ 0 \end{pmatrix}$$, and we run the algorithm with $Q$ and $x'$. One can see that $x \in \partial P$ if and only if $x' \in \partial Q$, therefore the algorithm will return correctly.

Thus, for a given $\varepsilon > 0$, the algorithm will terminate after at most $d+2 + 2(d+1)\ln \frac{\vol{E_1}}{\varepsilon}$ iterations, each iteration calls the supporting function $h_P$ and the extremal function $\sigma_P$ exactly once, and all the other linear algebra operations are performed in polynomial time with respect to $d$.

Now recall that for a rational number $x$, its size, denoted by $\langle x \rangle$, is the number of bits needed to represents $x$. For a vector $x = (x\_1, x\_1, .., x\_{d+1})^T \in \mathbb{Q}^d$, we denote $\langle x \rangle = \langle \max\_{i} \|x\_i\| \rangle = \max\_{i} \langle \|x\_i\| \rangle$. And we call the size of $P$, denoted $\langle P \rangle$, to be the maximum size of its vertices. 

The final catch is that weakly polynomial time is defined in $d$ and sizes $\langle P \rangle = D$ of coordinates that represent $P$'s vertices. What is left to prove, is that $\langle \vol{E_0} \rangle$ is also bounded by polynomial of $D$. Here, for the sake of brevity, we only present a sketch of the proof.
1. Since the coordinates of vertices have size bounded by $D$, the "area" of facets and the volume of $Q$ have size to be polynomial of $D$. This leads to the dihedral angles of $Q$ also having size to be polynomial of $D$.
2. Moving to normal polytopes, since the smallest distinguishable angle has size of polynomial of $D$, so are the length of $[x, y]\_{\mathbb{S}^d}$ for any $x, y \in S$ where $S \in \mathcal{S}(Q)$. It follows that the smallest ball $B\_{\mathbb{S}^d}(u, \theta)$ that contains $S$ must also has $\theta$ of size $D$.
3. Finally, recall that as $x$ tends to $\frac{\pi}{2}$, one has $\tan x = \frac{1}{\tan\left(\frac{\pi}{2} - x\right)} \approx \frac{1}{\frac{\pi}{2} - x}$, so in general, $\tan x$ has size of polynomial with respect to that of $x$. Projecting $S$ on $H$, it follows that $E_1 = B(0, r)$ has its radius $r$ whose size is of polynomial with respect to that of $\tan \theta$.

Therefore, the algorithm complexity is polynomial in $d$ and $\langle P \rangle$, which, together with the method using the relation between Optimisation oracle and Separation oracle as presented at the beginning, solves Relative position problem.

## 6. Concluding remarks

In this concluding remark, we wish to go beyond the scope of the two presented methods, and draw attention to the correspondence between methods in linear programming and ours. That is, the simplex method performs great in practice, but has poor worst-case performance, whereas the ellipsoid method has the theoretical guarantee, but involves unstable numerical operations and has slow convergence. Moreover, the simplex method both in linear programming and in our setting has its performance sensitive to the choice of pivot rule, and the ellipsoid method has similar convergence rate as expected from the fact that the proofs are closely related.

This begs the question if there exists an analogue of the interior-point method for Relative position problem, with potentially good performance in practice and a theoretical convergence guarantee. It is not obvious what would play the role of an interior point in this case, and we leave this direction as open for future research.

In the grand scheme of linear programming, each of the three methods opens up new research directions of their own: study of pivot rules for simplex method, study of cuts for ellipsoid method, and study of barriers for interior-point methods. With the correspondence given above, it is easy to see what directions for the study of Relative position problem. We have demonstrated that MPR algorithm can be outperformed by other pivot rules, but have yet to show any such rules that is guaranteed to terminate. We also have shown a possible representation for spherical ellipsoids, chosen for convenience and the fact that the spherical ellipsoids only represent cones, and it is these cones that are of our interests. Nonetheless, other representations are possible, such as intersections of $\mathbb{S}^d$ with elliptic cylinders, and this may improve numerical stability. We leave this hypothesis to be proven or disproven in the future.

And finally, as an analogue to linear programming, we ask if Relative position problem can be solved in _strongly_ polynomial time. Given the dependence upon the relationship between Optimisation oracle and Separation oracle, it is more difficult than that in linear programming, but even assuming $x \in P$, this question is interesting in its own right.

### References

{% bibliography --cited -f rna_polytopes %}