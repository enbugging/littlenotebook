---
layout: post
title: "Relative position problem: motivation"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\proj}[1]{\text{proj }(#1)}
\newcommand{\vol}[1]{\text{vol}(#1)}
$$

This is mainly taken from Chapter 4 and 5 of my Bachelor thesis at École Polytechnique under supervision of Prof. Sarah J. Berkemer and Prof. Yann Ponty, titled ["Polynomial-time parametric optimisation"](https://arxiv.org/abs/2304.14962)

## 1. Introduction

We have seen the motivation and the mathematical construction of RNA polytope [(Nguyen, 2024)](intro_to_RNA_polytopes.html), which simply amounts to translating the dynamic programming recursion equation into the new polytope algebra, where addition is taking the convex hull, and multiplication is taking the Minkowski sum of two polytopes.

Unfortunately, supporting both operations of the polytope algebra poses a real challenge: whilst DNA sequence alignment and base-pair counting model each have 3 features, Turner energy model has close to 8000. Suppose we can effectively construct the feature vector from a given pair of sequence and secondary structure, the large dimensionality will induce a bottleneck. For general dimension $d$, there are currently two approaches.

## 2. Computation in full dimension

One may attempt to compute $\mathcal{P}(q)$ in full, yet the challenge lies on how one represent the polytopes. Let $A$ and $B$ be two polytopes, then $A \oplus B$ and $A \otimes B$ are efficiently computable only when $A$ and $B$ are given in $\mathcal{H}$- and in $\mathcal{V}$-representation, respectively {% cite Tiwary2008 -f rna_polytopes %}.

We may attempt to maintain both representations: this is known as the Double Description method, whose dual version is known as Beneath-and-Beyond algorithm, but this approach is also not computationally feasible {% cite Bremner1999 -f rna_polytopes %}. In practice, even in small dimensions, this approach still causes problems as the length of $q$ increases: when we restrict ourselves to the study of multi-loops in Turner model, which corresponds to 3 features, Poznanovic et al. {% cite Poznanovic2021 -f rna_polytopes %} showed that computing $\mathcal{P}(q)$ when $q$ is a tRNA of length 50 nts took approximately 2 hours, and when $q$ is a 5S rRNA, it took on average 23 hours. They also reported that a sequence of length 175 nts would increased the computation time to a week, and for $q$ of length 354 nts would take more than 2 months. Thus, it cannot practically cover the effective length for which Turner model finds its applications, which is up to 700 nts {% cite Mathews1999a -f rna_polytopes %}, and not scalable for the data we have.

It is natural that as the length of $q$ increases, so will the complexity of $\mathcal{P}(q)$. On the one hand theory, the result above by Zuker and Sankoff provides some intuition. On the other hand, much less can be said about the asymptotic complexity of $\mathcal{P}(q)$, as it does not solely depend on the number of points $c(q, s)$ but also the _distribution_ of such points {% cite Har-Peled2011 -f rna_polytopes %}, which has not been well-studied. In practice, Poznanovic et al. {% cite Poznanovic2021 -f rna_polytopes %} demonstrated that an increase of less than 50 nts multiplies that number by a factor of 3.5. 

To confirm this state of affairs, we implement the modified Nussinov's algorithm naively in Python, where the two operations are supported by \texttt{SageMath} software system. A random sequence of 10 nts took approximately 0.5 seconds, but that of 100 nts took 20 minutes, and increasing the length by 50 nts extended the runtime to 1 hour, agreeing with observations by Poznanovic et al. {% cite Poznanovic2021 -f rna_polytopes %} that the bottleneck lies in the algorithmic aspect, i.e. the inherent difficulty of the computational geometry problem involved, and not of sequence. 

We may also compute both operations in the same representation, but Tiwary {% cite Tiwary2008 -f rna_polytopes %} showed that this problem is output-sensitive _strongly_ $\textsf{NP}$-hard. In other words, he showed that unless $\textsf{P} = \textsf{NP}$, there exists no algorithm to this problem whose complexity is polynomial with respect to the final polytope's complexity (which we recall to possibly be exponential of $q$'s length) and the maximum absolute value of coefficient.

Another idea is that if we represent polytopes as sets of vertices, as the dynamic programming algorithm executes, $A$ and $B$ are not sets of random points, and thus the employment of a general-purpose convex hull algorithm to compute $A \oplus B$ may not be necessary. In particular, we only need to _merge_ two polytopes, for which there exist such a method in case of 2- and 3-dimensional {% cite ORourke1998 -f rna_polytopes %}, summarised in Table 1 below.

<table text-align="center">
    <thead>
        <tr>
            <th rowspand="2">Dimension</th>
            <th colspan="2">Complexity</th>
        </tr>
        <tr>
            <th></th>
            <th>Merge 2 convex hulls</th>
            <th>Build from 2 convex hulls</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>$2$</td>
            <td>$O(n)$</td>
            <td>$O(n \log h)$</td>
        </tr>
        <tr>
            <td>$3$</td>
            <td>$O(n)$</td>
            <td>$O(n \log h)$</td>
        </tr>
    </tbody>
</table>

We remark that $n$ and $h$ in Table 1 concern the input and output of _individual_ operations in polytope algebra, and not those of the derived dynamic programming algorithms. In particular, for some energy models, it may be the case that the final polytope has few vertices, yet its construction is still time-consuming. And unfortunately, it is difficult to implement already in 3-dimensional {% cite ORourke1998 -f rna_polytopes %} and unclear how to generalise for higher dimensions, thus renders this approach intractable.

## 3. Dimensionality reduction: an useful heuristics of uncertain precision

A common approach is to consider a few features whilst fixing the remaining parameters, as explored by Poznanovic et al. {% cite Poznanovic2021 -f rna_polytopes %}. Unfortunately, we lose degrees of freedom in the process, which may make a pair $(q, s)$ learnable for the full model, but not for the restricted version. Or in other words, reducing dimension allows fewer signature to be learnable.

In our example earlier, suppose we ignore the feature $c_{GC}$, then we have only 9 signatures left, namely $(0, 0)$, $(0, 1)$, $(0, 2)$, $(1, 0)$, $(1, 1)$, $(1, 2)$, $(2, 0)$, $(2, 1)$, and $(2, 2)$. Drawing on the plane, one obtains a square grid, as shown in Figure 1. On the other hand, the signature $(1, 1)$ lies in the interior of the square and thus is not learnable.

![](/littlenotebook/assets/img/rna_polytopes/rna_polytope.png#center)
<div align="center">Figure 1a. RNA polytope in 3D $\mathcal{P}(q)$.</div>
![](/littlenotebook/assets/img/rna_polytopes/rna_polytope_2d.png#center)
<div align="center">Figure 1b. Its projection onto $c_{AU} - c_{GU}$ plane</div>

To see how this happens, let $d > 0$, $E$ be an energy model with a parameter vector $p = (p_1, p_2, ..., p_d)^T \in \mathbb{R}^d$, which assigns to each pair of sequence $q$ and secondary structure $s$ the energy $E(p, q, s) = \langle c(q, s), p \rangle$. We denote $\text{proj} : \mathbb{R}^{d+1} \mapsto \mathbb{R}^d$ the map projecting vectors of $\mathbb{R}^{d+1}$ onto $\mathbb{R}^d$ by omitting the last coordinate. Let $E'$ be a model whose energy function is $E'(p, q, s) = \langle \proj{c(q, s)}, \proj{p} \rangle$, $q$ be a sequence, $s$ be some secondary structure, $x = c(q, s)$, $x' = \proj{x}$.

Assuming we fix $p_d \neq 0$ and suppose, without loss of generality, that $p_d > 0$. Since if there exists a vector $r$ such that $E(r, q, s) = \max_{s'} E(r, q, s')$, then for any scalar $\lambda \geq 0$, one has $E(\lambda r, q, s) = \max_{s'} E(\lambda r, q, s')$, in fact one need not constraint the choices of parameters for $E$ to those of the form $$\begin{pmatrix} p' \\ p_d \end{pmatrix}$$ for $p' \in \mathbb{R}^{d-1}$, but in fact, one may choose any vector $r \in \mathbb{R}^d$ such that $r_d > 0$ and $E(r, q, s) = \max_{s'} E(r, q, s')$, and with $\lambda = \frac{p_d}{r_d}$, one obtains a satisfying parameter vector $p = \lambda r$. Unfortunately, there exists no such other vector $r$, for if $r_d = 0$ will lead to $p_d = 0$ for any choice of $\lambda$, and if $r_d < 0$, then choosing $\lambda = \frac{p_d}{r_d}$ results in $E(p, q, s) = \min_{s'} E(p, q, s')$. Thus, if $(q, s)$ is learnable for $E$ but only with parameter vectors $r$ such that $r_d < 0$, then $(q, s)$ is not learnable for $E'$, no matter how one varies other parameters.

In rough terms, for any non-zero parameter fixed, we lose a halfsphere of $\mathcal{S}(\mathcal{P})(q)$ so if we fix $k$ parameters, there remains only $2^{-k}$ part of the parameter space that we can explore. In case of Turner model where $k$ is more than 7500, one can see how this affects the performance. Unfortunately, this case is of our greatest interest, since specifying a parameter to be zero effectively ignores the corresponding feature entirely and does not bring any new information.

As to _how much_ we can lose, assuming that $\mathcal{P}(q)$ is the realisation of a random variable, Amenta and Ziegler {% cite Amenta1996 -f rna_polytopes %} showed that we are guaranteed to retain at least a certain portion of learnable structures regardless of $q$'s length. However, little is known about the precise bound.

## 3.  Relative position problem

We have reviewed in the previous section the two current approaches to the computational aspect, namely computing the full polytopes and/or considering its projection to a lower dimension. We showed that each approach has its own limits, but to conclude, suppose we wish to study Turner energy model, even if we combine both approach and consider only 3 features, as shown by Poznanovic et al. {% cite Poznanovic2021 -f rna_polytopes %} as well as by our experiments, the parametric analysis is still practically intractable.

Nonetheless, for some problems, it is not necessary to construct the whole polytope, especially when we concern only some properties of $\mathcal{P}(q)$ relevant to our pair $(q, s)$. Given an sequence $q$, instead of finding what secondary structures (or, to be more precise, those of what signatures) an energy model $E$ can predict, we can ask the following question:

"Suppose we observe a structure $s$ for a sequence $q$ in experiments, does there exist a parameter set $p$ of which $E$ can predict $s$?"

Mathematically speaking, we ask if there exists a parameter vector $p$ such that $\langle c(q, s), p \rangle = h_{\mathcal{P}(q)}(p)$, or equivalently, if $c(q, s)$ lies on the boundary of $\mathcal{P}(q)$, i.e. if $(q, s)$ is learnable for $E$. We call this the learnability problem, and note that although there is a difference amongst predicting the structure $s$, the signature $c(q, s)$, and the energy $\langle c(q, s), p \rangle$. But, if we can find a parameter $p$ for which $s$ is an optimal solution, Wuchty et al. {% cite Wuchty1999 -f rna_polytopes %} demonstrated an algorithm to recover _all_ secondary structures with energy close to $\langle c(q, s), p \rangle$. This plays a great role, considering that the majority of experimentally observed secondary structures have energy close to the minimum, and that in real world, an RNA may exhibit multiple suboptimal structures.

To generalise, we have a fast access to $\mathcal{P}(q)$ via the dynamic programming algorithm, which takes a parameter set $p$ and returns $h_{\mathcal{P}(q)}(p)$, or with traceback, $\sigma_{\mathcal{P}(q)}(p)$, and we wish to determine if $c(q, s) \in \partial P$. Thus, in a more general setting, we ask the following question called Relative position problem.

**Problem.** (Relative position problem) Given a polytope $P \subseteq \mathbb{R}^{d+1}$ represented by its supporting function $h_P$, its extremal function $\sigma_P$, and a point $x \in \mathbb{R}^{d+1}$. Determine if $x \not\in P$, $x \in \mathring{P}$, or $x \in \partial P$.

Note that by the relationship between Optimisation oracle and Separation oracle, we can effectively determine if $x \in P$ in weakly polynomial time, and thus hereinafter we shall assume $x \in P$.

It turns out that this roblem is fundamental in computation geometry, by its relation with collision detection problem, stated as follow:

"Given two polytopes $A$ and $B$ in $\mathbb{R}^{d+1}$, determine if they are disjoint."

Gilbert, Johnson, and Keerthi {% cite Gilbert1988 -f rna_polytopes %} showed that this problem can be reduced to determining of $A \ominus B$ contains the origin, where $A \ominus B = \\{a - b \mid (a, b) \in A \times B\\}$ denotes the Minkowski difference. They also showed that a given $p \in \mathbb{R}^{d+1}$, one has $h_{A \ominus B} (p) = h_A (p) - h_B(-p)$, thus although the Minkowski difference, much like the Minkowski sum, can be complicated (in the sense of high complexity), its supporting function is easy to compute. Finally, they demonstrated an algorithm for the following problem

Given a polytope $P$ represented by its supporting function, determine if $P$ contains the origin.

The key idea is Carathéodory's theorem, which states that for a given polytope $P \subseteq \mathbb{R}^{d+1}$, any point $x \in P$ can be written as a linear combination of $d+2$ points of $P$ (which one can relax to $d+2$ vertices of $P$). Their method, so-called GJK algorithm, aimed to find such a simplex, whose vertices lie on the boundary of $P$, that contains $x$. If no such simplexes can be found, $x \not\in P$, otherwise it is easy to check if x lies in the interior or on the boundary of $P$.

For our interest when the dimension is high, their method, or rather the subroutine, called Johnson's distance subalgorithm, to compute the distance from a point $x$ to a $d+1$-simplex on which relies GJK algorithm, has two main weaknesses.
- Firstly, it is not numerical stable enough, which was the reason why in the original paper {% cite Gilbert1988 -f rna_polytopes %}, the authors introduced a backup procedure.
- Secondly, and more importantly, for polytopes in $\mathbb{R}^{d+1}$, it requires computing the distance from $x$ to affine subspaces generated by all $2^{d+2}-1$ combinations of vertices.

There have been attempts to mitigate one or both the aforementioned issues: in particular, Cameron {% cite Cameron1997 -f rna_polytopes %} modified the order in which the combinations are checked, thus improving the performance in practice, but without any theoretical bound. More recently, Montanari et al. {% cite Montanari2017 -f rna_polytopes %} decreased the number of combinations to be checked to $2^{d+1}$, but no methods which require only a polynomial (with respect to $d$) number of checks are known. For most problems where GJK algorithm is used, $d$ is often fixed and small, typical $d = 2$, and thus one can afford such a number of checks. Unfortunately, our problem has high dimension, thus GJK algorithm is not applicable.

Snethen {% cite Snethen2008 -f rna_polytopes %} proposed another algorithm, called Minkowski Portal Refinement, hereinafter abbreviated as MPR algorithm, to the same problem, which relies on the same principle, but with some modifications which allow reducing the number of checks from $2^{d+1}$ to $d+1$.

To summarise the idea, one first finds and fixes a point $o$ in the interior, then find a _portal_ defined as a $d$-simplex $P'$ through $d$ points on the boundary, which oftentimes are some vertices of $P$. Such a portal is said to be satisfying if the line segment $[o, x]_{\mathbb{R}^{d+1}}$ intersects $P'$. If that be the case, let $x'$ be the orthogonal projection of $x$ onto $P'$, the point $y = \sigma_P(x - x')$ forms with $P'$ a $d+1$-simplex. Since $[o, x]$ passes through $P'$, if $x$ does not lie in the interior of $P'$, $[o, x]$ must pass through another facet of this $d+1$-simplex that contains $y$. This new facet plays the role of the new portal, and the algorithm continues to termination.

In the original paper {% cite Snethen2008 -f rna_polytopes %}, Snethen did not specify a definite method to find an initial portal. He proposed fixing the point $o$ first, then finding a portal, but this does not guarantee to work on the first attempt, and may need multiple tries. In another post [(Nguyen, 2024)](relative_position_problem_simplex_method.html), I propose another approach, where we find a non-degenerate $d+1$-simplex $Q$ whose vertices lie on the boundary of $P$ first, then choose $o$ as a point lying in the interior of $Q$, e.g. the centroid of $Q$ as it is convex.

Strangely, MPR algorithm is little mentioned in the literature, to which Neumayr and Otter {% cite Neumayr2017 -f rna_polytopes %} addressed by describing an improved version with handling of termination conditions and special cases in 3-dimensional space. To the extend of our knowledge, neither full descriptions in general dimension nor proofs of termination for MPR are known. Another post [(Nguyen, 2024)](relative_position_problem_simplex_method.html) shall present more detailed description of MPR algorithm in the context of our scheme, and aim at addressing the aforementioned issues.

Recently, Hornus {% cite Hornus2017 -f rna_polytopes %} demonstrated another approach to the problem with a method called Decision Sphere Search, hereinafter abbreviated as DSS algorithm.

To summarise, suppose $x$ lies on the boundary of $P$, and let $S$ be its corresponding normal spherical polytope. In one iteration, one feeds a vector $p \in \mathbb{S}^d$ and obtains the great circle defined by $C = \mathbb{S}^d \cap \\{y \mid \langle y, p \rangle = \langle x, p \rangle \\}$ which divides $\mathbb{S}^d$ into two halfspheres, namely $C_{+} = \mathbb{S}^d \cap \\{y \mid \langle y, p \rangle > \langle x, p \rangle \\}$ and $C_{-} = \mathbb{S}^d \cap \\{y \mid \langle y, p \rangle < \langle x, p \rangle \\}$. If $h_P(p) = \langle x, p \rangle$, we conclude that $x \in \partial P$. Else, $S$ must lie in the halfsphere $C_{-}$, whose centre is given by $p' = \frac{x - \sigma_P(p)}{\\| x - \sigma_P(p) \\|}$. This new vector $p'$ is then fed into the next iteration.

The collection of halfspheres from the execution defines a spherical polytope $S'$ that bounds $S$, i.e $S \subseteq S'$. At some point, either we come across a vector $p \in S$ which proves $x \in \partial P$, or prove that $S$ is necessarily empty by showing that $S'$ is empty. This approach has its own limitations, amongst which is the lack of bound on complexity that is independent of the volume $\vol{S}$ of $S$, which we shall overcome in another post [(Nguyen, 2024)](relative_position_problem_ellipsoid_method.html) with ideas from ellipsoid method in linear programming.

### References 
{% bibliography --cited -f rna_polytopes %}