---
layout: post
title: "On critical sets in stable matching"
author: "Nguyen Doan Dai"
categories: math
tags: math
---

$$
\newcommand{\NP}{\textsf{NP}}
\newcommand{\def}{\text{def}}
\newcommand{\d}{\text{d}}
\newcommand{\FPT}{\textsf{FPT}}
\newcommand{\APX}{\textsf{APX}}
\newcommand{\sgn}{\text{sgn}}
$$

This is a study during my second year at ENS Ulm, in the context of a research project supervised by Prof. Ana Busic.

## 1. Introduction

Let us take a real-life example: Uber application. There are passengers, and there are drivers. Passengers demand taxis; drivers provide service. There are a few categories of drivers: those with big cars, small cars, electric cars, handicap-friendly cars, etc. There are also a few categories of passengers: groups, individuals, those prefer electric cars for environmental reasons, those who need assisstance, those who need big cars because they have a lot of luggages, etc. At any given moment, a driver or a passenger of any type can arrive independent from each other. A driver can be matched with multiple passengers, and vice versa, but ultimately, we can only pair one driver with one passenger (group). And the whole is to devise a matching policy so that the expected waiting time for passengers, and for drivers once they are free, is finite.

One can encounter similar models elsewhere: organ donors, dating application, submitted articles versus reviewers, etc. In many cases, there is not even a clear two-side separation as in passengers versus drivers. This motivated the notion of stable bipartite matching model introduced by Caldentey, Kaplan, and Weiss {% cite Caldentey2009 -f on_critical_sets_of_stable_matching.bib %} under the assumption that passengers and drivers arrive independently. This assumption was later removed in the work of Busic, Gupta, and Mairesse {% cite Busic2013 -f on_critical_sets_of_stable_matching.bib %} and later generalised to general matching graphs by Mairesse and Moyal {% cite Mairesse2016 -f on_critical_sets_of_stable_matching.bib %}.

In their paper, Busic et al. derived a necessary condition for a bipartite matching model to be stabilisable, i.e. there exists a matching policy that makes the model stable. The condition is analogous to the condition given in Hall's marriage theorem to guarantee the existence of a perfect bipartite matching. Hall's marriage theorem is generalised by Tutte's theorem for general graphs; and similarly, Mairesse and Moyal {% cite Mairesse2016 -f on_critical_sets_of_stable_matching.bib %} generalised the aforementioned necessary condition for general matching graphs.

Now let us use the analogy of graph matching problem once again. For bipartite graphs, when Hall's condition does not hold, Hall's theorem says that there is no perfect matching, but it does not tell us the size of the maximum matching. In that regard, Ore {% cite Ore1955 -f on_critical_sets_of_stable_matching.bib %} introduced the notion of deficiency. Recall that for a bipartite graph $G = (A \cup B, E)$, the deficiency of $G$ with respect to $A$ is the maximum, over all subsets $U \subseteq A$, of $\|U\| - \|N_G(U)\|$. Then if $d$ is the deficiency of $G$, $G$ admits a matching of size $\|X\| - d$. In a way, the notion of deficiency tells us how close a bipartite graph is to having a perfect matching.

In this post, we introduce a similar notion for graph matching models. In rough words, the necessary condition is an inequality, and we want to calculate the minimum gap between the two sides. Moreover, we also want to identify a subset of vertices that realises this minimum. This corresponds, in a sense, to the most critical part of a matching model, i.e. the part that a smallest amount of pertubation can make the whole model unstable.

## 2. Preliminary 

### 2.1. Matching models 

Let $G = (V, E)$ be an undirected graph, called the matching graph. We equip it with a probability distribution $\mu$ on $V$, which denotes the probability of arrival of each type of item at an instance, i.e. at a given moment $t$, the probability that an item of type $u$ arrives is $\mu(u)$. We assume that items arrive independently.

Suppose at time $t$, an item of type $v$ arrives. If there exists an item of type $u$ in the buffer such that $(u, v) \in E$, then we may match them, after which $u$ is removed from the buffer. Otherwise, $v$ is added into the buffer. When there are multiple choices for $u$, we need to specify a matching policy to decide which item we match with $v$.

The matching model for bipartite graph is defined similarly, except that now we have a bipartite graph $G = (A \cup B, E)$ to which we equip with a distribution $\mu$ on $A \times B$, denoting that at a given time $t$, an element $a \in A$ and $b \in B$ arrives with probability $\mu((a, b))$, which are then added into the buffer. Then, if there are two vertex $u \in A$ and $v \in B$ in the buffer such that $(u, v) \in E$, we can match them. And likewise, when there are multiple choices, we specify a matching policy to decide which item we choose. Denote $\mu_A$ and $\mu_B$ the marginal of $\mu$ on $A$ and $B$.

### 2.2. Deficiency and critical sets

Busic et al. {% cite Busic2013 -f on_critical_sets_of_stable_matching.bib %} showed that a necessary condition for a bipartite matching model to be stablisable is that

$$
\begin{cases}
    \mu_A(U) < \mu_B(N(U)) & \forall U \subsetneq A \\
    \mu_B(V) < \mu_A(N(V)) & \forall V \subsetneq B \\
\end{cases}
$$

Similarly, Mairesse and Moyal {% cite Mairesse2016 -f on_critical_sets_of_stable_matching.bib %} showed that a necessary condition for a matcing model to be stabilisable is that for all independent sets $U \subset V$, we have
\\[
    \mu(U) < \mu(N(U)),
\\]
where $\mu(U) = \int_U d\mu = \sum_{u \in U} \mu(u)$, and likewise for $\mu(N(U))$, $\mu_A(U)$, $\mu_A(N(U))$, $\mu_B(V)$ and $\mu_A(N(V))$.

Analogous to the notion of deficiency in graph theory, for an independent set $U \subseteq V$, we define the deficiency of $U$ to be 
\\[
    \def_{G, \mu} (U) = \mu(N(U)) - \mu(U),
\\]
and the deficiency of matching model $(G, \mu)$ to be 
\\[
    \def(G, \mu) = \min_{U \subseteq V \text{ independent}} \def_{G, \mu} (U).
\\]

An independent set $U \subseteq V$ such that $\deg_{G, \mu} (U) = \def(G, \mu)$ is called a critical set.

The problem then is as follow: given $G$ and $\mu$, find a critical set.

## 3. Finding a critical set is tractable

A priori, there is no reason for this problem to be tracable: effectively, we are optimising over the set of independent sets. A similar problem, finding the maximum independent set, is well-known to be $\NP$-hard. That being said, even in the case of finding maximum independent set, we can still find approximation algorithms. One of the most well-known approach is to pose it an integer linear programming problem, then we relax the constraints to obtain a linear programming problem, for which we have efficient algorithms. In case of maximum independent set and vertex cover, we can obtain a half-integral optimal solution, which then provides a 2-approximation solution to the original problem.

Let us follow the same approach.  Assign to each vertex $u$ a variable $x_u \in \\{-1, 0, 1\\}$, where $x_u = -1$ means $u \in U$ and $x_u = 1$ means $u \in N(u)$. Then we can write 
\\[
    f(U) = \sum_{u \in V} c_u x_u.
\\]
Such an assignment $(x\_u)\_{u \in V}$ satisfies $x_u + x_v \geq 0$ for all $(u, v) \in E$; conversely, if an assignment $x = (x_u)\_{u \in V} \in \\{-1, 0, 1\\}^V$ satisfies $x_u + x_v \geq 0$ for all $(u, v) \in E$, then no two adjacent vertices $u$ and $v$ can have $x_u = x_v = -1$, and so the set $U = \\{u \mid x_u = -1\\}$ is independent, and the two conditions are equivalent.

Then the problem now is an integer programming problem $(1)$.

$$
\begin{equation}
    \begin{array}{ll@{}ll}
        \text{minimise}  & \displaystyle\sum_{u \in V} c_u x_u &\\
        \text{subject to}& x_u + x_v \geq 0, &(u, v) \in E &\\
        &x_u \in \{-1, 0,1\}, &u \in V
    \end{array}
\end{equation}
$$

Relaxing $(1)$ gives a linear programming problem $(2)$.

$$
\begin{equation}
    \begin{array}{ll@{}ll}
        \text{minimise}  & \displaystyle\sum_{u \in V} c_u x_u &\\
        \text{subject to}& x_u + x_v \geq 0, &(u, v) \in E &\\
        &-1 \leq x_u \leq 1, &u \in V
    \end{array}
\end{equation}
$$

Note that vertex cover problem can be posed as 

$$
\begin{equation}
    \begin{array}{ll@{}ll}
        \text{minimise}  & \displaystyle\sum_{u \in V} c_u x_u &\\
        \text{subject to}& x_u + x_v \geq 1, &(u, v) \in E &\\
        &x_u \in \{0, 1\}, &u \in V
    \end{array}
\end{equation}
$$

and its LP relaxation always has a half-integral optimal solution, thus we have reason to believe the same for our problem.

What I did not expect, however, is that our LP relxation always has _integral_ optimal solution.

**Lemma 1**. Problem $(2)$ is totally dual integral.

_Proof._ This is equivalent to saying that for all weight $c$, there exists an integral optimal solution. Indeed, fix $c$ and let $x = (x_u)\_{u \in V}$ be an optimal solution. Denote $$\sgn(y) = \begin{cases}
    1 & \text{ if } y > 0 \\
    -1 & \text{ if } y < 0 \\
    0 & \text{ otherwise}
\end{cases}$$ and $x' = (\sgn(x_u))\_{u \in V}$, I claim that $c^T x' \leq c^T x$, and in particular, $x' \in \\{-1, 0, 1\\}^V$ is an integral optimal solution.

Indeed, suppose the contrary that $c^T x' > c^T x$. Denote $\delta = (\delta_u)\_{u \in V}$ where $\delta_u = x\_u - x'\_u$, then one has $c^T \delta < 0$. Consider $x'' = x - \varepsilon \delta = (x_u - \varepsilon \delta_u)\_{u \in V}$ for some $1 > \varepsilon > 0$, we have that for all edge $(u, v) \in E$,

$$
\begin{aligned}
    x''_u + x''_v 
    & = x_u + x_v - \varepsilon (\delta_u + \delta_v) \\
    &  = x_u + x_v - \varepsilon (x_u - \sgn(x_u) + x_v - \sgn(x_v)) \\
    & = (1 - \varepsilon)(x_u + x_v) + \varepsilon (\sgn(u) + \sgn(v)), 
\end{aligned}
$$

where $x_u + x_v \geq 0$ gives that both terms are non-negative. Moreover, $\|x\'\'\_u\| = (1-\varepsilon)\|x_u\|$, so $x''$ is admissible. Finally, $c^T x'' = c^T x + \varepsilon c^T \delta < c^T x$, contradicting the fact that $x$ is optimal. <span style="float:right;">$\square$</span>

Then, it suffices to solve $(2)$ with your favourite method for linear programming problems. Using either ellipsoid method or interior-point method solves $(2)$ in polynomial time; the latter gives the complexity $O(n^7)$ where $n = \|V\|$ is the number of vertices.

## 4. Last remarks

### 4.1. Generalisation to hypergraphs

Lemma 1 no longer holds for general hyper-graph. Numerical experiments show that nonetheless, if a hyper-graph has rank $r$, then there exists an optimal rational solution whose denominators divide $(r-1)!$, suggesting that it is an approximation solution of $(1)$, e.g. via some rounding technique. Pursuing this direction as done with vertex cover problem should suggest that $(1)$ is in $\FPT$ and $\APX$-complete.
		
However, this generalised version does not seem to have any relevance, thus I decided not to pursue the matter.

### 4.2. Solving $(2)$ using combinatorial algorithms

It is a curious open question about how much one can improve the complexity from $O(n^7)$. The fact that $(2)$ is totally dual integral is equivalent to saying that the corresponding polytope is integral, thus in the process of solving $(2)$, the active set in simplex method is always integral.
		
In many cases, this corresponds to a combinatorial algorithm (for example, see a brief review by Black et al. {% cite Black2021 -f on_critical_sets_of_stable_matching.bib %}).

### 4.3. Streaming algorithms

It might be of independent interest to solve $(2)$ and identify a critical set in streaming fashion. Little has been done in this regard for other problems in literature, but there is one such study done for fractional set cover by Indyk et al. {% cite Indyk2017 -f on_critical_sets_of_stable_matching.bib %}.

## References

{% bibliography --cited -f on_critical_sets_of_stable_matching.bib %}
