---
layout: post
title: "Intractability of critical sets in hypergraph matching"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\NP}{\textsf{NP}}
\newcommand{\def}{\text{def}}
\newcommand{\d}{\text{d}}
\newcommand{\APX}{\textsf{APX}}
\newcommand{\sP}{\textsf{#P}}
$$

This is a study during my second year at ENS Ulm, in the context of a research project supervised by Prof. Ana Busic.

## 1. Introduction

In my previous post ([Nguyen, 2024](on_critical_sets_in_stable_matching)), we have discussed the problem of finding the critical sets in a stochastic matching network. In case of bipartite network, knowing such critical sets can lead to an assymptotically optimal randomised policy {% cite Busic2015 -f on_critical_sets_of_stable_matching.bib %}. One can then hope that such a policy might be generalisable to hypergraphs, for which the first step is to find the crtitical sets. On the other hand, as remarked in my previous post ([Nguyen, 2024](on_critical_sets_in_stable_matching)), the failure of linear programming approach hints that perhaps the problem is computationally intractable. At the time, I did not know if the question of tractability would be relevant, and my thank for Prof. Busic for drawing my attention to its relevance.

In this post, I will show that indeed the problem of critical sets is $\NP$-hard.

## 2. Problem statement

Recall that in case of an undirected graph $G = (V, E)$ equipped with weights $c : V \rightarrow \mathbb{R}\_{\geq 0}$, for $U \subseteq V$, denote $N(U) = \\{v \mid \exists u \in U, (u, v) \in E\\} \setminus U$ be its open neighbourhood. By abuse of notation, denote $c(U) = \sum_{u \in U} c(u)$, and we wish to find
\\[
    \min_{U \subseteq V \text{ independent}} c(N(U)) - c(U).
\\]

It is straightforward to generalise this problem to hyperpgraph, as one only needs to define the notion of open neighbourhoods and indepedent sets. Let $G = (V, E)$ be a hypergraph, for $U \subseteq V$, its open neighbourhood is defined as 
\\[
    N(U) = \left\\{ \bigcup_{e \cap U \neq \emptyset} e \right\\} \setminus U,
\\]
and $U$ is said to be independent if for all $e \in E$, we have $e \subsetneq U$. The problem of critical sets is then: given a hypergraph $G = (V, E)$ with vertex weight $c : V \rightarrow \mathbb{R}$, find a set
\\[
    U^* \in \arg \min_{U \subseteq V \text{ independent}} c(N(U)) - c(U).
\\]

## 3. Reductions

We reduce this problem from minimum vertex cover. Recall the problem of minimum vertex cover: given an undirected graph $G = (v, E)$, find a set $U \subseteq V$ of minimum cardinality such that any edge is incident to at least one vertex in $U$, i.e. for all $(u, v) \in E$, either $u$ or $v$ is in $E$ (or both). As popularly known, the problem is $\NP$-hard {% cite Karp1972 -f on_critical_sets_of_stable_matching.bib %}, and moreover, the problem is $\APX$-hard {% cite Papadimitriou1991 -f on_critical_sets_of_stable_matching.bib %} even for the case where the degrees are at most 3 {% cite Alimonti2000 -f on_critical_sets_of_stable_matching.bib %}.

Now let us construct a new hypergraph as follow. Denote $n = \|V\|$. For each edge $e$ in $E$, we associate a new auxiliary vertex $u_e$. Let $G' = (V', E')$ be a hypergraph with $V' = V \cup \\{ u_e \mid e \in E \\}$ and $E' = \\{ \\{u, v, u_e \\} \mid (u, v) \in E\\}$. Lastly, we define the weight $c : V' \rightarrow \mathbb{R}\_{\geq 0}$ by $c(u) = 1$ for all $u \in V$ and $c(u_e) = 2n+1$ for all $e \in E$.

Now consider $U \subseteq$ minimising $c(N(U)) - c(U)$. Write $U_V = U \cap V$ and $U_E = U \cap (V' \setminus V)$; also write $N_V = N(U) \cap V$ and $N_E = N(U) \cap (V' \setminus V)$. We have
\\[
    c(N(U)) - c(U) = \|N_V\| - \|U_V\| + (2n+1)(\|N_E\| - \|U_E\|).
\\]

Now if $U_E \neq V' \setminus V$, let $u_{e'} \in (V' \setminus V) \setminus U_E$, and let $U' = \\{u_{e'}\\} \cup U_E$, we have

\\[
    c(N(U')) - c(U') \leq n - (2n+1)(1 + \|U_E\|), 
\\]

and

$$
\begin{align*}
    c(N(U)) - c(U) 
    & = |N_V| - |U_V| - (2n+1)(|N_E| - |U_E|) \\
    & \geq - n - (2n+1)|U_E| \\
    & > n - (2n+1)(1 + |U_E|),
\end{align*}
$$

contradicting the minimality of $c(N(U)) - c(U)$. Hence, $U_E = E$, which implies $N_E = \emptyset$, and $N(U) \subseteq V$. Moreover, we have
\\[
    c(N(U)) - c(U) = -(2n+1)\|E\| + \|N(U)\|.
\\]

Now what are the constraints on $N(U)$? Let $e = (u, v) \in E$, we know that $u_e \in U$, so $u$ and $v$ cannot both be in $U$, and both are either in $U$ or in $N(U)$, thus at least one is in $N(U)$, meaning that $N(U)$ is a vertex cover of $G$. Conversely, any vertex cover $N$ of $G$ can be used to construct an admissible set for our problem, simply by taking $U = V' \setminus N$. Thus $G$ admits a minimum vertex cover of size $k$ if and only if
\\[
    \min_{U \subseteq V \text{ independent}} c(N(U)) - c(U) = -(2n+1)\|E\| + k, 
\\]

which completes our $\NP$-reduction.

**Remark 1**. This shows a stronger statement, that our problem remains $\NP$-hard when we restrict ourselves to 3-regular hypergraphs. A small modification can be used to show that our problem remains $\NP$-hard even in the unweighted version. For $e \in E$, we simply need to replace $u_e$ of weight $2n+1$ by $2n+1$ auxiliary vertices.

## References

{% bibliography --cited -f on_critical_sets_of_stable_matching.bib %}