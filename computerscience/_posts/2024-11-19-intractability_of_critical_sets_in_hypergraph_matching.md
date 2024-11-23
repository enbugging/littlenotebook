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
\newcommand{\L}{\textsf{L}}
\newcommand{\sP}{\textsf{#P}}
$$

This is a study during my second year at ENS Ulm, in the context of a research project supervised by Prof. Ana Busic.

## 1. Introduction

In my previous post ([Nguyen, 2024](on_critical_sets_in_stable_matching)), we have discussed the problem of finding the critical sets in a stochastic matching network. In case of bipartite network, knowing such critical sets can lead to an assymptotically optimal randomised policy {% cite Busic2015 -f on_critical_sets_of_stable_matching.bib %}. One can then hope that such a policy might be generalisable to hypergraphs, for which the first step is to find the crtitical sets. On the other hand, as remarked in my previous post ([Nguyen, 2024](on_critical_sets_in_stable_matching)), the failure of linear programming approach hints that perhaps the problem is computationally intractable. At the time, I did not know if the question of tractability would be relevant, and my thank for Prof. Busic for drawing my attention to its relevance.

In this post, I will show that indeed the problem of critical sets is $\NP$-hard and $\APX$-hard.

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

## 3. Polynomial-time reduction and $\NP$-hardness.

We reduce this problem from minimum vertex cover. Recall the problem of minimum vertex cover: given an undirected graph $G = (v, E)$, find a set $U \subseteq V$ of minimum cardinality such that any edge is incident to at least one vertex in $U$, i.e. for all $(u, v) \in E$, either $u$ or $v$ is in $E$ (or both). As popularly known, the problem is $\NP$-hard {% cite Karp1972 -f on_critical_sets_of_stable_matching.bib %}.

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

Now what are the constraints on $N(U)$? Let $e = (u, v) \in E$, we know that $u_e \in U$, so $u$ and $v$ cannot both be in $U$, and both are either in $U$ or in $N(U)$, thus at least one is in $N(U)$, meaning that $N(U)$ is a vertex cover of $G$. Conversely, any vertex cover $N$ of $G$ can be used to construct an admissible set for our problem, simply by taking $U = V' \setminus N$. Thus we have
\\[
    \min_{U \subseteq V \text{ independent}} c(N(U)) - c(U) = -(2n+1)\|E\| + \tau(G), 
\\]

which completes our polynomial time reduction.

**Remark 1**. This shows a stronger statement, that our problem remains $\NP$-hard when we restrict ourselves to 3-regular hypergraphs. A small modification can be used to show that our problem remains $\NP$-hard even in the unweighted version. For $e \in E$, we simply need to replace $u_e$ of weight $2n+1$ by $2n+1$ auxiliary vertices.

## 4. $\L$-reduction and $\APX$-hardness

$\APX$-hardness implies $\NP$-hardness, but since our $\L$-reduction is somewhat trickier and more technically heavy, I chose to include the polynomial-time reduction above nonetheless.

Our $\L$-reduction also starts with minimum vertex cover, which is known to be $\APX$-hard {% cite Papadimitriou1991 -f on_critical_sets_of_stable_matching.bib %} even for the case where the degrees are at most 3 {% cite Alimonti2000 -f on_critical_sets_of_stable_matching.bib %}.

Given a cubic connected graph $G = (V, E)$. Assume without loss of generality that $n = \|V\| \geq 2$. We construct the following auxiliary graph $G' = (V', E')$. $V'$ is obtained by combining two copies of $V$, namely $V$ and $\overline{V} = \\{ \overline{u} \mid u \in V\\}$, i.e. each vertex $u \in V$ has another copy $\overline{u} \in \overline{V}$. Then, for each edge $(u, v) \in E$, we have two hyper-edges in $E'$, namely $\\{u, v, \overline{u}\\}$ and $\\{u, v, \overline{v}\\}$. We assign to each vertex the following weight $c$: $c(u) = 1$ and $c(\overline{u}) = \lambda$ for all $u \in V$, where $\lambda$ is a constant to be determined later.

Let $U \subseteq V'$ be minimising $c(N(U)) - c(U)$. As usual, denote $U_V = U \cap V$, and $\overline{U_V} = U \cap \overline{V}$, and $N_V = N(U) \cap V$ and $\overline{N_V} = N(U) \cap \overline{V}$. We have
\\[
    c(N(U)) - c(U) = \|N_V\| - \|U_V\| + \lambda\left(\| \overline{N_V} \| - \| \overline{U_V} \|\right)
\\]

First suppose $\overline{N_V} \cup \overline{U_V} \neq \overline{V}$, let $\overline{u} \in \overline{V} \setminus (\overline{N_V} \cup \overline{U_V})$, and consider $U' = U \cup \\{\overline{u'}\\}$. Then $N(U') = N(U) \cup N(\\{\overline{u}\\})$. Since $G$ is cubic, $\|N(\\{\overline{u}\\})\| = 4$, and thus

$$
\begin{align*}
    c(N(U')) - c(U') 
    & = c(N(U) \cup N(\{\overline{u}\})) - c(U \cup \{u\}) \\
    & \leq c(N(U)) + c(N(\{\overline{u}\})) - c(U) - c(\{u\}) \\
    & \leq c(N(U)) - c(U) + 4 - \lambda \\
    & < c(N(U)) - c(U),
\end{align*}
$$

which will contradict the minimality of $U$ if $\lambda > 4$.

Then we know that $\overline{N_V} \cup \overline{U_V} = \overline{V}$, suppose that $\overline{N_V} \neq \emptyset$. Let $\overline{u} \in \overline{N_V}$ and consider $U' = (U \setminus \\{u\\}) \cup \\{\overline{u}\\}$ which is still an independent set.

Now this is where the analysis gets most tricky. To simplify our lives, we have the following lemma.

**Lemma 1.** For $A, B \subseteq V'$, one has $N(A \cup B) \subseteq N(A) \cup N(B)$.
_Proof._ Standard set manipulation gives

$$
\begin{align*}
    N(A \cup B)
    & = \left(\bigcup_{e \cap (A \cup B) \neq \emptyset} e \right) \setminus (A \cup B) \\
    & = \left(\bigcup_{(e \cap A) \cup (e \cap B) \neq \emptyset} e \right) \setminus (A \cup B) \\
    & = \left[\left(\bigcup_{e \cap A \neq \emptyset} e\right) \cup \left(\bigcup_{f \cap B \neq \emptyset} f\right) \right] \setminus (A \cup B) \\
    & = \left[\left(\bigcup_{e \cap A \neq \emptyset} e\right) \setminus (A \cup B) \right] \cup \left[ \left(\bigcup_{f \cap B \neq \emptyset} f\right) \setminus (A \cup B) \right] \\
    & \subseteq \left[\left(\bigcup_{e \cap A \neq \emptyset} e\right) \setminus A \right] \cup \left[ \left(\bigcup_{f \cap B \neq \emptyset} f\right) \setminus B \right] \\
    & = N(A) \cup N(B),
\end{align*}
$$

as desired. <span style="float:right;">$\square$</span>

Thus $N(U') \subseteq N(U) \cup N(\\{\overline{u}\\})$. Now note that since $G$ is cubic, $\|N(\\{\overline{u}\\}) \| = 4$, and by construction, $N(\\{\overline{u}\\}) \subset V$, hence

$$
\begin{align*}
    c(N(U')) - c(U')
    & \leq c(N(U)) + c(N(\{\overline{u}\})) - c(U \setminus \{u\} \cup \{\overline{u}\}) \\
    & = c(N(U)) + 4 - c(U) + c(\{u\}) - c(\{\overline{u}\}) \\
    & = c(N(U)) + 4 - c(U) + 1 - \lambda \\ 
    & < c(N(U)) - c(U), 
\end{align*}
$$

which also contradicts the minimality of $U$ if $\lambda > 5$.

Then we know that $\overline{U_V} = \overline{V}$, implying $N(U) \subseteq V$. Similarly as in the previous section, consider $(u, v) \in E$, we have that $\overline{u}, \overline{v} \in U$, so $u$ and $v$ cannot be both in $U$, i.e. at least one of them is in $N(U)$, meaning $N(U)$ is a vertex cover of $G$. Moreover, we now have that $U_V = V \setminus N_V = V \setminus N(U)$, thus
\\[
    c(N(U)) - c(U) = \|N(U)\| - \lambda (|V| + \|V \setminus N(U)\|) = (1 + \lambda) \|N(U)\| - 2\lambda \| V \|.
\\]
Since $U$ minimises $c(N(U)) - c(U)$, necessarily we have $\|N(U)\| = \tau(G)$. Again, this reasoning works for any $\lambda > 4$.

Next, we want to bound $\tau(G)$ in terms of $n$. There are many, many ways to do that, but here we include a folklore argument for simplicity: $G$ is cubic, thus each vertex in $N(U)$ has degree $3$, and hence there are at most $3 \tau(G)$ edges incident to $N(U)$. Howerver, $N(U)$ is a vertex cover, thus all edges of $G$ are incident to $N(U)$, which implies $|E| \leq 3 \tau(G)$, or
\\[
    \tau(G) \geq \frac{|E|}{3} \geq \frac{n-1}{3} \geq \frac{n}{6}, 
\\]
since $G$ is connected and $n \geq 2$.

Therefore, 
\\[
    c(N(U)) - c(U) \leq -(\lambda+1)\frac{n}{6} - 2\lambda n = -\frac{13\lambda+1}{6} n < 0.
\\]
Let $0 \leq \varepsilon \leq \min\left(1, \frac{1}{2C}\right)$ for some constant $C > 0$ to be chosen later, and $U' \subseteq V$ be some subset such that $c(N(U')) - c(U') \leq (1 - C \varepsilon) \[c(N(U)) - c(U)\]$, then we have
\\[
    c(N(U')) - c(U') \leq -(1 - C \varepsilon)\frac{13\lambda+1}{6} n, 
\\] 
which implies that in fact $\overline{V} \subseteq U'$. Indeed, $\frac{13\lambda+1}{6} > 2\lambda$, which implies 
\\[
    c(N(U')) - c(U') \leq -(1 - C \varepsilon)\frac{13\lambda+1}{6} n < -\lambda n. 
\\]
On the other hand, if there exists some $\overline{u} \in \overline{V} \setminus U'$, then $U' \subseteq (V' \setminus \\{\overline{u}\\})$ $ = V \cup (\overline{V'} \setminus \\{u\\})$, which implies, 
\\[
    c(N(U')) - c(U') \geq -c(U') \geq -(\lambda-1) n - n = -\lambda n.
\\]
And once we have $\overline{V} \subseteq U'$, an argument as above shows that $N' = N(U')$ is a vertex cover, thus we have

$$
\begin{align*}
    c(N(U')) - c(U') 
    & \leq (1-C \varepsilon) [c(N(U)) - c(U)] \\
    \Leftrightarrow (\lambda+1) |N'| - 2\lambda n 
    & \leq (1-C \varepsilon) [(\lambda + 1) \tau(G) - 2\lambda n] \\
    \Leftrightarrow (\lambda+1) |N'| 
    & \leq (1-C \varepsilon) (\lambda + 1) \tau(G) - 2\lambda C \varepsilon n \\
    \Leftrightarrow |N'| 
    & \leq (1-C \varepsilon) \tau(G) + \frac{2\lambda}{\lambda + 1} C \varepsilon n \\
    & \leq (1-C \varepsilon) \tau(G) + \frac{6\lambda}{\lambda + 1} C \varepsilon \tau(G) \\
    \Leftrightarrow \frac{|N'|}{\tau(G)} & \leq 1 + \epsilon C \frac{5\lambda - 1}{\lambda}.
\end{align*}
$$

Thus to recap, if we choose $\lambda > 5$ and $C \leq \frac{1}{5} < \frac{\lambda}{5\lambda - 1}$, then $\frac{\|N'\|}{\tau(G)} \leq 1 + \varepsilon$, and our $\L$-reduction is complete.



## References

{% bibliography --cited -f on_critical_sets_of_stable_matching.bib %}