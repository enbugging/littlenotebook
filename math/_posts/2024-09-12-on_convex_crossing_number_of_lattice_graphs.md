---
layout: post
title: "On convex crossing number of lattice graphs"
author: "Nguyen Doan Dai"
categories: math
tags: math
---

$$
\newcommand{\crn}[1]{\text{cr}(#1)}
\newcommand{\ccrn}[1]{\text{ccr}(#1)}
\newcommand{\id}{\text{id}}
\newcommand{\var}[1]{\text{var}(#1)}
$$

This was a study during my Bachelor thesis in my first year at ENS Ulm, carried out at Hamilton Institute, Maynooth University, under supervision of Prof. Damien Woods.

For the context, it was an attempt to understand quantitatively how pseudoknotted the SST lattice - a model of DNA computation - is, i.e. what is the minimum number of crossing if one draws the polymer graph of SST lattice, which can be reduced to the problem of convex crossing number.

## 1. Introduction

Let us recall the definition of crossing number. For an undirected graph $G = (V, E)$, a drawing of $G$ is a map from the vertices $u \neq v \in V$ to distinct points $x_u \neq x_v$ on the plane, and from edges $(u, v) \in E$ to simple continuous curves $\gamma_{u, v}$ connecting their two endpoints {% cite Pach -f probabilistic_method_on_convex_crossing_number.bib %}.

For any edge $(u, v) \in E$ and any vertex $w$ different from $u$ and $v$, $\gamma_{u, v}$ must not pass through $x_w$. If two curves $\gamma_1$ and $\gamma_2$ cross at a point different from their endpoints, then the set of intersection points must be finite. The crossing number of $G$, denoted $\crn{G}$, is then the minimum, over all the drawings, of the number of crossings in a drawing. One can safely assume then, that any two edges cross at most once {% cite Pach -f probabilistic_method_on_convex_crossing_number.bib %}.

Convex crossing number, denoted $\ccrn{G}$ adds one more constraint on the drawings: that all vertices must lie on the boundary of a convex set, wherein the edges must also lie. There has been some literature regarding convex crossing number. For instance, for the $n\times n$ mesh graph $M_{n, n}$, Shahrokhi et al. {% cite Shahrokhi2003 -f probabilistic_method_on_convex_crossing_number.bib %} showed an upper bound and a lower bound and together implied that the crossing number of $M_{n, n}$ is $\Theta(n^2 \log n)$. To show the lower bound, they used techniques from graph isoperimetric inequalities.

This note has two parts. Part 1 aims at generalising this result for $M_{m, n}$. Section 2 summarises the proof that $\ccrn{M_{m, n}} = \Theta(mn \log (mn))$, and proves the same for $T_{m, n}$. To my best knowledge, neither of these questions have been addressed.

Part 2 explores the same lower bound for the convex crossing number of $M_{n, n}$ using a probabilistic method. Section 3 outlines the setting for a probabilistic method; Section 3 describes the argument (which does _not_ work, for the reasons which will be explained in Section 4). Section 5 described an alternative proof of central limit theorem for crossings of a random graph randomly embedded, and Section 6 discusses with final remarks.

## 2. Convex crossing number of $M_{m, n}$

Let us re-visit the proof that $\ccrn{M_{n, n}} = \Theta(n^2 \log n)$ and see how it can be carried over to $M_{m, n}$. The upper bound follows from the following theorem by Shahrokhi et al. {% cite Shahrokhi2003 -f probabilistic_method_on_convex_crossing_number %}.

**Theorem 1.** {% cite Shahrokhi2003 -f probabilistic_method_on_convex_crossing_number %} If $G$ is drawn in the plane with $c$ crossings, then a convex drawing of $G$ with $O\left((c + \sum\_{v \in V} d\_v^2) \log n\right)$ crossings can be constructed.

Since $M_{n, n}$ is planar, one has that $\ccrn{M_{n, n}} = O(n^2 \log n)$, and the same argument applies to $M_{m, n}$ to give that $\ccrn{M_{m, n}} = O(mn \log (mn))$. On the other hand, the lower bound for $M_{n, n}$ comes from the following therem. For a function $f : \mathbb{N} \to \mathbb{R}_{\geq 0}$, we say that a graph $G$ satisfies an $f$-edge isoperimetric inequality if for any $k$-vertex subset $U$ of $V$ and $k \leq \frac{n}{2}$, there are at least $f(k)$ edges between $U$ and $V \setminus U$. Define the difference function of $f$, denoted by $\Delta_f$, as
\\[
    \Delta_f(i) = f(i + 1) - f(i)
\\]
for any natural number $i$.

**Theorem 2.** {% cite Shahrokhi2003 -f probabilistic_method_on_convex_crossing_number %} Assume that $G = (V, E)$ satisfies an $f$-edge isoperimetric inequality so that $\Delta_f$ is positive and decreasing until $\frac{n}{2}$, then
\\[
    \ccrn{G} \geq \frac{\|V\|}{8} \sum_{j = 1}^{\lfloor \frac{\|V\|}{2} \rfloor - 1} f(j) \left(\Delta_f(j) - \Delta_f(j+1)\right) - \frac{1}{2} \|E\| f(3) - \sum_{v \in V} d_v^2.
\\]

Then the bound for $M_{n, n}$ follows from the fact that it satisfies $f$-edge isoperimetric inequality for $f(x) = \sqrt{2x}$, as shown by Theorem 3 in a paper by Bollobás and Leader {% cite Bollobs1991 -f probabilistic_method_on_convex_crossing_number %}. The last two terms are of order $O(n^2)$, so it suffices to show that 
\\[
    \sum_{j = 1}^{\lfloor \frac{n^2}{2} \rfloor - 1} f(j) \left(\Delta_f(j) - \Delta_f(j+1)\right) = \Omega( \log n).
\\]
Simple calculation gives

$$
\begin{align*}
    f(x) \left(\Delta_f(x) - \Delta_f(x+1)\right)
    & = \sqrt{2x} \left(2\sqrt{2(x+1)} - \sqrt{2x} - \sqrt{2(x+2)} \right) \\
    & = 2\sqrt{x} \left(\frac{1}{\sqrt{x} + \sqrt{x+1}} - \frac{1}{\sqrt{x+1} + \sqrt{x+2}}\right) \\
    & = 2\sqrt{x} \left(\frac{\sqrt{x+2} - \sqrt{x}}{\left(\sqrt{x} + \sqrt{x+1}\right)\left(\sqrt{x+1} + \sqrt{x+2}\right)}\right) \\
    & \geq 2\sqrt{x} \left(\frac{\sqrt{x+2} - \sqrt{x}}{\left(2\sqrt{x+2}\right)^2}\right) = \frac{1}{2(x+2)} \frac{2\sqrt{x}}{\sqrt{x+2} + \sqrt{x}}\\
    & \geq \frac{1}{2(x+2)} \sqrt{\frac{x}{x+2}} = \Omega\left(\frac{1}{x}\right),
\end{align*}
$$

and we conclude with harmonic series. Well, it turns out that much work was done on edge isoperimetric inequalities from the paper by Bollobás and Leader. In particular, Ahlswede and Bezrukov {% cite Ahlswede1995 -f probabilistic_method_on_convex_crossing_number %} showed that for $M_{m, n}$, one can also take $f(x) = \sqrt{2x}$ (although in quite cryptic language; for a paraphrasing, see Lemma 3.4 in a paper by Rolim et al. {% cite Rolim1995 -f probabilistic_method_on_convex_crossing_number %}). Theorem 2 applies the same way, and we obtain that
\\[
    \ccrn{M_{m, n}} = \Theta(mn \log (mn)),
\\]
which resolves an open question by Schaefer {% cite Schaefer2013 -f probabilistic_method_on_convex_crossing_number %}.

**Remark 1.** Shahrokhi et al. {% cite Shahrokhi2004 -f probabilistic_method_on_convex_crossing_number %} also showed that if one takes $H$ to be one of the infinite square, hexagonal, or triangular lattices in the place and let $C$ be any bounded open convex domain, then for any $\lambda > 1$, denoting $H_\lambda$ the induced subgraph of $H$ on the vertex set $V(H) \cap \lambda C$, one has 
\\[
    \ccrn{H_\lambda} = \Omega\left(\|V(H_\lambda)\| \cdot \log \| V(H_\lambda) \|\right),
\\]
where $\Omega(\cdot)$ refers to the asymptotics as $\lambda$ tends to infinity. However, it is not clear how it will be translated to an asymptotics for $\Omega(M_{m, n})$.

**Remark 2.** Consider the torus graph $T_{m, n}$, whereof $M_{m, n}$ is a subgraph, thus one has that $\ccrn{T_{m, n}} \geq \ccrn{M_{m, n}}$, so the same lower bound applies. Now the upper bound no longer applies, or at least not trivially so; Kainen {% cite Kainen2011 -f probabilistic_method_on_convex_crossing_number %} showed that
\\[
    \ccrn{T_{m, n}} \leq \frac{m^3 + 3n^2}{2}
\\]
for even $m$ and $n$. However, there are generalisations of the Theorem 1 above, which are applicable in this case. For example, Pach and Tóth {% cite Pach2006 -f probabilistic_method_on_convex_crossing_number %} showed that if a family of graph $(G\_n)\_{n \in \mathbb{N}}$ whose member has genus at most $g$, then one has
\\[
    \crn{G} \leq c_g \sum_{v \in V} d_v^2
\\]
for some constant $c_g$ depending only on $g$. In particular, for $T_{m, n}$, one can take $g = 1$, thus giving $\crn{T_{m, n}} = \Omega(mn)$. Theorem 1 then implies that $\ccrn{T_{m, n}} = O(mn \log (mn))$, thus 
\\[
    \ccrn{T_{m, n}} = \Theta(mn \log (mn)).
\\]

## 3. Preliminary for a probabilistic approach

Consider $n\times n$ grid $G = M_{n, n} = (V, E)$, with $\|V\| = n^2$ vertices and $\|E\| = 2n^2 - 2n$ edges. Let $\sigma \in S_{n^2}$ be a permutation, and we arrangement the vertices of $V$ on the line according to $\sigma$, i.e. $\sigma(0)$, then $\sigma(1)$, and so on. Denote $c_n = \ccrn{M_{n, n}}$ its convex crossing number.

Recall that by taking $\sigma = \id_{\\{0, \cdots, n^2-1\\}}$, we obtain that $c_n = O(n^3)$. To obtain an elementary non-trivial lower bound, consider a diagram of $G$ with $c_n$ crossings. Each of them can be removed by removing an edge from $G$, thus after deleting at most $c_n$ edges, we obtain a planar graph with at least $2n^2 - 2n - c_n$ edges and $n^2$ vertices. From Euler's formula, one obtain an inequality that for any planar graph with $v$ vertices and $e$ edges must have $e \leq 3v$. In this case, we have $2n^2 - 2n - cn \leq 3n^2$, or $c_n \geq n^2 - 2n$.

Recall furthermore that there are $\frac{\|E\|(\|E\|-1)}{2}$ pairs of edges, but not all pairs can intersect. In particular, pairs of edges sharing one endpoint cannot cross, and all the other pairs of edges can. Instead, there are exactly
\\[
    \mathcal{E}\_n = \frac{\|E\|(\|E\|-1)}{2} - \frac{1}{2}\sum\_{v \in V} d\_v^2 = \Theta(n^4)
\\]
pairs of vertices that can cross in any drawings of $G$, where for $v \in V$, $d_v$ denotes the degree of $v$. 

Choosing $\sigma$ from $S_{n^2}$ uniformly random, and let $X\_1, X\_2, \cdots, X\_{\mathcal{E}\_n}$ be $\mathcal{E}\_n$ random variables, where $X\_i = 1$ if the $i^{\text{th}}$ pair of edges cross in a drawing of $G$ that has $c(\sigma)$ crossings, and $0$ otherwise. Let $X = \sum\_{i = 1}^{\mathcal{E}\_n} X_i$.

First, let us find $\mathbb{E}[X]$. Consider a pair of edges that can cross, $(a, b)$ and $(c, d)$. Since the crossing of these edges involves only themselves and no other edges, we can safely ignore other vertices. There are $24$ permutations of $a$, $b$, $c$, and $d$, but since we place them on a circle, we can fix the first vertex, say, $a$. Then there remain $6$ permutations, which are $(a, b, c, d)$, $(a, b, d, c)$, $(a, c, b, d)$, $(a, c, d, b)$, $(a, d, b, c)$, $(a, d, c, b)$. Only $2$ of which gives a crossing, thus $\mathbb{E}[X_i] = \frac{1}{3}$, and so $\mathbb{E}[X] = \mu_n = \frac{\mathcal{E}_n}{3}$.

## 4. Dependency graph and tail bound

**Definition 1.** Let $Y = \\{Y_{n, k} \mid 1 \leq k \leq N_n\\}$ be an array of random variables. We say that a graph $G$ is a _dependency_ graph for $Y$ if the following two conditions are satisfied:
- There exists a bijection between the random variables $Y_{n, k}$ and the vertices of $G$.
- If $V_1$ and $V_2$ are two disjoint sets of vertices of $G$ such that no edge of $G$ joins a vertex in $V_1$ and another in $V_2$, then the corresponding sets of random variables are independent.

In our case, consider a graph $\Gamma_n$ with $\mathcal{E}\_n$ vertices, and there is an undirected simple edge between vertex $i$ and $j$ ($i \neq j$) if and only if the $i^{\text{th}}$ and $j^{\text{th}}$ pair of edges share exactly one edge. Following our argument above, one see that $\Gamma_n$ is a dependency graph. Moreover, it is easy to see that $\Gamma_n$ is a subgraph of $K_{\|E\|} \square K_{\|E\|}$, which has the chromatic number $\|E\|$, so $\chi(\Gamma_n) = \chi_n \leq \|E\| = \Theta(n^2)$.

The crucial ingredient in this argument is a quantitative version of the famous Janson's dependency criterion.

**Theorem 3.** (Simplification of equation (2.11), Corollary 2.6, {% cite Janson2004 -f probabilistic_method_on_convex_crossing_number.bib %}) Let $Y_1, Y_2, \cdots, Y_n$ admit $\Gamma$ as a dependency graph with chromatic number $\chi$, and assume that $Y_i \sim Ber(p)$. Denote $Y = \sum_{i = 1}^n Y_i$, then
\\[
    \mathbb{P}(\mathbb{E}[Y] - Y \geq t) \leq \exp\left(-\frac{t^2\left(1 - \frac{\chi}{4n}\right)}{2\chi\cdot n \cdot p}\right).
\\]

In our case, this gives, 
\\[
    \mathbb{P}(X \leq c_n) = \mathbb{P}(\mu_n - X \geq \mu_n - c_n) \leq \exp\left(-\frac{(\mu_n - c_n)^2\left(1 - \frac{\chi_n}{4\mathcal{E}_n}\right)}{2\chi_n \cdot \mathcal{E}_n \cdot \frac{1}{3}}\right).
\\]

On the other hand, each vertex has at most 4 neighbours, so there are at most $5! = 120$ possible arrangement of that vertex and its neighbours. That gives $C_1^{n^2}$ collections of orders at local scale, each of which, if compatible, will determine one unique arrangement of vertices (note that $C_1 < 120$ since not all collections of local orders are compatible and thus they will not give arrangements). But, since we are arranging the vertices around the circles, each arrangement is over-counted $n^2$ times due to cyclic permutations, so there are at most $\frac{C_1^{n^2}}{n^2}$ configurations.

This means the probability that a drawing has minimal convex crossing number is at least $\frac{n^2}{C_1^{n_2}}$, implying
\\[
    \exp\left(-\frac{(\mu_n - c_n)^2\left(1 - \frac{\chi_n}{4\mathcal{E}_n}\right)}{2\chi_n \cdot \mathcal{E}_n \cdot \frac{1}{3}}\right) \geq \frac{n^2}{C_1^{n_2}}.
\\]
Taking logarithm gives
\\[
    \frac{(\mu_n - c_n)^2\left(1 - \frac{\chi_n}{4\mathcal{E}_n}\right)}{2\chi_n \cdot \mathcal{E}_n \cdot \frac{1}{3}} \leq n^2 \ln C_1 - 2\ln n.
\\]
Now note that $\left(n - \frac{\ln n}{n}\right)^2 = n^2 - 2\ln n + \frac{\ln^2 n}{n^2}$, so we can say that $n^2 \ln C_1 - 2\ln n \geq C_2 (n - \frac{\ln n}{n})^2$ for some $C_2 > 0$. On the other hand, $\chi_n \mathcal{E}_n = \Theta(n^6)$ and $1 - \frac{\chi_n}{4\mathcal{E}_n} = 1 - \Theta\left(\frac{1}{n^2}\right)$, so in summary, we can write
\\[
    \frac{(\mu_n - c_n)^2}{n^6} \leq C_3^2 \left(n - \frac{\ln n}{n}\right)^2 \Leftrightarrow \frac{\mu_n - c_n}{n^3} \leq C_3 \left(n - \frac{\ln n}{n}\right),
\\]
for some $C_3 > 0$.

Now, note that $\mu_n \sim \frac{1}{3} \mathcal{E}_n \sim \frac{2}{3} n^4$, so **if** $C_3 < \frac{2}{3}$, then the only way that this inequality is true, is that for $n$ large enough, 
\\[
    \frac{c_n}{n^2} - C_3 \frac{\ln n}{n} \geq \frac{\mu_n - C_3 n^4}{n^3} \geq 0,
\\]
which gives $c_n = \Omega(n^2 \ln n)$, as desired.

## 5. But does it work?

Let us hope for the following stronger inequality, since the term $1 - \frac{\chi_n}{4\mathcal{E}_n}$ can be seen as almost 1 for $n$ very big.

$$
\begin{align*}
    (\mu_n-c_n)^2 
    & \leq \frac{2}{3} \chi_n \mathcal{E}_n \left(n^2 \ln C_1 - 2\ln n\right) \\
    & \leq \frac{2}{3} \chi_n \mathcal{E}_n \left(n \sqrt{\ln C_1} - \frac{\ln n}{\sqrt{\ln C_1}}\right)^2 \\
    \Leftrightarrow \mu_n - c_n
    & \leq \frac{\sqrt{6}}{3} \sqrt{\chi_n \mathcal{E}_n} \left(n \sqrt{\ln C_1} - \frac{\ln n}{\sqrt{\ln C_1}}\right) \\
\end{align*}
$$

For the argument in section 3 to work, we would need 
\\[
    \mu_n \geq \frac{\sqrt{6}}{3} n \sqrt{\chi_n \mathcal{E}_n \ln C_1}.
\\]
Now, suppose we do have $\chi_n \sim \|E\| \sim 2n^2$ (which should be true, but I will not include further justification), then as $\mathcal{E} \sim 2n^4$, 
\\[
    \frac{\sqrt{6}}{3} n \sqrt{\chi_n \mathcal{E}_n \ln C_1} \sim n^4 \cdot \frac{2\sqrt{6\ln C_1}}{3},
\\]
so we would need $\frac{2\sqrt{6\ln C_1}}{3} < \frac{2}{3}$, i.e $\ln C_1 < \frac{1}{6}$, or $C_1 < e^{\frac{1}{6}} \approx 1.18136041...$. This provides a heuristic why this line of argument most likely will not work: one would not expect $C_1$ to be so small.

Can the dependence of the coefficients be improved? Janson mentioned in his paper {% cite Janson1988 -f probabilistic_method_on_convex_crossing_number.bib %} that there is little room left for improvement, so this crude bound, unfortunately, seems optimal. Perhaps more refined analysis is possible with the machinery of weighted dependency graph developed by Féray {% cite Feray2016 -f probabilistic_method_on_convex_crossing_number.bib %}, but it also appears to me that we have used much, if not all, of information available from the hypotheses.

## 6. Limiting distribution

It is known that the distribution of $X$ is asymptotically normal, as demonstrated by central limit theorem of Arenas-Velilla et al. {% cite Arenas-Velilla2023 -f probabilistic_method_on_convex_crossing_number.bib %}. In their work, Arenas-Velilla and colleagues used Stein's method to prove an explicit upper bound on the Kolmogorov distance to normal distribution. Here we will show another argument using Janson's normality criterion. It is only capable of demonstrating convergence in distribution, but I include it nonetheless.

**Theorem 4.** {% cite Janson1988 -f probabilistic_method_on_convex_crossing_number.bib %} Let $Y = \\{Y_{n, k} \| 1 \leq k \leq N_n\\}$ be an array of random variables such that for all $n$ and for all $1 \leq k \leq N_n$, the inequality $\|Y_{n, k}\| \leq A_n$ holds for some real number $A_n$, and that the maximum degree of a dependency graph is $\chi_n$.
    
Set $Y_n = \sum_{i = 1}^{N_k} Y_{n, k}$, $\mu_n = \mathbb{E}[Y_n]$, and $\sigma_n^2 = \var{Y_n}$. If there exists a natural number $\ell$ such that
\\[
    N_n \chi_n^{\ell - 1} \left(\frac{A_n}{\sigma_n}\right)^\ell \xrightarrow{n \to \infty} 0, 
\\]
then $Y_n \xrightarrow{} \mathcal{N}(\mu_n, \sigma_n^2)$, where the convergence is in distribution.

For our case, one can take $A_n = 1$, then for all $\ell > 2$,
\\[
    N_n \chi_n^{\ell - 1} \left(\frac{A_n}{\sigma_n}\right)^\ell = O(n^{4 + 2(\ell-1) - 3\ell}) = O(n^{2 - \ell}) = o(1).
\\]

## 7. Last remarks

It is thus a curious question whether this approach can be made to work. If so, it would appear to prove a non-trivial lower bound on the convex crossing number: that a graph of $n$ vertices and $m$ edges, under some right conditions (such as $2$-vertex-connected and all vertices have degree at least 3), will have the crossing number of order $\Omega(m\ln n)$. This will match all known lower and upper bounds.

On the other hand, even if this line of argument proves to be working in some cases, there are inherent limits to it. For instance, it is an unforgivable mistake to assume that only because of convergence in distribution, one can use normal distribution as the true distribution in calculation. For instance, if one were to believe in normal distribution, then by symmetry, one would reach the conclusion that the _maximum_ number of convex crossing for a grid is $\frac{2}{3} \mathcal{E}_n - O(n^2 \ln n)$. This is false, as demonstrated by the following two results that such the number is $\mathcal{E}_n - O(n^3)$.

**Theorem 5.** {% cite Chimani2018 -f probabilistic_method_on_convex_crossing_number.bib %} Let $G$ be a bipartite graph with $M$ pairs of independent edges and bipartite crossing number $k$. Then, the maximum number of convex crossings is $M - k$.

**Theorem 6.** {% cite Shahrokhi2000 -f probabilistic_method_on_convex_crossing_number.bib %} The bipartite crossing number of $M_{n, n}$ is $O(n^3)$.

### References

{% bibliography --cited -f probabilistic_method_on_convex_crossing_number.bib %}