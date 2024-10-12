---
layout: post
title: "Introduction to RNA polytopes"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\conv}{\text{conv}}
\newcommand{\vertex}{\text{vert}}
$$

This is mainly taken from Chapter 1, 2, and 3 of my Bachelor thesis at École Polytechnique under supervision of Prof. Sarah J. Berkemer and Prof. Yann Ponty, titled ["Polynomial-time parametric optimisation"](https://arxiv.org/abs/2304.14962)

## 1. Introduction

### 1.1. Importance of RNA
Since the discovery of ribosomal RNA and transfer RNA in the late 1950s, a variety of non-coding RNA families have been discovered, which, to name a few, includes long non-coding RNAs {% cite Statello2021 -f rna_polytopes %}, circular RNAs {% cite Liu2022 -f rna_polytopes %}, and antisense RNAs and CRISP RNAs {% cite Ashwath2023 -f rna_polytopes %}, each with potential for new therapies. Non-coding RNAs have shown to act as catalysts for chemical reactions, as gene regulators {%cite Statello2021 -f rna_polytopes -f rna_polytopes %}, and as DNA replicator {% cite Kuhnlein2021 -f rna_polytopes %}, thereby are believed to play a vital role in the origin of life. We now know that at least $76\%$ of human genome is transcribed, yet only 1.2% of which encodes proteins {% cite Lee2019 -f rna_polytopes %}, and most of the rest, i.e. non-coding RNAs, have their functions yet to be completed unveiled, and applications to be discovered. Thus, there raises the need to understand properties and functions of RNAs.

### 1.2. Importance of computational method
On the one hand, advances in sequencing techniques have resulted in an exponential growth for the total length of sequences for the past decade {% cite Katz2022 -f rna_polytopes %}. On the other hand, inferring structure of RNA by experimental methods, despite much progress in recent years, remains both arduous, time-consuming, and technically demanding. There have been attempts to combine experiment data with simulation, yet this approach has yet to overcome the limitation in both accuracy and throughput {% cite Liu2022 -f rna_polytopes %}. Although RNA has its tertiary structure largely determined by its secondary structure {% cite Tinoco1990 --file rna_polytopes %}, the same challenges persist, making prediction by computational methods being most viable approach, amongst which is free energy minimisation method.

### 1.3. Importance of energy model
Whilst delaying formal definition to later sections, we recall that an RNA sequence consists of nucleotides, commonly denoted by four letters $A$, $C$, $G$, and $U$. A free energy model of RNA then admits, amongst others, the following assumptions:
- An RNA exists as it is and folds as a whole, contrasting with co-transcription model where RNA folds as it is being transcribed.
- An RNA admits only the canonical base pairs $A-U$ and $G-C$, and the so-called wobble pairs $G-U$.
- An RNA folds in isolation, independent from the surrounding environment. 

A general thermodynamic framework consists of a model which assigns to a given pair of sequence and secondary structure a measure, e.g. energy, entropy, and/or enthalpy, with respect to the state where no nucleotides are paired. Then, under some assumptions such as pseudoknot-free and heuristics regarding the energy function, it is possible to decompose the total energy as those of substructures, whereby a dynamic programming algorithm can exploit the sub-problem hierarchy to minimise or maximise this energy, before tracing back to find an optimal solution efficiently. One may see that whilst the performance bottleneck may lie in the algorithm, the limit of prediction capability lies in the underlying energy model, whence comes the need for a biologically accurate energy model.

At the one end of the spectrum is the model counting number of base pairs as featured in work of Nussinov and Jacobson {% cite Nussinov1980a -f rna_polytopes %}, whose simplicity allows full analysis of the parameter space but severely limits the prediction capacity. At the other end for pseudoknot-free secondary structures is Turner energy model {% cite Mathews1999a -f rna_polytopes %} with over 7600 features, generally considered to be biologically realistic, but at the same time difficult to well-tune. A subset of parameters were measured by optical melting experiments, but a large part was derived by fitting to the experiment data.

Despite its comprehensiveness and even regarding pseudoknot-free RNA secondary structures, Turner model fails to predict accurately in many cases. Additionally, the ad-hoc energy function for multi-loops, originally derived for simplicity and algorithm derivation's sake, outperforms other more complicated and realistic alternatives {% cite Ward2017 -f rna_polytopes %}. These two phenomena beg the question if it is due to a suboptimal choice of parameters or the fundamental limit of the model itself. This line of study is often called parametric analysis, where one explores the parameters to observe what a model can predict and argue about its capability, and is of particular importance, for this model is also applied to the study of single-stranded DNA {% cite SantaLucia2004 -f rna_polytopes %}, which has both an essential role in virology {% cite Malathi2019 -f rna_polytopes %} and potential in therapeutics {% cite Avci-Adali2015 -f rna_polytopes %}. 

## 2. Notations and preliminary

In what follows, let $X = \mathbb{R}^{d+1}$ or $\mathbb{S}^d$ for $d > 1$. We shall consider the Euclidean space $\mathbb{R}^{d+1}$ equipped with the canonical inner product $\langle \cdot, \cdot \rangle$. For $m \leq d$, we define a $m$-plane to be a $m$-dimensional subspace of $\mathbb{R}^{d+1}$. Unless explicitly used for other purposes, uppercases letters represent polytopes, polyhedra, and hyperplanes, whilst lowercase letters represent affine points as column vectors, thus for $x, y \in \mathbb{R}^{d+1}$, we have $\langle x, y \rangle = x^T y$. We denote $[x, y]_{\mathbb{R}^n}$ as the line segment uniquely given by two points $x, y \in \mathbb{R}^{d+1}$, where the subscript shall be neglected whenever the context is apparent.

Except what is defined here and with possibly different notations, the definitions may be found in Ziegler's _Lectures on Polytopes_ {% cite Ziegler1995 --file rna_polytopes %} and in Ratcliffe's _Foundations of Hyperbolic Manifolds_ {% cite Ratcliffe2007 --file rna_polytopes %}.

### 2.1. RNA sequences and secondary structures

An RNA sequence $q = q\_1 q\_2 ... q\_n$ is a string of nucleotides, each of which is one of the four letters (bases) $A$, $C$, $G$, and $U$. Then, given such a sequence, a base pair can form between two distinct positions $i < j$, and denoted by the unordered pair $\\{i, j\\}$. Hereinafter, we shall admit only the canonical base pairs and the wobble pairs, i.e. one must have $\\{s\_i, s\_j\\} = \\{A, U\\}$, $\\{G, C\\}$, or $\\{G, U\\}$. However, each base can belonged to at most one base pairs.

A set of $k$ such unordered pairs $s = \\{\\{i\_\ell, j\_\ell\\} \mid \ell=1,2,...k\\}$ is called a secondary structure, and said to be compatible if the base pairing are all admitted. Two base pairs $\\{i, j\\}$ and $\\{k, l\\}$ such that $i < k < j < l$ or $k < i < l < j$ are said to be crossing, and thus form what is called a pseudo-knot, and shall not be concerned in this thesis. Given an RNA sequence $q$ and a compatible secondary structure $s$, this assumption allows one to decompose a secondary structure into various smaller structures of some types, and thus may use dynamic programming to calculate certain features of the pair $(q, s)$. For instance, we may count the number $n\_{ij}$ of base pairs of type $i-j$, which can be summarised by a triplet $c(q, s) = (c\_{AU}, c\_{GC}, c\_{Gu}) \in \mathbb{N}^3$, which we call a feature vector or a signature.

Continue with this example, we may impose an energy model over these features, by choosing a parameter set $p = (p\_{AU}, p\_{GC}, p\_{GU}) \in \mathbb{R}^3$. Then, for a pair $(q, s)$, we have the energy to be $E = \langle c(q, s), p \rangle = c\_{AU} p\_{AU} + c\_{GC} p\_{GC} + c\_{GU} p\_{GU}$. And moreover, given $q$ and $p$, one can find a compatible structure $s$ minimising $E$ by Nussinov's algorithm {% cite Nussinov1980a --file rna_polytopes %}: in particular, let $f(i, j)$ be the minimum energy over the subsequence $q' = q\_i q\_{i+1} \cdots q\_{j-1} q_j$, one has
<div align="center">
$$
	f(i, j) = \max
	\begin{cases}
		f(i+1, j) \\\
		f(i, j-1) \\\
		f(i-1, j+1) + p_{\{q_i, q_j\}} \text{ if } q_i \text{ and } q_j \text{ can form a base pair} \\\
		\max_{i < k < j} f(i, k) + f(k+1, j)
	\end{cases}.
$$
</div>
together with initialisation $f(i, i) = f(i, i-1) = 0$ for all $i$. Then, given the memoisation table $f$, one can trace back to find an optimal and compatible secondary structure $s$.

Note that this formulation of Nussinov's algorithm has inherent ambiguity, and in particular, a secondary structure can be traced back in different ways, but this does not change the final output and thus of our concern. Likewise, one signature can correspond with multiple secondary structures, and depending on the choice of parameters, it is possible that multiple signatures are optimal.

### 2.2. Convex set
A set $K \subseteq X$ is convex if for any two points $x, y \in K$, one has $[x, y] \in K$. It is clear that the intersection of two convex sets are convex, and $X$ is convex, hence any set $S \subseteq X$ admits a minimal convex set containing it, called the convex hull of $S$ and denoted as $\conv{S}$.
\\[
    \conv{S} = \bigcap \left\\{K' \mid S \subseteq K' \subseteq X, K' \text{convex}\right\\}.
\\]

Now we consider $X = \mathbb{R}^{d+1}$, $K \subseteq \mathbb{R}^{d+1}$ compact and convex, and for $y \in \mathbb{R}^{d+1}$, we define the supporting function of $K$ at $y$ as
\\[
    h_K(y) = \sup_{x \in K} \langle x, y \rangle, 
\\]
where one can show the maximum is attained on the boundary of $K$. For computational purpose, we define an extremal function $\sigma_K$ as an oracle returning such a vector $x \in \partial K$ at which the minimum is attained, i.e. $\langle \sigma_K(y), y \rangle = h_K(y)$, chosen arbitrarily if there exists many.
Let $A$ and $B$ be two convex sets, we define Minkowski sum of $A$ and $B$ as $A + B = \{a + b \mid (a, b) \in A \times B\}$. It can be proven that $A + B$ is convex, and if $A$ and $B$ are compact, then so is $A + B$.

### 2.3. Polytopes
For a finite set of points $P = \\{x\_1, x\_2, ..., x\_n\\} \subseteq X$, its convex hull $\conv{P}$ is a convex polytope, whose dimension is defined to be that of the minimal plane (i.e. great spheres in case $X = \mathbb{S}^d$) containing $P$. Conversely, any convex polytope can be represented as the convex hull of such a set $P$, giving the so-called $\mathcal{V}$-representation. By abuse of notation, we shall refer to $P$ as the polytope $\conv{P}$ whenever the context is apparent. It is clear that $\conv{P}$ is compact, and in what follows, we only consider convex polytopes, thus omit the word "convex" when speaking of polytopes. 

A face $F$ of $P$ is said to be $k$-dimensional if it is contained in a minimal $k$-plane. The faces of $P$ are themselves polytopes, with dimension ranging from 0, the _vertices_, to 1, the _edges_, and up to $d$, the _facets_, and $d+1$, which is $P$. Two faces are said to be adjacent if their intersection is non-empty. We denote $\vertex{P}$ to be the set of vertices of $P$.

In case $X = \mathbb{R}^{d+1}$, $P$ admits a normal fan $\mathcal{N}(P)$ associating to each face $F$ of $P$ the set of vector $y$ such that $h_P(y)$ is attained by only and any point $x \in F$, which can be shown to be a cone. For each cone $N \in \mathcal{N}(P)$, we define the normal spherical polytope $S = \mathbb{S}^d \cap N$. By convexity of $P$, $\mathcal{N}(P)$ is a complete fan, so the collection of such polytopes, denoted $\mathcal{S}(P)$, covers $\mathbb{S}^d$.

By abuse of notation, we denote $h\_{\conv{P}}$ and $\sigma\_{\conv{P}}$ by $h_P$ and $\sigma_P$ whenever they are well-defined, respectively. Let $A$ and $B$ be two polytopes, then one can show that $A + B$ is also a polytope.

By Minkowski-Weyl theorem, a polytope can be represented either as a set of vertices, also known as $\mathcal{V}$-representation, or as a bounded intersection of some halfspaces, which in $\mathbb{R}^{d+1}$ amounts to specifying a normal vector for each facet, also known as $\mathcal{H}$-representation. Thus, we shall call the complexity of a polytope $P$ to be the total number of its vertices and facets.

## 3. RNA Polytopes

### 3.1. Original motivation

The parametric analysis problem originally arose from sequence alignment algorithms where it was unclear how to specify the parameters. Thus, one wishes to study what parameter sets will give rise to what alignment, i.e. a way to classify parameter sets. In one of the earliest attempts, Fitch and Smith {% cite Fitch1983 -f rna_polytopes %} considered a 2-parameter function to measure similarity in Needleman-Wunsch algorithm, and studied two short sequences derived from mRNA of chicken $\alpha$- and $\beta$-hemoglobin, where they identified 11 possible optimal solutions. Their method involved computing alignment for a number of parameter sets to identify a region in which any parameter set would yield the same alignment. They identified the region and moved to its neighbours, sequentially searched through the parameter space until no such region could be found. Nonetheless, this approach requires heavy computation, redundant alignments, and overall poses difficult as to argue about the regions, e.g. proving that no other regions can exist.

One crucial observation is by the discrete nature of alignments (or, in case of RNA, secondary structures), there exist necessarily finitely many possible optimal solutions, even if one considers _all_ uncountably many possible parameter sets. This number is further reduced since we consider not the alignments (resp. secondary structures) themselves, but some features thereof.

In our example, we focus only on possible values of $c_{AU}$, $c_{GC}$, and $c_{GU}$. A brute-force approach reveals that for $q = UAUUCUGAUG$, there are 67 compatible secondary structures, whence arise 15 possible signatures, thus there can be no more than 15 possible optimal structures even if we consider all parameter sets. These signatures are given below.
<div align="center">
	$(0, 0, 0)$, $(0, 0, 1)$, $(0, 1, 0)$, $(0, 1, 1)$, $(0, 0, 2)$, <br>
	$(1, 0, 0)$, $(1, 0, 1)$, $(1, 1, 0)$, $(1, 1, 1)$, $(1, 0, 2)$, <br>
	$(2, 0, 0)$, $(2, 0, 1)$, $(2, 1, 0)$, $(2, 1, 1)$, $(2, 0, 2)$
</div>

Unfortunately, this approach is not generalisable as the length of $q$ grows: intuitively, a longer sequence is expected to have more possible structures, and if nucleotides are uniformly distributed, then this number can grow exponentially. Indeed, Zuker and Sankoff {% cite ZUKER1984 -f rna_polytopes %} demonstrated the following theorem.

**Theorem 1.** {% cite ZUKER1984 -f rna_polytopes %} Let $q$ be a sequence of length $n$, whose bases are given by i.i.d random variables with the probability of occurrence for $A$, $G$, $C$, $U$ to be $a$, $g$, $c$, and $u$ respectively. Denote $p = 2(au + gc)$, $\alpha = \left(\frac{1 + \sqrt{1+4\sqrt{p}}}{2}\right)^2$, and $H = \frac{\alpha(1+4\sqrt{p})^{\frac{1}{4}}}{2\sqrt{\pi}p^\frac{3}{4}}$, then as $n$ tends to infinity, one has the expected number of compatible secondary structures $E(n)$ to be
\\[
    E(n) \sim Hn^{-\frac{3}{2}}\alpha^n.
\\]
In particular, $p = \frac{1}{4}$ for the case where all nucleotides can occur with equal probability gives $\alpha = \frac{1+\sqrt{3}}{2} = 1.866...$.

Much less is known about the number of possible signatures, as it is necessarily model-dependent. For our model counting the number of base pairs for each type and in general for any energy model, it is expected that the small number of features will greatly reduces the possible signatures, for the value of each feature is bounded by the length $\|q\|$ of $q$. Therefore, assuming the number of features is fixed, the number of possible signatures will be bounded by a polynomial of $\|q\|$. But we remind that it is one thing to compute _the number of signatures_, it is another thing to compute _the signatures_ themselves: even in our model where signatures necessarily have integer coordinates, a priori there is no viable way to know which signatures our model can predict.

### 3.2. Reduction to polytopes

The next crucial observation is we need not care about all signatures, but only those whom our model can predict. Mathematically, suppose a parameter $p$, a model $E$ given by a set of features given by a function $c(\cdot, \cdot)$, and a sequence $q$, a dynamic programming algorithm will then compute
\\[
	\max_{s'} \langle c(q, s'), p \rangle,
\\]
where $s'$ ranges over all compatible secondary structures. This motivates us to define the RNA polytope of $q$ (in the model $E$) as
\\[
	\mathcal{P}(q) = \conv \\{c(q, s) \mid \text{secondary structure } s\\},
\\]
then a structure $s$ that our model can predict necessarily satisfies $\langle c(q, s), p \rangle = h_{\mathcal{P}(q)}(p)$, or equivalently $c(q, s) \in \partial \mathcal{P}(q)$. We call such a pair $(q, s)$ _learnable_ (for model $E$), and we may restrict ourselves to studying only the boundary of $P$.

Back to our example, from the list of 15 possible signatures, we have the polytope $\mathcal{P}(q)$ and its spherical normal polytopes $\mathcal{S}(\mathcal{P}(q))$ shown below.
![](/littlenotebook/assets/img/rna_polytopes/rna_polytope.png#center)
<div align="center">Figure 1a. Example of an RNA polytope, for $q = UAUUCUGAUG$ ...</div>
![](/littlenotebook/assets/img/rna_polytopes/rna_normal_polytope.png#center)
<div align="center">Figure 1b. ... and its spherical normal polytope.</div>
We have three remarks:
- Given a parameter set $p$, finding the signature it will predict corresponds to finding the spherical polytope it belongs to, thus this gives us a complete classification of parameter sets.
- In our example, it happens to be the case that all signatures lie on the boundary of $\mathcal{P}(q)$, but as the length increases, this phenomenon is not to be expected in general. Whilst we found no results concerning the complexity of $\mathcal{P}(q)$ as it depends not only on the number of possible signatures, but also _their distribution_ {% cite Har-Peled2011 -f rna_polytopes %}, which has not been well-studied for any model. Regarding the sequent alignment problem, amongst $O(n^d)$ signatures, where $n$ and $d$ denote the total length of two sequences and the dimension, Gusfield et al. {% cite Gusfield1994 -f rna_polytopes %} showed that only $O\left(n^{d\frac{d-1}{d+1}}\right)$ points lie on the boundary for $d = 2$, and Pachter and Sturmfels {% cite Pachter2004 -f rna_polytopes %} showed the same for general $d$.
- Finally, this phenomenon of unlearnability is not restricted to high dimension, and in fact it is usually the opposite: for a given sequence, lower dimensions allow fewer signatures to be learnable. One of my other posts [(Nguyen, 2024)](relative_position_problem_motivation.html) in particular shows how it can happen for $d = 2$.

In the early 1990s, various methods were proposed to construct systematically the decomposition of parameter space with less computations, such as that by Fernandez-Baca and Srinivasan {% cite Fernandez-Baca1991 -f rna_polytopes %} (amongst others {% cite Vingron1994 -f rna_polytopes %}). Their algorithm involves finding an initial convex set of points in the polytope $P$, then gradually extending the set until it matches the boundary of $P$. This scheme is applicable to general dimension, with complexity to be polynomial of that of $P$ and $\sigma_P$.

Nonetheless, the methods until then treated $\sigma\_P$ and $h_P$ as black boxes, which both introduced unnecessary overheat and was not entirely satisfied in the context of bioinformatics in general, where many of the algorithms are dynamic-programming-based. Examples include Needleman-Wunsch and Smith-Waterson algorithm for sequence alignment, Nussinov's and Zuker's algorithm for RNA secondary structure, and Fitch's algorithm for phylogenetic tree construction. Since the optimal solution is constructed by solving subproblems, one may ask if it is possible to construct $\mathcal{P}(p)$ in a similar fashion by modifying the original dynamic programming scheme. One then may hope to find a more efficient mean to construct the polytopes, and even gather information about it without explicit constructions.

### 3.3. Dynamic programming in polytope algebra

In the context of graphical models, Pachter and Sturmfels {% cite Pachter2004 Pachter2005 -f rna_polytopes %} demonstrated that an alternative approach to construct the polytope from the dynamic programming scheme is possible. In particular, they showed that for any models whose optimal value is a linear combination of the parameters which can be found by a max-sum decomposition, there is a natural extension of the dynamic programming to compute the corresponding polytope of the model. In particular, recall the (max) tropical semi-ring given by $(\mathbb{R}\cup\{-\infty\}, \max(\cdot, \cdot), +)$, one can consider the same max-sum decomposition but in polytope algebra, given by $(\{\text{convex polytopes}\}, \oplus, \otimes)$, where for any two given polytopes $A$ and $B$, one defines

\\[
	A \oplus B = \conv\\{A \cup B\\}, A \otimes B = A + B.
\\]

![](/littlenotebook/assets/img/rna_polytopes/A.png#center)
<div align="center">Figure 2a. $A$.</div>
![](/littlenotebook/assets/img/rna_polytopes/B.png#center)
<div align="center">Figure 2b. $B$</div>
![](/littlenotebook/assets/img/rna_polytopes/AoplusB.png#center)
<div align="center">Figure 2c. $A \oplus B$</div>
![](/littlenotebook/assets/img/rna_polytopes/AotimesB.png#center)
<div align="center">Figure 2d. $A \otimes B$</div>

A more thorough treatment can be found in original articles, but for sake of simplicity, in this thesis, we shall derive the scheme in a more intuitive fashion, and without machinery used by Pachter and Sturmfels, such as Newton polytopes.

First, we note that the Nussinov's algorithm in our example uses only two operations, namely addition and taking the max. It turns out to be the paradigm for many other dynamic programming algorithm, that Tendeau formalised in definition {% cite Tendeau1998 -f rna_polytopes %}, and thus in this case, instead of computing in full the dynamic programming to trace out the boundary of $\mathcal{P}(q)$, we only need to keep track of how the polytopes evolve as the algorithm is executed.

Given two terms $a(p)$ and $b(p)$ (which stands for values of form $f(i, j)$ for some $i$ and $j$ in our cases) whose exact values depend on the choice of parameter set $p$ at the beginning of execution, to which we associate two polytopes $A$ and $B$. Recalling the definition of such polytopes, we have that $a(p) = h\_A(p)$, and likewise, $b(p) = h\_B(p)$. Therefore, suppose we associate to the term $\max(a(p), b(p))$ a polytope $C$, then one must have 
\\[
    h\_C(p) = \max(a(p), b(p)) = \max(h\_A(p), h\_B(p))
\\]
for all $p$; similarly, to the term $a(p) + b(p)$ the corresponding polytope $D$ must satisfy 
\\[
    h\_D(p) = a(p) + b(p) = h_A(p) + h\_B(p).
\\]

On the other hand, consider the polytope $C' = A \oplus B$ and $D' = A \otimes B$, and some $p \in \mathbb{R}^{d+1}$, then by construction, one has
\\[
	h_{C'}(p) = \max_{x \in C'} \langle x, p \rangle = \max_{x \in A \cup B} \langle x, p \rangle = \max\left(\max_{x \in A} \langle x, p \rangle, \max_{x \in B} \langle x, p \rangle \right) = \max(h_A(p), h_B(p)),
\\]
and similarly, 
\\[
	h_{D'}(p) = \max_{x \in D'} \langle x, p \rangle = \max_{(a, b) \in A \times B} \langle a + b, p \rangle = \max_{(a, b) \in A \times B} \langle a, p \rangle + \langle b, p \rangle = \max_{a \in A} \langle x, p \rangle + \max_{b \in B} \langle b, p \rangle,
\\]
so we conclude that $h\_{C'} = h\_C$ and $h\_{D'} = h\_D$. As a convex set is determined by its supporting function, we conclude that $C = A \oplus B$ and $D = A \otimes B$.

With a mathematical foundation laid down, we can at least in theory carry out the modified dynamic programming to study the parametric space. For instance, going back to our example, we can modify the Nussinov's algorithm to polytope algebra, yielding the following recurrent relation, where $P\_{AU} = (1, 0, 0)^T$, $P_{GC} = (0, 1, 0)^T$, $P_{GU} = (0, 0, 1)^T$, and $F(i, j)$ denotes the polytope associated with the $q\_i q\_{i+1} ... q\_{j-1} q\_j$.

<div align="center">
$$
F(i, j) = \bigoplus
\begin{cases}
	F(i+1, j) \\
	F(i, j-1) \\
	F(i-1, j+1) \otimes P_{\{q_i, q_j\}} \text{ if } q_i \text{ and } q_j \text{ can form a base pair} \\
	\bigoplus_{i < k < j} F(i, k) \otimes F(k+1, j)
\end{cases}.
$$
</div>

Similar to original Nussinov's algorithm mentioned above, this recurrent equation has redundancy, which can be overcome by alternative formulations. But this shall not change the final output, nor alter the complexity up to a constant factor.

### Goal

Using polytopes, one can in theory prove if a given parameter set is learnable (meaning that there exist some secondary structure optimal with respect to this parameter set) by deciding if the corresponding point lies on the boundary of the RNA polytope. And suppose that it is is the case, one can even decide how robust it is, i.e. suppose $s$ is an optimal secondary structure with respect to $q$, how much one must deviate from $q$ so that $s$ is no longer optimal. Both problems are studied in subsequent posts ([Nguyen, 2024](relative_position_problem_motivation.html); [Nguyen, 2024](RNA_robustness_problem_telescoping_method)).

### References

{% bibliography --cited --file rna_polytopes %}