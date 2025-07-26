---
layout: post
title: "An Ising model in RNA folding (part 2)"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\APX}{\textsf{APX}}
\newcommand{\sP}{\textsf{#P}}
$$

This is part 2/3 of my Bachelor thesis during my first year at Ulm, carried out at Hamilton Institute, Maynooth University and under supervision of Prof. Damien Woods, titled "Computational complexity of an Ising model in pseudoknotted nucleic acid folding", to be published.

## 5. Complexity of $$\textsf{SPIN}_\Delta$$

All reductions used to prove complexity hardness results regarding $$\textsf{SPIN}$$ and $$\textsf{MFE}$$ rely on the following two essential observations. In this section, it is always assumed without loss of generality that $$S \in \mathcal{S}$$.

Firstly, only the spins of paired nucleotides are taken into account in $$\Delta_1$$ or $$\Delta_2$$. This means we can remove all the unpaired nucleotides, and thus split the original strand $$q$$ into a number of contiguous subsequences $$q_1, q_2, \cdots, q_m$$, where all nucleotides within $$q_i$$ for a given $$i$$ are paired (not necessarily with each other). Denote by $$\{a_i, a_i + 1, \cdots, b_i\}$$ the positions of $$q_i$$'s nucleotides in $$q$$, we can rewrite $$\Delta_1$$ as follow

$$
	\Delta_1(S, \sigma) = \alpha \sum_{(i, j) \in S} \sigma_i \sigma_j - \sum_{k = 1}^m \sum_{i = a_k}^{b_k-1} \sigma_i \sigma_{i+1}.
$$


Secondly, and more importantly, as discussed in the derivation of Ising model,  $$\alpha$$ balances the strength of two interactions between being pairing and adjacency. Intuitively, if $$\alpha$$ is small enough, it is always more favourable that adjacent nucleotides have the same spin _regardless_ of how they are paired. This intuition is made precise in the following lemma.

**Lemma 4.** Let $$S$$ be any secondary structure, and $$\sigma^S = \min_{\sigma} \Delta_1 (S, \sigma)$$. If $$\alpha < \frac{2}{n}$$, then for all $$1 \leq k \leq m$$ and $$a_k \leq i, j \leq b_k$$, $$\sigma^S_i = \sigma^S_j$$. The same holds for $$\Delta_2$$.

_Proof._ Suppose the contrary, let $$1 \leq k \leq m$$ and $$a_k \leq \ell \leq b_k$$ such that $$\sigma^S_{\ell+1} \neq \sigma^S_\ell$$. Suppose that at least half the nucleotides between $$a_k$$ and $$b_k$$ have the spin $$\delta$$. Let $$\sigma_j = \begin{cases} \delta & \text{ if } a_k \leq j \leq b_k \\ \sigma^S_j & \text{ otherwise} \end{cases}$$, then

$$
\begin{align*}
    \Delta_1(S, \sigma) - \Delta_1(S, \sigma^S) 
    & = \alpha \sum_{a_k \leq i \leq b_k \mid (i, j) \in S} (\delta - \sigma^S_i) \sigma_j - \sum_{i = a_k}^{b_k - 1} (1 - \sigma^S_i \sigma^S_{i+1}),
\end{align*}
$$

where

$$
\begin{align*}
    \left|\alpha \sum_{a_k \leq i \leq b_k \mid (i, j) \in S} (\delta - \sigma^S_i) \sigma^S_j\right|
    & \leq \alpha \sum_{a_k \leq i \leq b_k \mid (i, j) \in S} |(\delta - \sigma^S_i) \sigma^S_j| \\
    & = \alpha \sum_{a_k \leq i \leq b_k \mid (i, j) \in S} |\delta - \sigma^S_i| \\
    & = 2\alpha | \{ i \mid a_k \leq i \leq b_k, \sigma^S_i = -\delta \} \\
    & \leq \alpha (b_k - a_k + 1) \leq \alpha n < 2,
\end{align*}
$$

and as $$\sigma^S_i \sigma^S_{i+1} \leq 1$$, we have $$\sum_{i = a_k}^{b_k - 1} (1 - \sigma^S_i \sigma^S_{i+1}) \geq (1 - \sigma^S_\ell \sigma^S_{\ell+1}) = 2$$, which implies $$\Delta_1(S, \sigma) - \Delta_1(S, \sigma^S) < 0$$, as desired.

Entirely analogously, we have

$$
    \Delta_2(S, \sigma) - \Delta_2(S, \sigma^S) 
    = \alpha \sum_{a_k \leq i \leq b_k \mid (i, j), (i+1, j-1) \in S} (\delta - \sigma^S_i) \sigma_j - \sum_{i = a_k}^{b_k - 1} (1 - \sigma^S_i \sigma^S_{i+1}),
$$

where

$$
\begin{align*}
    \left|\alpha \sum_{a_k \leq i \leq b_k \mid (i, j), (i+1, j-1) \in S} (\delta - \sigma^S_i) \sigma^S_j\right|
    & \leq \alpha \sum_{a_k \leq i \leq b_k \mid (i, j), (i+1, j-1) \in S} |(\delta - \sigma^S_i) \sigma^S_j| \\
    & = \alpha \sum_{a_k \leq i \leq b_k \mid (i, j), (i+1, j-1) \in S} |\delta - \sigma^S_i| \\
    & = 2\alpha | \{ i \mid a_k \leq i \leq b_k, \sigma^S_i = -\delta \} \\
    & \leq \alpha (b_k - a_k + 1) \leq \alpha n < 2,
\end{align*}
$$

and the rest follows similarly. <span style="float:right;">$$\square$$</span>

### 5.1. $$\sP$$-hardness of $$\textsf{SPIN}_{\Delta_1}$$

We reduce to $$\textsf{SPIN}_{\Delta_1}$$ from $$\textsf{MAXCUT}$$ problem, which is $$\sP$$-hard {% cite Jerrum1993 -f rna_ising_model %}. Let $$G = (V, E)$$ be a graph, a cut is a subset $$U \subseteq V$$. Denote $$\overline{U} = V \setminus U$$ and by $$E(U, \overline{U})$$ the set of edges with one endpoint in $$U$$ and the other in $$\overline{U}$$.

**Problem $$\textsf{MAXCUT}$$.** Given a graph $$G = (V, E)$$ as input, output $\max_U \| E(U, \overline{U}) \|$.

Given a graph $$G = (V, E)$$, we number the vertices $$V = \{v_1, v_2, \cdots, v_{\|V\|}\}$$. Denote $$d_i = \deg v_i$$ the degree of $$v_i$$, we number $$v_i$$'s neighbours as $$v_{i, 1}, v_{i, 2}, \cdots, v_{i, d_i}$$. To $$v_{i, j}$$, we assign the index $$I_{i, j} = i-1 + \sum_{k = 1}^{i-1} d_i + j-1$$.

Consider a nucleic acid sequence of length $n = 2\|E\| + \|V\| - 1 = \|V\| - 1 + \sum_{i = 1}^{\|V\|} d_i$, we construct the following secondary structure $$S$$: for an edge $$(v_i, v_j) \in E$$ and suppose that $$v_i = v_{j, \ell}$$ and $$v_j = v_{i, k}$$, we have $$(I_{i, k}, I_{j, \ell}) \in S$$. It is clear that $\|S\| = \|E\|$, and each cut $$U$$ defines a spin configuration $$\sigma$$, say, by assigning $$\sigma_{I_{i, j}} = 1$$ if and only if $$v_i \in U$$.

Now here is the intuition: since each nucleotide can bond with another one other, to have arbitrary degree for vertices in $$G$$, we encode the neighbourhood of each vertex as a contiguous subsequence in the nucleic acid. These subsequences are separated by one nucleotide and all nucleotides within one are paired, thus by choosing $$\alpha < \frac{2}{n}$$ and Lemma 4, in an optimal spin configuration each subsequence has all its nucleotides having the same spin. This partitions $$V$$ into $$U$$ and $$\overline{U}$$.

More rigorously speaking, by construction, we have that the contiguous subsequences are $$P_i = \{I_{i, 1}, I_{i, 2}, \cdots, I_{i, d_i}\}$$. Let $$\sigma^S$$ be an optimal spin configuration, i.e. $$\Delta_1(S, \sigma^S) = \min_{\sigma} \Delta_1 (S, \sigma)$$. By Lemma 4, we have that for all $1 \leq k \leq \|V\|$ and $$1 \leq i, j \leq d_k$$, we have that $$\sigma^S_i = \sigma^S_j = \delta_i$$, thus we can unambiguously define the spin $$\delta_i$$ as the spin of nucleotides in $$P_i$$. Let $$U^S = \{v_i \mid \delta_i = 1\}$$, then $$\overline{U^S} = \{v_i \mid \delta_i = -1\}$$. Therefore, 

$$
    \sum_{(i, j) \in S} \sigma_i \sigma_j = |E \setminus E(U^S, \overline{U^S})| - |E(U^S, \overline{U^S})| = |E| - 2|E(U^S, \overline{U^S})|.
$$


On the other hand, 

$$\sum_{k = 1}^{|V|} \sum_{i = I_{k, 1}}^{I_{k, d_k}-1} \sigma_i \sigma_{i+1} = \sum_{k = 1}^{|V|} \sum_{i = I_{k, 1}}^{I_{k, d_k}-1} 1 = 2|E| - |V|,$$

thus 

$$
\begin{align*}
	\Delta_1(S, \sigma^S) 
	& = \alpha \sum_{(i, j) \in S} \sigma_i \sigma_j - \sum_{k = 1}^{|V|} \sum_{i = I_{k, 1}}^{I_{k, d_k}-1} \sigma_i \sigma_{i+1} \\
	& = \alpha (|E| - 2|E(U^S, \overline{U^S})|) - 2|E| + |V|, 
\end{align*}
$$

and as $$\Delta_1(S, \sigma^S) = \min_{\sigma} \Delta_1(S, \sigma)$$, we have that $\|E(U^S, \overline{U^S})\| = \max_U \|E(U, \overline{U})\|$.

Thus to conclude, not only solving $$\textsf{SPIN}_{\Delta_1}$$ gives a solution to $$\textsf{MAXCUT}$$, there is a one-to-one correspondence between maximum cuts of $$G$$ and optimal spin configuration of $$S$$. Indeed, each optimal spin configuration corresponds to one unique cut, and vice versa.

### 5.2. $$\APX$$-hardness of $$\textsf{SPIN}_{\Delta_1}$$

As for $$\APX$$-hardness, whilst it appears to be at least complicated, if not impossible, to give an approximation-preserving reduction from the above construction, we rely on the fact that $$\textsf{MAXCUT}$$ remains $$\APX$$-hard for cubic graphs {% cite Alimonti2000 -f rna_ising_model %}. This gives a strengthened version of Lemma 4. In this subsection, $$G = (V, E)$$ is a simple cubic graph, i.e. $$d_v = 3$$ for all $$v \in V$$.

**Lemma 5.** Let $$S$$ be any secondary structure such that paired nucleotides form contiguous subsequences of length at most $$3$$, and $$\sigma^S = \min_\sigma \Delta_1(S, \sigma)$$. If $$\alpha < 1$$, then for all $$1 \leq k \leq m$$ and $$a_k \leq i, j \leq b_k$$, $$\sigma^S_i = \sigma^S_j$$.

_Proof._ Consider three consecutive paired nucleotides at position $$i-1$$, $$i$$, and $$i+1$$ for some $$i$$. If $$\sigma^S_{i-1} = \sigma^S_i = \sigma^S_{i+1}$$ then there is nothing to be done. Otherwise, exactly two nucleotides amongst them have the same spin. 
	
Consider $$\Delta_1$$, changing the spin of the remaining nucleotide reduces the energy by at least $$2$$ via its contribution to the term $$-\sum_{j, j+1 \in P} \sigma_j \sigma_{j+1}$$, and increases the energy by at most $$2 \alpha$$ via its contribution to the term $$\alpha \sum_{(j, k) \in S} \sigma_j \sigma_k$$. If $$\alpha < 1$$, then this change is more favourable. The same argument applies when there are only two adjacent paired nucleotides. <span style="float:right;">$$\square$$</span>

**Remark 6.** Upon a closer inspection of the proof of Lemma 4, one sees that it suffices to have $$\alpha < \frac{2}{\max_k (b_k - a_k + 1)}$$, which means we can choose $$\alpha < \frac{2}{3}$$ for example. Lemma 5 provides a strengthened version, and captures an intuition that is not evident in Lemma 4. Perhaps it will be beneficial to make this intuition precise for future works.

Thus in view of this lemma or the succeeding remark, we can take $$\alpha = \frac{1}{2}$$ in $$\Delta_1$$ for the construction above corresponding to $$G$$.

In this construction, each vertex $$u$$ corresponds to exactly three nucleotides at position $$I_{u, 1}$$, $$I_{u, 2}$$, and $$I_{u, 3}$$. Suppose we have an algorithm that returns a spin configuration $$\sigma$$ giving a $$\gamma$$-approximation to $$\textsf{SPIN}_{\Delta_1}$$. In view of Lemma 5, it is always better that $$\sigma_{I_{u, 1}} = \sigma_{I_{u, 2}} = \sigma_{I_{u, 3}}$$, and amongst three values, we always have at least two equal, thus we can change the spin of at most one nucleotide to obtain a new spin configuration giving the energy no worse than that of the original. Thus we may assume that $$\sigma_{I_{u, 1}} = \sigma_{I_{u, 2}} = \sigma_{I_{u, 3}}$$ for all $$u \in V$$. This gives rise to a cut $$U_{\sigma}$$.

Denoting $$C = \max_{U'} E(U', \overline{U'})$$, noting that $2\|E\| = 3\|V\|$ as $$G$$ is cubic, calculation as above gives

$$
\begin{align*}
	\min_{\sigma'} \Delta_1(S, \sigma') 
	& = \frac{1}{2} (|E| - 2C) - 2|E| + |V| \\
	& = \frac{1}{2} \left(\frac{3}{2}|V| - 2C\right) - 6|V| + |V| \\ 
	& = -C - \frac{17}{4}|V|, 
\end{align*}
$$

and

$$
	\Delta_1(S, \sigma) = -|E(U_{\sigma}, \overline{U_{\sigma}})| - \frac{17}{4}|V|.
$$


Now, by an elementary argument of a probabilistic method, note that we have a bound $C \geq \frac{\|E\|}{2} = \frac{3}{4} \|V\|$. For all $$0 \leq r \leq 1$$ and  $$\gamma \geq 1 - \frac{3}{20}r$$, we have

$$
    \geq 1 - \frac{3}{20}r \Leftrightarrow 20 \gamma \geq 20 - 3r \Leftrightarrow 3(\gamma - 1 + r) \geq 17(1 - \gamma)
$$


$$
    \Rightarrow (\gamma - 1 + r) C \geq \frac{3(\gamma - 1 + r)}{4} |V| \geq \frac{17}{4} (1 - \gamma) |V|,
$$

which, together with $$\Delta_1(S, \sigma) \leq \gamma \min_{\sigma'} \Delta_1(S, \sigma')$$, implies

$$
\begin{align*}
	|E(U_{\sigma}, \overline{U_{\sigma}})| + \frac{17}{4}|V|
	& \geq \gamma C + \gamma \frac{17}{4}|V|\\
	& = (1 - r) C + (\gamma - 1 + r) C + \gamma \frac{17}{4}|V| \\
	& \geq (1 - r) C + \frac{17}{4}|V| \\
	\Rightarrow |E(U_{\sigma}, \overline{U_{\sigma}})| & \geq (1 - r) C
\end{align*}
$$

Thus we have an $$\textsf{L}$$-reduction and $$\textsf{SPIN}_{\Delta_1}$$ is $$\APX$$-hard.

### 5.3. Complexity of $$\textsf{SPIN}_{\Delta_2}$$

The reduction goes more or less the same, except now instead of one bond for each edge of $$G$$, we have a stack, i.e. two bonds of form $$(i, j)$$ and $$(i+1, j-1)$$.

More formally, let $$G = (V, E)$$ be a simple graph, we number the vertices $$V = \{v_1, v_2, \cdots, v_{\|V\|}\}$$. Denote $$d_i = \deg v_i$$ the degree of $$v_i$$, we number the neighbours of $$v_i$$ as $$v_{i, 1}, v_{i, 2}, \cdots, v_{i, d_i}$$. To $$v_{i, j}$$, we assign two indices $$I_{i, j} = i-1 + 2\sum_{k = 1}^{i-1} d_i + 2(j-1)$$ and $$J_{i, j} = I_{i, j}+1$$.

Consider a nucleic acid sequence of length $$n = 4\|E\| + \|V\| - 1$$, we construct the following secondary structure $$S$$: for an edge $$(v_i, v_j) \in E$$ with $$i < j$$ and suppose that $$v_i = v_{j, \ell}$$ and $$v_j = v_{i, k}$$, we have $$(I_{i, k}, J_{j, \ell}), (J_{i, k}, I_{j, \ell})\in S$$. It is clear that $$\|S\| = \|E\|$$, and each cut $$U$$ defines a spin configuration $$\sigma$$, say, by assigning $$\sigma_{I_{i, j}} = 1$$ if and only if $$v_i \in U$$.

Now each edge $$(v_i, v_j) \in E$$ corresponds to the stack formed by $$(I_{i, k}, J_{j, \ell})$$ and $$(J_{i, k}, I_{j, \ell})$$ and vice versa. Indeed, consider a stack formed by $$(i, j)$$ and $$(i+1, j-1)$$, then since the contiguous subsequences of vertices are separated by one unpaired nucleotide, $$i$$ and $$i+1$$ must belong to a subsequence of some vertex, and similarly do $$j$$ and $$j-1$$. As $$G$$ is simple, this stack corresponds to some edge. The proof of $$\sP$$- and $$\APX$$-hardness go analogously.

### References

{% bibliography --cited -f rna_ising_model %}