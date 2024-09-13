---
layout: post
title: "An Ising model in RNA folding: Complexity of $\\textsf{MFE}$"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\NP}{\textsf{NP}}
\newcommand{\APX}{\textsf{APX}}
\newcommand{\sP}{\textsf{#P}}
\newcommand{\A}{\textsf{A}}
\newcommand{\T}{\textsf{T}}
\newcommand{\G}{\textsf{G}}
\newcommand{\C}{\textsf{C}}
$$

## 6. Complexity of $\textsf{MFE}_\Delta$
### 6.1. $\sP$-hardness of $\textsf{MFE}\_\{\Delta\_1\}$ and of $\textsf{MFE}\_\{\Delta\_2\}$

The idea is similar as in the previous section, except now it is hard to construct a sequence that admits a specific secondary structure as its optimal. To get around this constraint, we construct a specific alphabet $\Sigma$.

Let $G = (V, E)$ be an undirected graph, we first number the edges of $G$, i.e. $E = \\{(u_i, v_i) \mid 1 \leq i \leq \|E\|\\}$. Then, consider the alphabet $\Sigma = \\{s\_i, t\_i \mid 1 \leq i \leq \|E\|\\} \cup \\{\\#\\}$, where the only allowed pairings are between $s_i$ and $t_i$ for some $i$, and the base $\\#$ pairs with no bases.

Let $v$ be a vertex, and suppose $i_1, i_2, \cdots, i_{d_v}$ are the indices of edges having $v$ as one of the endpoints. For $1 \leq j \leq d_v$, if $v = u_{i_j}$, then we take $s_{i_j}$, else we take $t_{i_j}$. Concatenating all these nucleotides, we obtain a sequence of nucleotides $q_v$. Finally, we consider the sequence $q = \overline\{q_1 \\# q_2 \\# \cdots \\# q_\{\|V\|\}\}$ and $\alpha < \frac{2}{n}$ where $n = \|q\|$. By Lemma 4, we have that given any secondary structure, in an optimal spin configuration, for all $1 \leq i \leq \|V\|$, all nucleotides in $q_i$ have the same spin. This partitions $V$ into two subsets $U$ and $\overline{U} = V \setminus U$, comprising vertices whose subsequence is assigned spin $1$ and  $-1$, respectively. 

Now consider an optimal secondary structure $S^\*$ along with an optimal spin configuration $\sigma^\*$, i.e. $\Delta_1(S^\*, \sigma^\*) = \min_{S, \sigma} \Delta_1(S, \sigma)$, which gives rise to a cut $U^\*$. $S^\*$ will not include pairs $(i, j)$ having $\sigma^\*_i = \sigma^\*_j$, since it is more beneficial not to include such pairs. Thus, every $(i, j) \in S^\*$ has $\sigma^\*_i = -\sigma^\*_j$, meaning they correspond to an edge between two vertices, one in $U^\*$ and another in $\overline{U^\*}$. Conversely, any such edge also corresponds to two nucleotides which can be paired, and favourably so since they have opposite spin.

Computation as in previous section gives
\\[
	\sum_{(i, j) \in S^\*} \sigma^\*_i \sigma^\*_j = -\|S^\*\| = - \|E(U^\*, \overline{U^\*})\|,
\\]
and

$$
\begin{align*}
	\Delta_1(S^*, \sigma^*) 
	& = \alpha \sum_{(i, j) \in S^*} \sigma^*_i \sigma^*_j - \sum_{k = 1}^{|V|} \sum_{i = I_{k, 1}}^{I_{k, d_k}-1} \sigma^*_i \sigma^*_{i+1} \\
	& = - \alpha |E(U^*, \overline{U^*})| - 2|E| + |V|, 
\end{align*}
$$

thus as $\Delta_1(S^\*, \sigma^\*) = \min_{S, \sigma} \Delta_1(S, \sigma)$, we have that $\|E(U^\*, \overline{U^\*})\| = \max_U \|E(U, \overline{U})\|$, i.e. $U$ is a maximum cut.

In fact, more can be said, since as in previous section, there is a one-to-one correspondences between maximum cuts and optimal spin configurations $\sigma^\*$ (meaning that there exists a secondary structure $S^\*$ such that $\Delta_1(S^\*, \sigma^\*) = \min_{S, \sigma} \Delta_1(S, \sigma)$). But, in our case, each possible pair $(s_i, t_i)$ appear _exactly once_, thus given an optimal spin configurations $\sigma^\*$, the secondary structure $S^\*$ such that $\Delta_1(S^\*, \sigma^\*) = \min_{S, \sigma} \Delta_1(S, \sigma)$, which must exist by definition, is unique. This means the reduction is parsimonious, and we have $\sP$-hardness of $\textsf{MFE}_\{\Delta_1\}$.

The $\sP$-hardness proof of $\textsf{MFE}\_{\Delta\_2}$ goes similarly as above, where the base pairs are now replaced with the stacks analogously to how we modified the $\sP$-hardness proof of $\textsf{SPIN}\_{\Delta\_1}$ to obtain that of $\textsf{SPIN}\_{\Delta\_2}$.

### 6.2. $\APX$-hardness of $\textsf{MFE}_{\Delta}$

Here we argue only for $\Delta_1$, with some modification as argued with $\APX$-hardness of $\textsf{SPIN}_\{\Delta_2\}$ to obtain one for $\textsf{MFE}\_\{\Delta\_2\}$.

The idea is to use Lemma 5 and argue as done for $\textsf{SPIN}_\{\Delta_1\}$. Consider a cubic graph $G = (V, E)$ and the construction of nucleic acid sequence as in the previous subsection. In view of Lemma 5, we can choose $\alpha = \frac{1}{2}$.

Suppose we have an algorithm that return a secondary structure $S$ and a spin configuration $\sigma$ giving a $\gamma$-approximation to $\textsf{MFE}\_\{\Delta\_1\}$. In view of Lemma 5, it is always better that $\sigma_{I_{u, 1}} = \sigma_{I_{u, 2}} = \sigma_{I_{u, 3}}$, and amongst three values, we always have at least two equal, thus we can change the spin of at most one nucleotide to obtain a new spin configuration giving the energy no worse than that of the original. Thus we may assume that $\sigma_{I_{u, 1}} = \sigma_{I_{u, 2}} = \sigma_{I_{u, 3}}$ for all $u \in V$. This gives rise to a cut $U_{\sigma}$.

Denoting $C = \max_{U'} E(U', \overline{U'})$, noting that $2\|E\| = 3\|V\|$ as $G$ is cubic, calculation as above gives $\min_{\sigma'} \Delta_1(S, \sigma') =  -C - 2\|V\|$ and $\Delta_1(S, \sigma) = -\|E(U_{\sigma}, \overline{U_{\sigma}})\| - 2\|V\|$. Entirely analogously as in the $\APX$-hardness proof of $\textsf{SPIN}\_\{\Delta\_1\}$, calculation gives that if $\gamma \geq 1 - \frac{3}{11}r$ then $\|E(U\_{\sigma}, \overline\{U\_{\sigma}\})\| \geq (1 - r) C$, as desired.

### 6.3. $\NP$-hardness of $\textsf{MFE}_\{\Delta_2\}$ with bounded $\Sigma$

In this subsection, we consider alphabet of Watson-Crick pairing $\Sigma = \\{\A, \T, \G, \C\\}$. Our reduction is almost the same to that of Lyngsø for BPS model. Recall $\textsf{BINPACKING}$ problem.

**Problem $\textsf{BINPACKING}$.** Given positive integers $B$ and $C$, and $k$ positive integers $a_1, \cdots, a_k$ as input, output $\textsf{YES}$ if there exists a partition of $\\{1, 2, \cdots, k\\}$ into $B$ sets $I_1, I_2, \cdots, I_B$ such that for all $1 \leq i \leq B$, one has $\sum_{j \in I_i} a_j \leq C$. Otherwise, output $\textsf{NO}$.

Given an instance of $\textsf{BINPACKING}$ problem, denote $A = \sum_{i = 1}^{k} a_i$, and consider the following sequence
\\[
    q = \C^{a_1} \A \C^{a_2} \A \cdots \A \C^{a_k} \left(\A \G^C\right)^B.
\\]

Consider a secondary structure, it is a graph formed by $\C$'s and $\G$'s along with the bonds between $\C$'s and $\G$'s. This graph is clearly bipartite and the only possible contiguous subsequences of paired nucleotides are those of $\C$'s or $\G$'s, thus an optimal spin configuration is, e.g. assigning $1$ to $\C$'s and $-1$ to $\G$'s.

Moreover, if $\sum_{i = 1}^{k} a_i > B \cdot C$ then the answer is obviously \textsf{NO}, thus we can assume  $\sum_{i = 1}^{k} a_i \leq B \cdot C$, which means that an optimal secondary structure will have all $\C$'s paired. Similarly, given a subsequence of form $\G^C$, in order to maximise its contribution to the term $\sum_{i, i+1 \in P} \sigma_i \sigma_{i+1}$, it must be that all paired nucleotides form a contiguous subsequence of $\G^C$.

Then, 

$$
\begin{align*}
	\Delta_2(S, \sigma) 
	& = \alpha \sum_{(i, j) \in S \mid (i+1, j-1) \in S} \sigma_i \sigma_j - \sum_{i, i+1 \in P} \sigma_i \sigma_{i+1} \\
	& = \alpha \Delta'_2 (S) - A + k,
\end{align*}
$$

thus an secondary structure $S^\*$ and a spin configuration $\sigma^\*$ such that $\Delta_2(S^\*, \sigma^\*) = \min_{S, \sigma} \Delta_2(S, \sigma)$ will have $\Delta'\_2(S^\*) = \min\_\{S\} \Delta'\_2(S)$. Arguing as in the reduction of Lyngsø {% cite Lyngsø2004 -f rna_ising_model %}, we have that $\min_{S, \sigma} \Delta_2(S, \sigma) = (1 + \alpha) (k - A)$ if and only if the answer to the given $\textsf{BINPACKING}$ problem is $\textsf{YES}$.

By strong $\NP$-hardness {% cite Garey1979 -f rna_ising_model %}, we can assume that $B$, $K$ and $a_i$'s are bounded by a polynomial of the input size, then the length of $q$ is $\|q\| = A + k - 1 + K(B+1)$ also bounded by the input size. This completes our reduction.

### 6.4. $\textsf{MFE}_\{\Delta_1\}$ is in $\|\Sigma\|-W[1]$

In the previous reduction, the size of alphabet $\Sigma$ grows with that of graph $G$, which prompts the question of whether a similar reduction can be made with an alphabet of bounded size. Here, we show that it is unlikely, as $\textsf{MFE}_\{\Delta_1\}$ admits a $W[1]$-algorithm with $\|\Sigma\|$ as the parameter.

Consider two types of nucleotides, i.e. two characters in the alphabet $\Sigma$, which can be paired, i.e. $s = s_i$ and $t = t_i$ for some $i$. In the sequence $q$, consider the sets of indices corresponding to occurrences of $s$ and $t$, denoted by $I_{s}$ and $I_{t}$ respectively. All the $s - t$ pairings are between an index in $I_{s}$ and another in $I_{t}$.

The key idea here is we can ``shuffle" the nucleotides in $I_s$ to minimise their contribution to the term $-\sum_{i, i+1 \in P} \sigma_i \sigma_{i+1}$ without affecting the term $\alpha \sum_{(i, j) \in S} \sigma_i \sigma_j$, since we know that they pair with $I_t$, and know precisely the locations of type-$t$ nucleotides. This reduces the number of possibilities we need to consider in order to find an optimal spin configuration.

To put this more rigorously, indices in $I_s$ form some number of contiguous subsequences in $q$, so we write $I_s = \bigcup_{i = 1}^A [l_i, r_i]$ for some $A$. Consider one of such subsequences, say $[l_k, r_k]$ for some $1 \leq k \leq A$. Its contribution to $\Delta_1(S, \sigma)$ are of two form: either to $\alpha \sum_{(i, j) \in S} \sigma_i \sigma_j$, or to $\sum_{i, i+1 \in P} \sigma_i \sigma_{i+1}$. If we permute the nucleotides in $[l_k, r_k]$, then the term $\alpha \sum_{(i, j) \in S} \sigma_i \sigma_j$ remains the same. To the term $\sum_{i, i+1 \in P} \sigma_i \sigma_{i+1}$, its contribution is precisely 

$$
	\sigma_{l_k - 1} \sigma_{l_k} \unicode{x1D7D9}_{l_k-1 \in P} + \sigma_{r_k} \sigma_{r_k + 1} \unicode{x1D7D9}_{r_k \in P} + \sum_{i = l_k}^{r_k-1} \sigma_i \sigma_{i+1}.
$$

Fixing $\sigma_{l_k}$ and $\sigma_{r_k}$, the first two terms do not change. On the other hand, 
\\[
	\sum_{i = l_k}^{r_k-1} \sigma_i \sigma_{i+1} = r_k - l_k + 1 - 2\|\\{l_k \leq i \leq r_k - 1 \mid \sigma_i \neq \sigma_{i-1}\\}\|,
\\]
thus maximising the left-hand size is minimising
\\[
    G_k = \|\\{l_k \leq i \leq r_k - 1 \mid \sigma_i \neq \sigma_{i-1}\\}\|.
\\]
We have a few cases
- If $\delta = \sigma_{l_k} \neq \sigma_{r_k}$, then $G_k \geq 1$, and $G_k = 1$ can be obtained by rearranging the nucleotides so that we can divide $[l_k, r_k]$ into two subsegments, one with all spins being $\delta$ and another with all spins being $-\delta$. More formally speaking, there exists $j$ such that for all $l_k \leq \ell \leq j$, $\sigma_\ell = \delta$, and for all $j+1 \leq \ell \leq r_k$, $\sigma_\ell = -\delta$.
- If $\delta = \sigma_{l_k} = \sigma_{r_k}$ and there is no index $l_k \leq j \leq r_k$ such that $\sigma_j \neq \delta$, then it is clear that $\sigma_j = \delta$ for all $l_k \leq j \leq r_k$.
- If $\delta = \sigma_{l_k} = \sigma_{r_k}$ and there is at least one $l_k \leq j \leq r_k$ such that $\sigma_j \neq \delta$, then $G_k \geq 2$, and $G_k = 2$ can be obtained similar to the first case.

Thus, without loss of generality, we can assume that there exists a critical index $l_k \leq j < r_k$ such that for all $l_k \leq \ell \leq j$, $\sigma_\ell = \sigma_{l_k}$, and for all $j < \ell < r_k$, $\sigma_\ell = -\sigma_{l_k}$. Note that this is true regardless of whichever nucleotide $s$ being considered. Furthermore at the end, we can assign arbitrary spins to unpaired nucleotides without affecting the energy, thus this does not affect the final result.

With this observation, we have a dynamic programming scheme. Denote by  
\\[
	\text{DP}(i, \delta\_l, \delta\_r, (u^{1}\_c)\_{c \in \Sigma}, (u^{-1}\_c)\_{c \in \Sigma})
\\]
the minimum free energy of $q$'s prefix formed by the first $i$ contiguous subsequences, each comprises nucleotides of one single types. Moreover, suppose the $i^{\text{th}}$ contiguous subsequence is the segment $[l, r]$, then we require that $\sigma_l = \delta_l$ and $\sigma_r = \delta_r$. For the character $c \in \Sigma$, we require that there remain $u^\{1\}\_c$ (resp. $u^\{-1\}\_c$) _unpaired_ nucleotides with spins assigned to $1$ (resp. $-1$). 

The transition is straightforward. Suppose the $i^{\text{th}}$ subsequence is  formed by nucleotide of type $s \in \Sigma$ which can only be paired with $t$. We choose a critical index $l_k \leq j < r_k$. The term $\sum_{i = l_k}^{r_k-1} \sigma_i \sigma_{i+1}$ is

$$
	A_{\delta_l, \delta_r, j} = 
	\begin{cases}
		r - l - 1 & \text{if } \delta_l \neq \delta_r \\
		r - l + 1 & \text{if } \delta_l = \delta_r \text{ and } j = r_k \\
		r - l - 3 & \text{otherwise} \\
	\end{cases}
$$

The number of newly created unpaired $s$'s with spin $\delta_l$ is $$\Delta u^{\delta_l}_s = j - l + 1 + \unicode{x1D7D9}_{\delta_l = \delta_r}$$ and similarly that of newly created unpaired $s$'s with spin $-\delta_l$ is $\Delta u^{-\delta_l}_s = r - l + 1 - u^{\delta_l}$. Then, we assign

$$
\begin{align*}
	\text{DP}(i, \delta_l, \delta_r, (v^{1}_c)_{c \in \Sigma}, (v^{-1}_c)_{c \in \Sigma}) \leftarrow \min 
	\left[\text{DP}(i, \delta_l, \delta_r, (v^{1}_c)_{c \in \Sigma}, (v^{-1}_c)_{c \in \Sigma}), \right. \\
	\left. \text{DP}(i-1, \varepsilon_l, \varepsilon_r, (u^{1}_c)_{c \in \Sigma}, (u^{-1}_c)_{c \in \Sigma}) - A_{\delta_l, \delta_r, j} - \varepsilon_r \cdot \delta_l \right.\\
	\left. - \alpha(\max(u^1_s + \Delta u^1_s, u^{-1}_t)) + \max(u^{-1}_s + \Delta u^{-1}_s, u^1_t))\right]
\end{align*}
$$

for $(v^{1}\_c)\_{c \in \Sigma}$ and $(v^{-1}\_c)\_{c \in \Sigma}$ defined by 

$$
	v^\delta_c = 
	\begin{cases}
		u^\delta_c & \text{if } c \neq s, t \\
		u^\delta_s + \Delta u^\delta_s - u^{-\delta}_t & \text{if } u^\delta_s + \Delta u^\delta_s \geq u^{-\delta}_t \\
		0 & \text{if } u^\delta_s + \Delta u^\delta_s < u^{-\delta}_t \\
		u^{-\delta}_t - (u^\delta_s + \Delta u^\delta_s) & \text{if } u^\delta_s + \Delta u^\delta_s \leq u^{-\delta}_t \\
		0 & \text{if } u^\delta_s + \Delta u^\delta_s > u^{-\delta}_t \\
	\end{cases}
$$

for $\delta \in \\{1, -1\\}$. Iterating through all $(u^{1}\_c)\_{c \in \Sigma}$ and $(u^{-1}\_c)\_{c \in \Sigma}$ for $i-1$, we compute $\text{DP}(i, \delta\_l, \delta\_r, \textsf{\*}, \textsf{\*})$. Suppose the sequence $q$ is formed by $Q$ contiguous subsequences, each of which comprise nucleotides of the same type, then 

\\[
	\min\_{S, \sigma} \Delta\_1(S, \sigma) = \min\_{\delta\_l, \delta\_r, (u^{1}\_c)\_{c \in \Sigma}, (u^{-1}\_c)\_{c \in \Sigma}}\text{DP}(Q, \delta\_l, \delta\_r, (u^{1}\_c)\_{c \in \Sigma}, (u^{-1}\_c)\_{c \in \Sigma}),
\\]

thus giving an algorithm to solve $\textsf{MFE}\_\{\Delta_1\}$. Its complexity is $O(n^{2\|\Sigma\|+1})$, implying that $\textsf{MFE}\_\{\Delta_1\} \in \|\Sigma\|-W[1]$, as desired. Note that this scheme also works for not-strictly complementary alphabets, i.e. given a character $s$, there might be two or more characters which can bind to it. In this case, the transition is somewhat more complicated but follows a similar line of argument.

For strictly complementary alphabets, an easy optimisation is possible: note that amongst the four numbers $u^1\_s$, $u^{-1}\_s$, $u^1\_t$, and $u^{-1}\_t$, exactly two of them are always $0$, thus instead of storing the values of all four, one only needs to store those of two, along with to binary variables indicating which two of the four being non-zero. This gives a dynamic programming algorithm with complexity $O(n^{\|\Sigma\|+1})$, e.g. $O(n^5)$ for Watson-Crick pairing, which is more practical.

## 7. Concluding remarks

We introduce Ising-based models for nucleic acid folding, and show various computational hardness results for the problem of evaluation a secondary structure and finding the minimum free energy. In one sentence, these models are computationally hard, if not intractable. Given the simplicity of base pair counting model and the the relevance of base pair stacking model, it seems to to be difficult to introduce angles to our current energy models whilst keeping the prediction tractable.

Whilst we have not shown the same for the problem of computing the partition function, in view of these theorems, it is reasonable to expect similar computational hardness results. On one hand, with careful check, one finds that the dynamic programming scheme introduced in Section 6 is non-redundant, so it is straightforward to derive another scheme to calculate the partition function for Ising BPC model. On the other hand, the hardness of partition function for Ising BPS model, whilst being expected, remains to be studied. It is worth noting that the problem of partition function for the original BPS model, to our best knowledge, has not been studied either.

Last but not least, we remind that as said in Section 1, it is imperative to take at least some of the geometric constraints into account. These results then beg the question of how to effectively introduce them into nucleic acid folding.

### References

{% bibliography --cited -f rna_ising_model %}