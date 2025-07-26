---
layout: post
title: "A slightly modified proof of $\\textsf{NP}$-hardness of RNA interaction problem"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\NP}{\textsf{NP}}
\newcommand{\A}{\textsf{A}}
\newcommand{\T}{\textsf{T}}
\newcommand{\G}{\textsf{G}}
\newcommand{\C}{\textsf{C}}
\newcommand{\U}{\textsf{U}}
$$

_Abstract_. This note is devoted to a slightly modified proof of $$\NP$$-hardness for $$\textsf{RNA-RNAi}$$ under Nussinov counting model, along the lines of that by Alkan et al. {% cite Alkan2006 -f modifled_NP_hardness_proof_of_rna_interaction %}, with all assertions rigorously proven. Then, I outline how the proof can be extended to reduce the alphabet used in the encoding of $$\NP$$-hard problems, a priori somewhat big, to the usual alphabet of four nucleotides $$\A$$, $$\U$$, $$\G$$ and $$\C$$.

## 1. Introduction 

### 1.1. Why do we care about RNAs?

Ribonucleic acid (RNA) is one of the three major biological macro-molecules, along with protein and DNA, that are essential to all known forms of life, and allow us to manipulate biological entities at molecular level, paving ways to potential medical research, biotechnology, and synthetic biology. In case of RNA, recent applications include message RNA vaccines {% cite Pardi2018 -f modifled_NP_hardness_proof_of_rna_interaction %}, RNAi-based therapeutics {% cite Setten2019 -f modifled_NP_hardness_proof_of_rna_interaction %}, RNA nanostructured devices {% cite Ohno2019, Grabow2014 -f modifled_NP_hardness_proof_of_rna_interaction %} to control cell fate {% cite Shibata2017 -f modifled_NP_hardness_proof_of_rna_interaction %}.

It is worth to mention that we have even yet to start exploring RNA's full potential, let alone to understand its capacity or to use it for our cause. Since the discovery of ribosomal RNA and transfer RNA in the late 1950s, a variety of non-coding RNA families have been discovered, which, to name a few, includes long non-coding RNA {% cite Statello2021 -f modifled_NP_hardness_proof_of_rna_interaction %}, circular RNA {% cite Liu2022 -f modifled_NP_hardness_proof_of_rna_interaction %}, and antisense RNA and CRISP RNA {% cite Ashwath2023 -f modifled_NP_hardness_proof_of_rna_interaction %}, each with potential for new therapies. Non-coding RNAs have shown to act as catalysts for chemical reactions, as gene regulators {% cite Statello2021 -f modifled_NP_hardness_proof_of_rna_interaction %}, and as DNA replicator {% cite Kuhnlein2021 -f modifled_NP_hardness_proof_of_rna_interaction %}, thereby are believed to play a vital role in the origin of life. We now know that at least $$76\%$$ of human genome is transcribed, yet only 1.2\% of which encodes proteins {% cite Lee2019 -f modifled_NP_hardness_proof_of_rna_interaction %}, and most of the rest, i.e. non-coding RNAs, have their functions yet to be completed unveiled, and applications to be discovered. Thus, there raises the need to understand properties and functions of RNAs.

On the one hand, advances in sequencing techniques have resulted in an exponential growth for the total length of sequences for the past decade {% cite Katz2022 -f modifled_NP_hardness_proof_of_rna_interaction %}. On the other hand, the nucleotides sequences alone do not determine RNAs' functionalities, as there have been evidences of RNAs having the same function yet not sharing the sequence {% cite Klosterman2002 -f modifled_NP_hardness_proof_of_rna_interaction %}. A possible explanation is that different nucleotide sequences can fold into similar structures, who ultimately determine RNA functionalities. For example, RNA viruses have a high mutation rate, and distant families of viruses show little resemblance in terms of nucleotide sequences, but instead show high conservation of secondary and tertiary structures.

Inferring structure of RNA by experimental methods, despite much progress in recent years, remains both arduous, time-consuming, and technically demanding. There have been attempts to combine experiment data with simulation, yet this approach has yet to overcome the limitation in both accuracy and throughput {% cite Liu2022 -f modifled_NP_hardness_proof_of_rna_interaction %}. Although RNA has its tertiary structure largely determined by its secondary structure {% cite Tinoco1990 -f modifled_NP_hardness_proof_of_rna_interaction %}, the same challenges persist, making prediction by computational methods being most viable approach.

### 1.2. What are RNA structures?

RNA is a chain of nucleotides (also known as bases) amongst the four bases - Adenine, Guanine, Uracil, Cytosine - and can be represented as a string on a 4-letter alphabet: $$\A$$, $$\U$$, $$\C$$, $$\G$$. This chain forms the so-called primary structure of RNA. 

Then, hydrogen bonds between nucleotides fold the RNA into what we call its secondary structure. For an RNA $$s = (s_i)_{1 \leq i \leq n}$$, a secondary structure is a set of _unordered_ pairs of indices $$q \subseteq \{1, ..., n\}^2$$. Such a secondary structure is compatible with $$s$$ if for any $$(i, j) \in q$$, then 
- $$(j, i) \in q$$, meaning the pairs are unordered;
- $$\{s_i, s_j\} \in \{\{\A, \U\}, \{G, C\}, \{G, U\}\}$$, meaning the base pair follows Watson-Crick base pairing model, with the exception of $$\G$$-$$\U$$, also known as wobble pairing; 
- for any $$(i', j') \in q$$, then $$\{i, j\} \cap \{i', j'\} = \emptyset$$ or $$\{i, j\} = \{i', j'\}$$, meaning each base is paired with at most one other base.

For a given secondary structure $$q$$, two base pairs $$(i, j)$$ and $$(i', j')$$ are said to be crossing if $$i < i' < j < j'$$ or $$i' < i < j' < j$$, thus form what is called a pseudoknot. If

The asymmetry of these bonds and the possibility of bonds between three or more nucleotides fold the RNA even further into its three-dimensional shape, that we call tertiary structure. In this report, we only concern ourselves with secondary structure for its algorithmic problems simple enough to be analysed. A natural problem is predicting secondary structures given sequences.

**Problem.** ($$\textsf{RNA Folding}$$) Given a sequence $$s$$, determine its secondary structure.

Like proteins, two or more RNAs can also have inter-molecular bonds, thus forming complexes, and thus another natural question is predicting how two RNAs can bond with each other, in addition to folding into themselves. Let $$s$$ and $$r$$ be two RNA sequences (e.g., an antisense RNA and its target), a joint secondary structure between $$s$$ and $$r$$ is a set of pairings where each nucleotide is paired with at most one another nucleotide either from $$s$$ or $$r$$, following Watson-Crick pairing and wobble pairing. In this case, two base pairs $$(i, j)$$ between $$s_i$$ and $$s_j$$ and $$(i', j')$$ between $$s_{i'}$$ and $$s_{j'}$$ satisfying $$i < i' < j < j'$$ or $$i' < i < j' < j$$ form what is called an _internal_ pseudoknot (of $$s$$). If, instead we have two base pairs $$(i, j)$$ between $$s_i$$ and $$s_j$$ and $$(i', j')$$ between $$r_{i'}$$ and $$r_{j'}$$ satisfying $$i < i' < j < j'$$ or $$i' < i < j' < j$$, then they form what is called an _external_ pseudoknot.

**Problem.** ($$\textsf{RNA-RNA interaction}$$ or $$\textsf{RNA-RNAi}$$) Given two sequences $$s$$ and $$r$$, determine their joint secondary structure.

### 1.3. Challenges of determining RNA structure

The problem of determining DNA secondary structure was solved first, since it is, up to a large extent, formed by base pairing between two strands of nucleotides. Thus predicting it simply amounts to aligning two strands of nucleotides, which can be done with, e.g., Needleman-Wunsch algorithm {% cite Needleman1970 -f modifled_NP_hardness_proof_of_rna_interaction %} or Smith-Waterman algorithm {% cite Smith1981 -f modifled_NP_hardness_proof_of_rna_interaction %}. Historically, the problem of determining protein secondary structure was largely solved with great accuracy, simply because protein secondary structure can be categorised into three motifs: $$\alpha$$-helix, $$\beta$$-sheet, and random coil. Determining a segment has which motif can be done by aligning it with multiple segments of amino acids with known structure, e.g., determined experimentally {% cite Pirovano2010 -f modifled_NP_hardness_proof_of_rna_interaction %}.

Determining RNA secondary structure, however, poses significantly more difficult challenges, some of which can be attributed to the following reasons:

- Unlike DNA, RNA structure is formed by a sequence folding into itself, thus have more degrees of freedom. And unlike peptide bonds between amino acids in proteins, hydrogen bonds of nucleotides acids are weaker, thus RNA is much more unstable, and its structure can change depending on variables of environments, e.g. temperature. It is also a challenge when determining RNA structure experimentally.
- The lack of RNA structural data, which stems from difficulties in determining RNA structure via experiments.
- RNA is transcribed from DNA with mRNA, and it does not wait for the whole nucleotide sequence to be transcribed before folding. Instead, the RNA folds as it is being transcribed, in a process that we call co-transcriptional folding. It is shown that RNA structure is not determined by its static folding alone, but also its folding pathway, i.e. how it folds co-transcriptionally {% cite Lai2013 -f modifled_NP_hardness_proof_of_rna_interaction %}.

### 1.4. Goal and outline

This note is devoted to a slightly modified proof of $$\NP$$-hardness for $$\textsf{RNA-RNAi}$$ under Nussinov counting model, along the lines of that by Alkan et al. {% cite Alkan2006 -f modifled_NP_hardness_proof_of_rna_interaction %}, with all assertions rigorously proven. Then, I outline how the proof can be extended to reduce the alphabet used in the encoding of $$\NP$$-hard problems, a priori somewhat big, to the usual alphabet of four nucleotides $$\A$$, $$\U$$, $$\G$$ and $$\C$$.

## 2. Preliminaries

### 2.1. Minimum Free Energy

Despite its instability, experiments showed that RNA secondary structure tends to follow the laws of thermodynamic, giving birth to Minimum Free Energy approach. In particular, since the pioneer work of Anfinsen {% cite Anfinsen1973 -f modifled_NP_hardness_proof_of_rna_interaction %}, it is generally accepted that given an energy model of secondary structures, RNA generally folds into a structure with minimum free energy. Even if RNA can dynamically change its structure, Matthews et al. {% cite Mathews1999 -f modifled_NP_hardness_proof_of_rna_interaction %} shows that most experimentally measured RNA structures have their free energy in the vicinity of the minimum, and moreover, Wuchty et al. {% cite Wuchty1999 -f modifled_NP_hardness_proof_of_rna_interaction %} proposed an algorithm to iterate through _all_ such possible structure, given a vicinity.

For formality's sake, denote $$\mathcal{S}$$ and $$\mathcal{Q}$$ the space of sequences and secondary structures, an RNA free energy model comprises of two parts: a collection of structural features, denoted by $$\boldsymbol{f} := (f_1, f_2, ..., f_p) : \mathcal{S} \times \mathcal{Q} \rightarrow \mathbb{R}^p$$, a parameter set $$\boldsymbol{\theta} \in \mathbb{R}^p$$ (although the functions $$f_i$$ are commonly integer-valued), where $$p$$ denotes the number of features. Then, the energy associated with a sequence $$s$$ and a secondary structure $$q$$ compatible with $$s$$ is given by 

$$
	E(q) = \langle \boldsymbol{f}(s, q), \boldsymbol{\theta} \rangle,
$$

and predicting a structure of $$s$$ amounts to finding a $$q^*$$ compatible with $$s$$ such that

$$
	\boldsymbol{f}(s, q^*) = \Delta(s) = \min_{q \text{ compatible with } s} E(q).
$$


### 2.2. A first example: Nussinov algorithm
For instance, the first example was considered by Nussinov {% cite Nussinov1980 -f modifled_NP_hardness_proof_of_rna_interaction %}, where $$p = 3$$, and for a given $$q$$, $$f_1(q)$$, $$f_2(q)$$, and $$f_3(q)$$ are the number of $$\A$$-$$\U$$ pairs, $$\G$$-$$\C$$ pairs, and $$\A$$-$$\U$$ pairs, respectively. For simplicity sake, she also considered $$\boldsymbol{\theta} = (-1, -1, -1)$$. The model thus counts the number of base pairs, and prediction amounts to maximising it.

To predict secondary structure, a common assumption in Minimum Free Energy approach is that no pseudoknots are allowed. This is because then we can divide the RNA into multiple small segments and predict the structure of each segment. In particular, denote $$(s_k)_{i \leq k \leq j}$$ by $$s[i..j]$$, and $$f(i, j)$$ be the maximum number of base pairs in $$s[i..j]$$. When $$i \geq j$$, either $$s[i..j]$$ is empty or has one nucleotides, thus cannot form any base pairs, so $$f(i j) = 0$$. When $$i < j$$, we have three cases:
- either $$i$$ or $$j$$ is not a part of the pairing: then $$f(i, j) = f(i+1, j)$$ or $$f(i, j) = f(i, j-1)$$.
- $$(i, j)$$ forms a base pair: then $$f(i, j) = f(i+1, j-1) + 1$$;
- $$(i, j)$$ do not form a base pair: then there exists $$i \leq k < j$$ such that no base pairs crosses $$k$$, and thus $$f(i, j) = \max_{i \leq k < j} f(i, k) + f(k+1, j)$$.

We have
<div align="center">
$$
	f(i, j) = \max
	\begin{cases}
		f(i+1, j) \\
		f(i, j-1) \\
		f(i-1, j+1) + p_{\{q_i, q_j\}} \text{ if } q_i \text{ and } q_j \text{ can form a base pair} \\
		\max_{i < k < j} f(i, k) + f(k+1, j)
	\end{cases}.
$$
</div>
where the answer is $$f(1, n)$$ and the table of $$f$$ is calculated in $$O(n^3)$$.

Then, given the memoisation table $$f$$, one can trace back to find an optimal and compatible secondary structure. Note that such a structure is not unique. For instance, let $$s = \U\A\U\U\C\U\G\A\U\G$$, we have the following table, and some optimal structures.

|      | $$\U$$ | $$\A$$ | $$\U$$ | $$\U$$ | $$\C$$ | $$\U$$ | $$\G$$ | $$\A$$ | $$\U$$ | $$\G$$ |
|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| $$\U$$ |   0  |   1  |   1  |   1  |   1  |   2  |   2  |   3  |   4  |   4  |
| $$\A$$ |   -  |   0  |   0  |   0  |   0  |   1  |   2  |   2  |   3  |   3  |
| $$\U$$ |   -  |   -  |   0  |   0  |   0  |   1  |   2  |   2  |   3  |   3  |
| $$\U$$ |   -  |   -  |   -  |   0  |   0  |   1  |   2  |   2  |   3  |   3  |
| $$\C$$ |   -  |   -  |   -  |   -  |   0  |   1  |   2  |   2  |   3  |   3  |
| $$\U$$ |   -  |   -  |   -  |   -  |   -  |   0  |   0  |   1  |   1  |   1  |
| $$\G$$ |   -  |   -  |   -  |   -  |   -  |   -  |   0  |   1  |   1  |   1  |
| $$\A$$ |   -  |   -  |   -  |   -  |   -  |   -  |   -  |   0  |   1  |   1  |
| $$\U$$ |   -  |   -  |   -  |   -  |   -  |   -  |   -  |   -  |   0  |   0  |
| $$\G$$ |   -  |   -  |   -  |   -  |   -  |   -  |   -  |   -  |   -  |   0  |

<div align="center">Table 1. Memoisation table for $s = \U\A\U\U\C\U\G\A\U\G$.</div>

![](/littlenotebook/assets/img/rna_interaction/rna.png##center)
<div align="center">Figure 1. An optimal structure for $s$.</div>

### 2.3. Summary of some $$\NP$$-hardness results.

Compared with $$\textsf{RNA Folding}$$, the problem $$\textsf{RNA-RNAi}$$, despite being no less important, to my best knowledge, has fewer $$\NP$$-completeness results, which can be attributed to the fact that unlike $$\textsf{RNA Folding}$$, even in the absence of (both internal and external) pseudoknots, $$\textsf{RNA-RNAi}$$ is already $$\NP$$-hard. In particular, 

- Alkan et al. {% cite Alkan2006 -f modifled_NP_hardness_proof_of_rna_interaction %} showed that $$\textsf{RNA-RNAi}$$ is $$\NP$$-hard under Nussinov generalised model and general BPS model with fixed weight function $$\gamma$$.
- Mneimneh {% cite Mneimneh2009 -f modifled_NP_hardness_proof_of_rna_interaction %} around the same time, showed that $$\textsf{RNA-RNAi}$$ is $$\NP$$-hard in a more general setting where a weight function $$\gamma : (s \cup r)^2 \to \mathbb{R}$$ is fixed and the energy is given by 

$$
    E(q) = \frac{1}{2}\sum_{(i, j) \in q} \gamma(i, j).
$$

where the factor $$\frac{1}{2}$$ is to avoid double counting. In fact, he showed that $$\NP$$-hardness remains even if $$\gamma$$ is constant, i.e. the weights are uniform. Later, the author extends this result to three or more RNA sequences {% cite Mneimneh2015 -f modifled_NP_hardness_proof_of_rna_interaction %}.

## 3. Polynomial reduction to $$\NP$$-complete problems

### 3.1. Choice of $$\NP$$-hard problem, construction of the encoding

First, we consider a generalisation of $$\textsf{RNA-RNAi}$$, where the sequences $$s$$ and $$r$$, instead of being over the alphabet $$\{\A, \U, \G, \C\}$$, will be considered over an extended alphabet $$\{a, b, c, d, e, f, u, p, q, x, y\}$$. For convenience's sake, for a character $$u$$, denote $$\overline{u}$$ its complementary character, and we set $$\overline{a} = b$$, $$\overline{c} = d$$, $$\overline{e} = f$$, $$\overline{p} = q$$, $$\overline{x} = y$$, and $$\overline{n} = n$$. By abuse of notation, for a sequence $$s$$, we denote $$\overline{s}$$ its complementary sequence.

The idea of the proof is to reduce from the problem of longest common subsequence of multiple binary strings $$\textsf{mLCS}$$, which is known to be $$\NP$$-complete {% cite Maier1978 -f modifled_NP_hardness_proof_of_rna_interaction %}. We recall the problem $$\textsf{mLCS}$$: given a set of $$m$$ strings $$s_1, \cdots, s_m$$ of equal length $$n$$ over the alphabet $$\{0, 1\}$$ and an integer $$k$$, decide if there exists a string $$t$$ of length $$k$$ that is a common substring of all $$n$$ given strings. Duplicating the first string if necessary, without loss of generality, we may assume that $$m$$ is even and $$m \geq 4$$.

Given an instance of $$\textsf{mLCS}$$ of $$2m$$ strings for some $$m \geq 1$$, let $$s_{i, j} = p$$ if the $$j^{\text{th}}$$ character of $$s_i$$ is $$0$$, and $$x$$ otherwise. For $$1 \leq i \leq 2m$$:
- if $$i$$ is odd, then we construct

$$
    A_i = a s_{i, 1} a s_{i, 2} a \cdots a s_{i, n} a
$$

and  

$$
    B_i = \overline{A_i} = b \overline{s_{i, 1}} b \overline{s_{i, 2}} b \cdots b \overline{s_{i, n}} b
$$

- if $$i$$ is even, then we construct

$$
    A_i = a \overline{s_{i, n}} a \overline{s_{i, n-1}} a \cdots a \overline{s_{i, 1}} a
$$

and 

$$
    B_i = \overline{A_i} = b s_{i, n} b s_{i, n-1} b \cdots b s_{i, 1} b
$$


Next, we define

$$
	s = (c u^k A_1 d c) (c A_2 A_3 d c) (c A_4 A_5 d c) \cdots (c A_{2m-2} A_{2m-1} d c) (c A_{2m} u^k d) d^m,
$$

and 

$$
	r = (e B_1 B_2 f e) (e B_3 B_4 f e) (e B_5 B_6 f e) \cdots (e B_{2m-3} B_{2m-2} f e) (e B_{2m-1} B_{2m} f) f^{m-1}.
$$


Simple calculation shows that $\|D_i\| = \|E_i\| = 2n+1$ for all $$i$$, so $\|s\|$ and $\|r\|$ are in polynomial of $$n$$, $$m$$, and $$k$$. We will show below that this encoding is correct.

### 3.2. An alternative construction and two properties

The explicit construction above makes it easy to argue about the length of $$s$$ and $$r$$, but it is not so convenient to work with. Instead, we will consider an alternative construction, as follow. For the sake of convenience, in the context of this subsection only, we will denote by $$n$$ an integer that does not necessarily have the same meaning as above.

Let $$S_1, S_2, \cdots, S_n$$ not containing $$c$$ or $$d$$, $$C_1 = c S_1 d$$, and $$C_{i+1} = c S_{i+1} d c C_i d$$ for $$2 \leq i \leq n$$.

First, it is easy to see that for $$S_{m+1} = u^k A_1$$, $$S_{m+1-i} = A_{2i-2} A_{2i-1}$$ for $$2 \leq i \leq m$$, and $$S_{1} = A_{2m} u^k$$ then $$C_{m+1} = s$$. Similarly, replacing $$c$$ with $$e$$ and $$d$$ with $$f$$, for $$S_{m+1-i} = B_{2i-1} B_{2i}$$ for all $$1 \leq i \leq m$$ then $$C_{m} = r$$. So indeed this gives an alternative construction to $$s$$ and $$r$$.

Next, we have an easy observation.

**Lemma 1.** For any $$1 \leq i \leq n$$, in $$C_i$$, all $$c$$'s can be paired with $$d$$'s without introducing pseudoknots.

_Proof._ Since $$S_1, S_2, \cdots, S_n$$ do not contain $$c$$ or $$d$$, we can assume that $$S_i$$ is empty for all $$i$$. Then, this follows easily by induction. <span style="float:right;">$$\square$$</span>

The following observation, however, plays a crucial role later, and motivates this rather unnatural construction.

**Lemma 2.** Let $$1 \leq i \leq n$$.
- In $$cC_id$$ (notice the additional $$c$$ and $$d$$ at the beginning and at the end), the first $$c$$ _must be paired_ with the last $$d$$.
- In $$C_i$$, for any $$1 \leq j \leq i$$, a nucleotide of $$S_j$$ can only be paired with another nucleotide of $$S_j$$. Moreover, an nucleotide not in $$C_i$$ cannot be paired with a nucleotide in $$C_i$$.

_Proof._ 
- We argue by induction. For $$i = 1$$, we have $$cC_1 d = ccS_1dd$$, and the the statement is true because we forbid pseudoknots.
		
    Now suppose the statement holds for $$i-1$$, we have $$cC_i d = c c S_i d c C_{i-1} d d$$. By induction hypothesis on $$i-1$$, in $$c C_{i-1} d$$, the first $$c$$ must be paired with the last $$d$$, so we can reduce to string $$c c S_i d d$$, and by the forbiddance of pseudoknots, the statement follows.
		
- We argue by induction. For $$i = 1$$, there is nothing to prove, so now suppose the statement holds for $$i-1$$. We have that $$C_i = c S_i d c C_{i-1} d$$. Let us focus on the first $$d$$, and see with which $$c$$ it can be paired.
	-  if $$d$$ is paired with a $$c$$ lying before it, then there is only one such $$c$$, namely the first $$c$$.
			
        This bond splits $$C_i$$ into two halves, $$cS_i d$$ and $$c C_{i-1} d$$. By the previous point of the lemma, in $$c C_{i-1} d$$, the first $$c$$ must be paired with the last $$d$$. The statement then follows by the forbiddance of pseudoknots and by induction hypothesis.
    - if $$d$$ is paired with a $$c$$ lying after it. By induction hypothesis, $$d$$ cannot be paired with any nucleotide in $$C_{i-1}$$, so it can only be paired with the $$c$$ next to it.
			
        In this case, by a similar argument, the first $$c$$ of $$C_i$$ cannot be paired with any nucleotide in $$C_{i-1}$$, so it can only be paired with the last $$d$$. By the forbiddance of pseudoknots, the second point of the statement follows.
        
        Finally, that leaves $$S_i$$ and $$C_{i-1}$$. By induction hypothesis, no nucleotide of $$S_i$$ can be paired with any nucleotide in $$C_{i-1}$$, so the nucleotides of $$S_i$$ can only be paired with themselves. The first point of the statement follows.
<span style="float:right;">$$\square$$</span>

### 3.3. Proof that if $$\textsf{mLCS}$$ has solution, then $$s$$ and $$r$$ are perfectly paired.

Suppose the answer to $$\textsf{mLCS}$$ instance is yes. Let $$t$$ be a common substring of all $$s_i$$'s, and denote $$t_i$$ be $$s_i$$ after deleting $$t$$. Moreover, denote $$A^t_i$$ (resp. $$B^t_i$$) the nucleotides in $$A_i$$ (resp. $$B_i$$) corresponding with $$t$$, and $$A^s_i$$ (resp. $$B^s_i$$) the nucleotides in $$A_i$$ (resp. $$B_i$$) corresponding with $$t_i$$. Let us see how $$s$$ and $$r$$ can be paired up.
- By Lemma 1, all $$c$$'s can be paired with all $$d$$'s. Similarly, all $$e$$'s can be paired with all $$f$$'s.
- Now we pair $$u^k$$ with $$A^t_1$$. Note that they all lie to the left of the bond $$c-d$$, so no internal pseudoknots are created. Then, we pair $$A^s_1$$ with $$B^s_1$$.
- We then pair $$B^t_1$$ with $$B^t_2$$, and continue. More generally, we pair $$B^t_{2i-1}$$ with $$B^t_{2i}$$, $$B^s_{2i}$$ with $$A^s_{2i}$$, $$A^t_{2i}$$ with $$A^t_{2i+1}$$, then $$A^s_{2i+1}$$ with $$B^s_{2i+1}$$.
- Finally, after we pair $$B^s_{2m}$$ with $$A^s_{2m}$$, we pair $$A^t_{2m}$$ with $$u^k$$.
Thus, all the nucleotides are paired.

### 3.4. Proof that if $$s$$ and $$r$$ are perfectly paired, $$\textsf{mLCS}$$ has solution.

First, since we forbid external pseudoknots, the $$i^{\text{th}}$$ $$a$$ in $$s$$ can only be paired with the $$i^{\text{th}}$$ $$b$$ in $$r$$.

By Lemma 2, if we consider internal bonds, then the nucleotides in $$A_1$$ can only be paired with $$u^k$$, and we denote $$T$$ the positions in the string $$s_1$$ corresponding to these nucleotides. The remaining nucleotides of $$A_1$$ then must be paired with the corresponding of $$B_1$$, and then the nucleotides corresponding to $$T$$ in $$B_1$$ must be paird with those in $$B_2$$, by Lemma 2.

The argument continues. More generally, the nucleotides corresponding to $$T$$ in $$B_{2i-1}$$ must be paired with those in $$B_{2i}$$. Then, the nucleotides not corresponding to $$T$$ in $$B_{2i}$$ must be paired with those in $$A_{2i}$$. Then since $$A_{2i}$$ can only be with $$A_{2i+1}$$, the nucleotides corresponding to $$T$$ in $$A_{2i}$$ must be paired with those in $$A_{2i+1}$$. The remaining of $$A_{2i+1}$$ then must be paired with the corresponding in $$B_{2i+1}$$.

Finally, the nucleotides corresponding to $$T$$ in $$A_{2m}$$ must be paired with $$u^k$$. $$T$$ encodes the same positions in all strings $$s_i$$, and at these positions, all the characters in $$s_i$$'s are equal, thus $$T$$ is a solution to $$\textsf{mLCS}$$.

Verifying all nucleotides are paired can be done in polynomial time as $$\|s\|$$ and $$\|r\|$$ are of polynomial length of $$n$$, $$m$$, and $$k$$, so in fact, we have that checking if a pairing is complete is $$\NP$$-complete. The reduction is complete, and we have that $$\textsf{RNA-RNAi}$$ is $$\NP$$-hard.

## 4. Comparison with the proof by Alkan et al.

Reading the original proof by Alkan et al. {% cite Alkan2006 -f modifled_NP_hardness_proof_of_rna_interaction %}, we will find that it differs from our proof here. In particular, they assumed that $$m$$ is odd and at least $$1$$. Then, their construction of $$s$$ and $$r$$ are given by

$$
	s = u^k A_1 (c^1 A_2 A_3 d^1) (c^2 A_4 A_5 d^2) \cdots (c^m A_{2m} A_{2m+1} d^m) 
$$

and

$$
	r = (e^1 B_1 B_2 f^1) (e^2 B_3 B_4 f^2) \cdots (e^m B_{2m-1} B_{2m} f^m) B_{2m+1} u^k.
$$


So why do we consider a different proof?

The arguments in that by Alkan et al. {% cite Alkan2006 -f modifled_NP_hardness_proof_of_rna_interaction %} are exactly the same as ours, and in particular, they claim that
> Nucleotides $$c$$, $$d$$ only occur in $$s$$ and form valid bonds only with each other. **Because we do not allow internal pseudoknots, each $$c^i$$ block will be bonded with the $$d^i$$ block.** Similarly, nucleotides $$e$$, $$f$$ only occur in $$r$$ and form valid bonds only with each other. Again, because we do not allow internal pseudoknots, each $$e^i$$ block will be bonded with the $$f^i$$ block.

But the authors gave no proof of the fact that "each $$c^i$$ block will be bonded with the $$d^i$$ block," and I found it difficult to prove. Still, the idea of using repeated $$c$$'s and $$d$$'s plays a crucial role later, as we try to reduce to canonical nucleotides.

## 5. Toward canonical nucleotides

In this proof, we have used 11 nucleotides, 10 of which form 5 complementary pairs, namely $$a-b$$, $$c-d$$, $$e-f$$, $$x-y$$, $$p-q$$, in addition to a wild-card nucleotide $$u$$ which can be paired with $$x$$, $$y$$, $$p$$, or $$q$$. Aiming at biological relevance, in this section, we will reduce them to four canonical nucleotides $$\A$$, $$\U$$, $$\G$$, and $$\C$$.

### 5.1. 5 pairs to 4 pairs

Our first observation is that we can replace $$e$$ by $$c$$ and $$f$$ by $$d$$. Indeed, looking at our construction, we see that for any $$1 \leq i \leq 2m$$, if $$i$$ is odd (resp. even), then there are no nucleotides between $$A_i$$ and $$A_{i+1}$$ (resp. $$B_i$$ and $$B_{i+1}$$). By $$a-b$$ pairings between $$A_i$$ and $$B_i$$ and the forbiddance of external pseudoknots, this means that supposed all nucleotides are paired, a $$c$$ or a $$d$$ in $$s$$ cannot be paired with another $$c$$ or $$d$$ in $$r$$, and vice versa.

### 5.2. 4 pairs to 3 pairs

This is perhaps the most difficult part. The idea is using a new type of nucleotide pairing, e.g. $$g-h$$ such that $$u$$ can be paired with $$c$$, $$d$$, $$g$$, or $$h$$. then using $$g-h$$ and $$c-d$$, we can encode $$x-y$$ and $$p-q$$. To do this, we have to modify our construction.

Given an instance of $$\textsf{mLCS}$$ of $$2m+1$$ strings for some $$m \geq 1$$, let $$S_{i, j} = gcc$$ if the $$j^{\text{th}}$$ character of $$s_i$$ is $$0$$, and $$cgc$$ otherwise. For a string $$S$$, denote $$S^T$$ the reverse of $$S$$. For $$1 \leq i \leq m$$, we define

$$
	A_{2i} = a c S_{2i, 1} a c S_{2i, 2} a \cdots a c S_{2i, n} a, 
$$


$$
	A_{2i+1} = a \overline{S_{2i+1, n}}^T d a \overline{S_{2i+1, n-1}}^T d a \cdots a \overline{S_{2i+1, 1}}^T d a, 
$$


$$
	B_{2i-1} = b c S_{2i-1, n}^T b c S_{2i-1, n-1}^T b \cdots b c S_{2i-1, 1}^T b, 
$$

and 

$$
	B_{2i} = b \overline{S_{2i, 1}} d b \overline{S_{2i, 2}} d b \cdots b \overline{S_{2i, n}} d b.
$$

Finally, we define (without $$c$$ or $$d$$ as above)

$$
	A_1 = a \overline{S_{1, n}}^T a \overline{S_{1, n-1}}^T a \cdots a \overline{S_{1, 1}}^T a,
$$


$$
	A_{2m} = a S_{2m, 1} a S_{2m, 2} a \cdots a S_{2m, n} a,  
$$

and we set  

$$
	s = (c u^{3k} A_1 d c) (c A_2 A_3 d c) (c A_4 A_5 d c) \cdots (c A_{2m-2} A_{2m-1} d c) (c A_{2m} u^k d) d^m,
$$

and 

$$
	r = (c B_1 B_2 d c) (c B_3 B_4 d c) (c B_5 B_6 d c) \cdots (c B_{2m-3} B_{2m-2} d c) (c B_{2m-1} B_{2m} d) d^{m-1}.
$$


In view of Lemma 2, in order to prove that our construction is correct, we have to show that suppose all nucleotides are paired, then for $$1 \leq i \leq m$$
1. in $$c A_{2i} A_{2i+1} d$$, the first $$c$$ must be paired with the last $$d$$;
2. $$S_{2i, j}$$ must be paired with $$\overline{S_{2i+1, j}}^T$$.
The first point is straightforward, as when we ignore $$g$$, $$h$$, $$a$$,s and $$b$$, then for some $$\ell$$, $$c A_{2i} A_{2i+1} d = c^\ell d^\ell$$, thus forces the first $$c$$ must be paired with the last $$d$$.

The second point is an extension of the first point, and can be proved similarly. It suffices to show that the first $$c$$ of substring $$S_{2i, j}$$ must be paired with the last $$d$$ of substring $$\overline{S_{2i+1, j}}^T$$, and we argue by induction. The case $$j = n$$ is the first point above. Then assuming this is true for all $$\ell > j$$, ignoring all the $$a$$'s, consider the substring

$$
c S_{2i, j} \left(c S_{2i, j+1} \cdots c S_{2i, n} \overline{S_{2i+1, n}}^T d \cdots \overline{S_{2i+1, j+1}}^T d\right) \overline{S_{2i+1, j}}^T d,
$$

where in the parentheses, the first $$c$$ must be paired with the last $$d$$, so we can ignore the part inside the parentheses. Then this reduces to the first point, and the induction is complete. A similar statement goes analogously for $$B_{2i-1}$$ and $$B_{2i}$$.

Then, if all nucleotides are paired, $$S_{2i, j}$$ is perfectly paired with $$\overline{S_{2i+1, j}}^T$$, meaning $$S_{2i, j} = S_{2i+1, j+1}$$, and the remaining of the argument in subsection \ref{subsecition: RNA-RNAi to mLCA} follows similarly.

### 5.3. 3 pairs to 2 pairs

We notice that using $$g-h$$ and $$c-d$$, we can also encode $$a-b$$, by replacing $$a$$ with $$A = ccg$$ and $$b$$ with $$\overline{A}$$. We also have to modify our construction as follow: for $$1 \leq i \leq m$$, we define 

$$
	A_{2i} = c A c S_{2i, 1} c A c S_{2i, 2} a \cdots c A c S_{2i, n} c A c, 
$$


$$
	A_{2i+1} = d A d \overline{S_{2i+1, n}}^T d A d \overline{S_{2i+1, n-1}}^T d A d \cdots A d \overline{S_{2i+1, 1}}^T d A d, 
$$


$$
	B_{2i-1} = c B c S_{2i-1, n}^T b C b S_{2i-1, n-1}^T b \cdots c B c S_{2i-1, 1}^T B c, 
$$

and 

$$
	B_{2i} = d B \overline{S_{2i, 1}} d B d \overline{S_{2i, 2}} d B d \cdots d B \overline{S_{2i, n}} d B d.
$$

We also have to redefine $$A_1$$ and $$A_{2m}$$, as follow:

$$
	A_1 = A \overline{S_{1, n}}^T A \overline{S_{1, n-1}}^T A \cdots A \overline{S_{1, 1}}^T A,
$$


$$
	A_{2m} = A S_{2m, 1} A S_{2m, 2} A \cdots A S_{2m, n} A,  
$$

The proof goes exactly the same as above, except now the $$A$$ in $$c A c S_{2i, j}$$ cannot be paired with the $$A$$ in $$\overline{S_{2i+1, j}}^T d A d$$, thus must be paired with the corresponding $$B$$.

### 5.4. Wild-card nucleotide

This last part starts by noticing that assuming all nucleotides are paired, in $$s$$, the only nucleotides paired with $$u$$ are the ones in $$\overline{S_{1, j}}^T$$ and in $$S_{2m, j}$$ for $$1 \leq j \leq n$$, which comprise of only $$d$$ and $$h$$. 

Now we can let $$c = \U$$, $$d = \A$$, $$g = \C$$, $$h = \G$$, and $$u = \U$$ (one can compare the nucleotide chain $$c-d-u-h-g$$ with $$\U-\A-\U-\G-\C$$). This unfortunately introduces a new possible type of bond between $$c$$ and $$h$$, and we need to prove that this does not affect our proof.

To do this, we make one final modification to our construction of $$s$$ and $$r$$. Instead of assigning $$A = ccg$$, we assign $$A = cccgc$$. For $$1 \leq i \leq 2m$$ and $$1 \leq j \leq n$$, if the $$j^{\text{th}}$$ character of $$i^{\text{th}}$$ is $$0$$, we set $$S_{i, j} = cgcccc$$ instead of $$gcc$$, and if the character is $$1$$, we set $$S_{i, j} = ccgcc$$ instead of $$cgc$$. We also replace $$u^{3k}$$ with $$u^{5k}$$.

Now consider $$c \overline{S_{1, j}}^T c$$, which is either $$cdhdddc$$, $$cddhddc$$, or $$cdddhdc$$, and we see clearly that there can be no bound between $$c$$ and $$h$$, assuming all nucleotides are paired.

There can be bound between $$c$$ and $$h$$ in case of bonding $$S_{2i, j}$$ with $$\overline{S_{2i+1, j}}^T$$, we need to verify that even if it be the case, we still have that $$S_{2i, j}$$ and $$\overline{S_{2i+1, j}}^T$$ are perfectly paired if and only if $$S_{2i, j} = S_{2i+1, j}$$. Table 2 proves exactly that, and our reduction is complete.

<table style="text-align:center;">
 <tr>
  <td colspan="3" rowspan="3">&nbsp;</td>
  <td colspan="3">$$S_{2i,j}$$</td>
 </tr>
 <tr>
  <td>$$cgccc$$</td>
  <td>$$ccgcc$$</td>
  <td>$$cccgc$$</td>
 </tr>
 <tr>
  <td>$$\U\C\U\U\U$$</td>
  <td>$$\U\U\C\U\U$$</td>
  <td>$$\U\U\U\C\U$$</td>
 </tr>
 <tr>
  <td rowspan="3">$$\overline{S_{2i+1,j}}^T$$</td>
  <td>$$dhddd$$</td>
  <td>$$\A\G\A\A\A$$</td>
  <td>$$\textsf{_____}$$</td>
  <td>$$\textsf{__*__}$$</td>
  <td>$$\textsf{___*_}$$</td>
 </tr>
 <tr>
  <td>$$ddhdd$$</td>
  <td>$$\A\A\G\A\A$$</td>
  <td>$$\textsf{_*___}$$</td>
  <td>$$\textsf{_____}$$</td>
  <td>$$\textsf{___*_}$$</td>
 </tr>
 <tr>
  <td>$$dddhd$$</td>
  <td>$$\A\A\A\G\A$$</td>
  <td>$$\textsf{_*___}$$</td>
  <td>$$\textsf{__*__}$$</td>
  <td>$$\textsf{_____}$$</td>
 </tr>
</table>

<div align="center">Table 2. Possible bonding between $S_{2i, j}$ and $\overline{S_{2i+1, j}}^T$, with $\textsf{_}$ denotes a match and $\textsf{*}$ denotes a mismatch.</div>

### References

{% bibliography --cited --file modifled_NP_hardness_proof_of_rna_interaction %}