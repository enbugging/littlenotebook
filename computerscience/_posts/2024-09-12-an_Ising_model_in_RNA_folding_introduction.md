---
layout: post
title: "An Ising model in RNA folding: Introduction"
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
\newcommand{\U}{\textsf{U}}
$$

This is part 1 of my Bachelor thesis during my first year at Ulm, carried out at Hamilton Institute, Maynooth University and under supervision of Prof. Damien Woods, titled "Computational complexity of an Ising model in pseudoknotted nucleic acid folding", to be published.

## 1. Introduction
Since the early studies on nucleic acids, there has been much effort going into harvesting its power. With RNAs, we have been achieved new ways to make vaccine {% cite Pardi2018 -f rna_ising_model %}, edit genome {% cite Deltcheva2011 -f rna_ising_model %}, and make nanostructured devices {% cite Ohno2019, Kai2021, Grabow2014 -f rna_ising_model %} to control cell fate {% cite Shibata2017 -f rna_ising_model %}. With DNA, we have found treatments for some chronic disease {% cite Sussman2024 -f rna_ising_model %}, done computations {% cite Fan2020 -f rna_ising_model %}, and stored data {% cite Ceze2019 -f rna_ising_model %}. And yet, there is so much more that we do not know about them. We now know that at least 76\% of human genome is transcribed, yet only 1.2\% of which encodes proteins {% cite Lee2019 -f rna_ising_model %}, and most of the rest, i.e. non-coding RNAs, have their functions yet to be completed unveiled {% cite Statello2021 -f rna_ising_model %}, and applications to be discovered. For example, it is only recently that tRNAs, most commonly known as encoding of amino acid sequences, can assemble into a DNA replicator {% cite Kuhnlein2021 -f rna_ising_model %}. Given how relatively recent these achievements are, it is expected that there is more potential to explore.

Notwithstanding, it is not fair to say that these successes are entirely product of this century, as they are built upon decades of studying nucleic acids, e.g. how they form, how their structures are related to their functions, how to predict their structures, and how to design nucleic acids that form a particular target structure.

We now understand that they start out as sequences of nucleic acids, also known as bases, e.g. $\A$, $\T$, $\G$, $\C$ for DNA, and $\A$, $\U$, $\G$, $\C$ for RNA, lying on phosphate backbones. The hydrogen bonds between the bases form what is known as secondary structure. Many of such structural elements form recognisable three-dimensional motifs such as helices, major and minor grooves, and quadruplexes, together forming what is called a tertiary structure. Finally, several such structures can, in some situations, self-assemble into larger complexes called quaternary structures.

On the one hand, advanced in nucleic acid sequencing have allowed exponential growth in the number of nucleic acid sequences we determined {% cite Katz2022 -f rna_ising_model %}. On the other hand, the next step is significantly more difficult. We know that the structure of a nucleic acid also plays critical roles in gene regulation {% cite Mandal2004 -f rna_ising_model %}, protein-DNA recognition {% cite Rohs2009 -f rna_ising_model %}, DNA replication {% cite Sims1980 -f rna_ising_model %}, and other processes. Moreover, the sequence alone does not determine a nucleic acid's functionalities, as there have been evidences of RNA polymers having the same function yet not sharing the same sequence {% cite Klosterman2002, Tamura2004 -f rna_ising_model %}. Yet despite much progress, even determining secondary structure remains arduous, time-consuming, and technically demanding from an experimental point of view {% cite VonLohneysen2024 -f rna_ising_model %}.

An alternative approach is introducing a thermodynamic model which associates to each secondary structure a free energy, and the ones minimising free energy are considered most favourable, thus most likely. The models can be simple as maximising the number of paired nucleotides {% cite Nussinov1980 -f rna_ising_model %}, or complicated as a function of over 7600 parameters {% cite Mathews1999 -f rna_ising_model %}, but they all consider secondary structures as mere pairs of bases respecting some pairing, e.g. Watson-Crick pairing or wobble pairing.

But, these models are certainly not sufficient for practical purposes. One of the reasons is that they consider few geometric constraints, e.g. the minimum length of hairpins in case of nearest neighbour model {% cite Mathews1999 -f rna_ising_model %}, which is problematic because in order for two bases to pair, they must be close and align with each other. A predicted secondary structure by these models thus might not be realistic; this problem cannot be ignored, as it is known that these phenomena do affect the ways wherein nucleic acids can bind together {% cite Wang1979 -f rna_ising_model %}, thus one must take it into consideration when designing nucleic acids {% cite Dietz2009, Ke2009, Woo2011 -f rna_ising_model %}.

One can certainly consider a more detailed thermodynamic model with more geometric constraints, eventually achieving a physically and chemically realistic, if not accurate model of nucleic acids. The possible constraints are the strength of chemical bonds, volume exclusion, i.e. two atoms or bases cannot occupy the same position in space, and persistent length {% cite Hagerman1988 -f rna_ising_model %} i.e. two bases sufficiently far apart can be oriented independently in space. A model then can be simple as deriving how exactly bases can pair {% cite Kimchi2019 -f rna_ising_model %}, to an abstract 3D model of nucleic acids such as oxDNA {% cite Sengar2021 -f rna_ising_model %}, to molecular dynamics simulation using force field {% cite Liebl2023 -f rna_ising_model %}, but they all have speed as a bottleneck.

The problem is then how to find a good trade-off between performance and accuracy, preferably augmenting existing thermodynamic models with some notions of geometric constraints whilst maintaining tractability of MFE problem. The goal of this work is thus two-fold.

Firstly, we augment two energy models by allowing nucleotides to have angles with respect to the backbone, and to simplify this notion, each base is oriented either up or down, analogous to the spin in Ising model. Note that this is the only similarity with Ising model in statistical physics, and in particular, we have neither conditions on the spins on the boundary, nor notion of magnetisation. The energy models are base-pair counting by Nussinov {% cite Nussinov1980 -f rna_ising_model %}, chosen for its simplicity, and base-pair stacking by Lyngsø {% cite Lyngsø2004 -f rna_ising_model %}, chosen for its featuring of stack, which is an essential feature of nearest neighbour model {% cite Mathews1999 -f rna_ising_model %}.

Secondly, we study the computational tractability of problems associated with these models, namely the problem of calculating the energy of a secondary structure, and that of finding the minimum free energy. In particular, we show that in almost all cases, these models are computationally hard, suggesting the same if more geometric constraints are ever considered.

### 1.1. Related works

As we already noted, whilst we have learned a lot about geometry of nucleic acids {% cite Hartmann1996, Keating2011 -f rna_ising_model %}, little effort has been made to take this into account in _de novo_ thermodynamic prediction of structures. In particular, no work has been done on Ising-like model for nucleic acids, although similar models have been considered for proteins {% cite Istrail2009 -f rna_ising_model %}, where $\NP$-hardness has been shown for various models, whilst $\APX$- and $\sP$-hardness have not been considered. However, the majority of these protein folding models are lattice-based, e.g. HP model, which differ significantly from nucleic acid secondary structure model, thus rendering these results inapplicable.

Beyond Ising-like models, there have been some attempts to augment current energy models with some geometric constraints. For example, in nature, RNA folds as it is transcribed, a process called co-transcriptional folding, which has been shown to play an important role in its structure formation {% cite Lai2013 -f rna_ising_model %}. To this end, Proctor and Meyer augment Turner model to incorporate this effect {% cite Proctor2013 -f rna_ising_model %}. In a different direction, Deigan et al. {% cite Deigan2009 -f rna_ising_model %} augments Turner model with affinity of each nucleotide as determined experimentally. One should note, however, that these attempts are to refine Turner model with additional parameters, and none of them significantly change the model. In particular, as one starts out with Zuker's algorithm, tractability is not a question.

### 1.2. Organisation
The next section gives quick reminder of essential notions related to the RNA folding problem. Then, Section 3 introduces our Ising model and derives energy functions analogous to what was considered by Nussinov and Jacobson {% cite Nussinov1980 -f rna_ising_model %}, and by Lyngsø {% cite Lyngsø2004 -f rna_ising_model %}. Section 4 introduces two relevant algorithmic problems and states the main computational complexity results. Subsequent sections give proofs or outlines thereof, before Section 7 ends with some remarks.

## 2. Preliminaries

### 2.1. Nucleotide alphabet
We consider in this work a general notion of nucleotides (also called bases), where each of them is a character in an alphabet $\Sigma = \\{s_i, t_i \mid 1 \leq i \leq N\\}$, and the only allowed pairings are between $s_i$ and $t_i$.

For instance, Watson-Crick pairing correspond to $N = 4$, $\Sigma = \\{\A, \T, \C, \G\\}$, and $s_1 = \A$, $t_1 = \T$, $s_2 = \C$ and $t_2 = \G$. Note that in case of RNA, we also have wobble pairs, i.e. $\U$ can pair with either $\G$ or $\A$, but we shall omit them in this work for simplicity.

### 2.2. Secondary structure
A nucleic acid can be encoded as a string $q = \overline{q\_1, q\_2, \cdots, q\_n}\in \Sigma^*$, which may admit a secondary structure $S = \\{(\ell\_i, r\_i)\\}\_{i = 1}^s \in \\{1, 2, \cdots, n\\}^2$, who then must satisfy the following conditions in order to be considered as _admissible_.
- For all $1 \leq i \leq s$, we have that $\ell_i < r_i$. Moreover, $(q_{\ell_i}, q_{r_i}) = (s_j, t_j)$ or $(q_{\ell_i}, q_{r_i}) = (t_j, s_j)$ for some $1 \leq j \leq N$.
- A nucleotide can pair with at most one other, i.e. for all $1 \leq i, j \leq s$, if $\ell_i = \ell_j$ or $r_i = r_j$, then $i = j$.
For future convenience, we denote by $P = \\{i \mid \exists j, (i, j) \in S \lor (j, i) \in S\\}$ the positions of paired nucleotides.

### 2.3. Base-pair counting model and base-pair stacking model

Base-pair counting model, hereinafter abbreviated as BPC, was the first to be considered in RNA folding problem, in the paper by Nussinov and Jacobson {% cite Nussinov1980 -f rna_ising_model %}, where $\Delta'_1(S) = -\|S\|$.

On the other hand, base-pair stacking model, hereinafter abbreviated as BPS, was introduced by Lyngsø {% cite Lyngsø2004 -f rna_ising_model %}, where he showed that finding MFE in this model is $\NP$-hard. Here, the energy function is defined as $\Delta'_2(S) = -\|\\{(i, j) \in S \mid (i+1, j-1) \in S\\}\|$. Intuitively speaking, two bonds of the form $(i, j)$ and $(i+1, j-1)$ form what is called a stack, and this model aims at maximising the number of stacks.

Stack is an essential feature of celebrated Turner model {% cite Mathews1999 -f rna_ising_model %}, thus BPS model can be considered not only as a simplification but also a precursor thereof, and any negative computational complexity result for it strongly suggests the same for Turner model, which is our motivation.

## 3. Ising model for nucleic acids

Analogous to Ising model in physics, we first need a notion of spin. Given a sequence $q$ of length $n$, the nucleotide at position $i$ is assigned a discrete variable $\sigma_i \in \\{1, -1\\}$. The collection of such assignment is called a _spin configuration_ $\sigma = (\sigma_i)_{i = 1}^n \in \\\{1, -1\\\}^n$.

Next, a bond should be symmetric, meaning that if there is a pair $(\ell_i, r_i) \in S$, then it should make no differences whether we have $q_{\ell_i} = s_j$ and $q_{r_i} = t_j$, or $q_{\ell_i} = t_j$ and $q_{r_i} = s_j$. Thus like in the original Ising model, if a bond is possible, we favour the bond between nucleotides of opposite spins. To put it rigorously, we assign an energy $\alpha \sigma_i \sigma_j$ to the pair $(i, j) \in S$, where $\alpha > 0$ is an absolute constant.

On the other hand, we also prefer that nucleotides in proximity _and paired_ with some nucleotides not necessary amongst themselves. In other words, if the bases at position $i$ and $i+1$ are paired with those of position $j$ and $k$ respectively, i.e. $(i, j), (i+1, k) \in S$ for some $j$ and $k$, then we prefer that $\sigma_i = \sigma_{i+1}$, and thus we assign an energy $-\sigma_i \sigma_{i+1}$.

The motivation for this is less clear (aside from making the MFE problems non-trivial), but there is some from chemistry. In reality, the bond between two nucleotides is the hydrogen bond between the corresponding nitrogenous bases, which are attached to phosphate group. These groups are linked by ester bond, together forming the backbone of nucleic acid. As such, two consecutive paired nucleotides $i$ and $i+1$ cannot form arbitrary dihedral angles {% cite Hartmann1996 -f rna_ising_model %}. In the special cases where $i$ is paired with $j$ whilst $i+1$ is paired with $j-1$, it might be a part of a helix structure, thus have their angles fixed.

On the other hand, we do not consider $\sigma_i\sigma_{i+1}$ for all $i$, because the interaction above is, to some degree, local. In particular, if we see nucleotides as vectors from the phosphate group to the base, then two nucleotides of distance larger than persistent length have their directions uncorrelated {% cite Hagerman1988 -f rna_ising_model %}. Thus we choose the aforementioned definition to abstractise this effect.

Finally, we introduce the following energy function as sum of the aforementioned contributions 
\\[
    \Delta_1(S, \sigma) = \alpha \sum_{(i, j) \in S} \sigma_i \sigma_j - \sum_{i, i+1 \in P} \sigma_i \sigma_{i+1}.
\\]
One can see it, and in particular the first term, as an analogy to $\Delta'\_1$ in BPC model, by writing $\Delta'\_1 (S) = \sum_{(i, j) \in S} (-1)$. Thus we call this Ising BPC model. It is now clear that $\alpha$ plays the role of ratio between the strength of two interaction, the choice of which we delay until later, as there is no reason to suggest any particular value.

In this fashion, it is straightforward to define an analogy, which we call Ising BPS model, to $\Delta'\_2$ in BPS model as follow
\\[
    \Delta\_2(S, \sigma) = \alpha \sum\_{(i, j) \in S \mid (i+1, j-1) \in S} \sigma\_i \sigma\_j - \sum\_{i, i+1 \in P} \sigma\_i \sigma\_{i+1}.
\\]

## 4. Statements of results

Given our two energy functions, it is straightforward to define the problem of finding minimum free energy.

**Problem $\textsf{MFE}\_\{\Delta\}$ ($\Delta \in \{\Delta\_1, \Delta\_2\}$).** Given $\alpha$ and a nucleotide sequence $q$ as input, output $\min_{S, \sigma} \Delta(S, \sigma)$.

Another interesting problem is given a secondary structure $S$, evaluating the free energy of $S$. In the classical setting, the energy is entirely determined by $S$, whose evaluation is done in polynomial time. However, in our new setting, the energy also depends on the choice of spin configuration $\sigma$, and evaluating the energy of $S$ amounts to an optimisation problem over all possible choices of $\sigma$, which a priori is non-trivial. We thus introduce it as a separate problem.

**Problem $\textsf{SPIN}\_\{\Delta\}$ ($\Delta \in \{\Delta\_1, \Delta\_2\}$).** Given $\alpha$, a nucleotide sequence $q$, and a secondary structure $S$ as input, output $\min_{S, \sigma} \Delta(S, \sigma)$ if $S \in \mathcal{S}$, and $\infty$ otherwise.

Our results are summarise in the following theorems.

**Theorem 1.** $\textsf{SPIN}\_{\Delta\_1}$ and $\textsf{SPIN}\_{\Delta\_2}$ are $\sP$- and $\APX$-hard.

**Theorem 2.** $\textsf{MFE}_{\Delta_1}$ is $\APX$-hard and in $\|\Sigma\|-W[1]$. When $\Sigma$ is unbounded, it is $\sP$-hard.

**Theorem 3.** $\textsf{MFE}_{\Delta_2}$ is $\APX$-hard, $\NP$-hard when $\Sigma$ is bounded, and $\sP$-hard when it is unbounded.

### References

{% bibliography --cited -f rna_ising_model %}