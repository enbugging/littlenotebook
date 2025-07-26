---
layout: post
title: "Robustness problem in RNA folding: telescoping method"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

$$
\newcommand{\vertex}{\text{vert}}
$$

This is mainly taken from Chapter 6 of my Bachelor thesis at École Polytechnique under supervision of Prof. Sarah J. Berkemer and Prof. Yann Ponty, titled ["Polynomial-time parametric optimisation"](https://arxiv.org/abs/2304.14962)

## 1. Introduction

Along side with the problem of parameter learnability (and more generally relative position problem, see my other posts [(Nguyen, 2024)](relative_position_problem_motivation.html)), one can study how robust a parameter set is. To briefly see why this is relevant, it is believed that the environment surrounding an RNA can alter how the nucleotides pair with each other, effectively modifying the parameter of energy model. In that sense, robustness measures how stable an RNA structure is in changes of the surrounding: a stable RNA is more resistant to mutations and undesirable functions or defects thereof, which is a desirable property and supported by evolution.

For notations, readers are referred to my introduction on the topic [(Nguyen, 2024)](intro_to_RNA_polytopes.html). Formally speaking, let $$P$$ be a polytope, $$p \in \mathbb{R}^{d+1}$$, we define the robustness of $$p$$ (with respect to $$P$$), denoted $$\theta_P(p)$$, to be the supremum of angle $$\theta$$ such that for all $$p' \in \mathbb{R}$$ satisfying $$\cos^{-1} ((p, p')) = \frac{\langle p, p'\rangle}{\| p \|  \| p' \|} \leq \theta$$ and for all $$x \in P$$, if $$\langle x, p \rangle = h_P(p)$$ then $$\langle x, p' \rangle = h_P(p')$$, or equivalently, if $$\langle x, p \rangle < h_P(p)$$ then $$\langle x, p' \rangle < h_P(p')$$.

Since if $$\langle x, p \rangle = h_P(p)$$ then for any $$\lambda \geq 0$$, one has $$\langle x, \lambda p \rangle = h_P(\lambda p)$$, it is thus customary to consider $$p \in \mathbb{S}^d$$, and the notion of robustness coincides with the minimum distance from $$p$$ to another normal spherical polytope $$S \in \mathcal{S}(P)$$ not containing $$p$$, which since there are only finitely many polytopes all of which are closed, the supremum is attained.

Then, consider the RNA polytope $$\mathcal{P}(q)$$ of a sequence $$q$$, and we know that a secondary structure $$s$$ is optimal for $$q$$ with respect to a parameter set $$p$$, one can define the robustness of $$p$$ as $$\theta_{\mathcal{P}(q)} (p)$$. The challenge is how to compute this efficiently, e.g. without resorting to compute the entire RNA polytope, in view of many computational hardness involved (see review on the topic [(Nguyen, 2024)](relative_position_problem_motivation.html)).

Here, I present a method whereby, similarly to modification of the original dynamic programming algorithm from tropical algebra to polytope algebra as presented by Pachter and Sturmfels {% cite Pachter2004 Pachter2005 -f rna_polytopes %}, one can modify the algorithm further to retrieve only a certain part of the polytope $$\mathcal{P}(q)$$. The key idea is that for any two polytopes $$A$$ and $$B$$, a given quadrant of $$A \oplus B$$ or of $$A \otimes B$$ is determined only by the corresponding quadrants of $$A$$ and $$B$$, as demonstrated in these figures below.

![](/littlenotebook/assets/img/rna_polytopes/A.png#center)
<div align="center">Figure 1a. $F_D(A)$.</div>
![](/littlenotebook/assets/img/rna_polytopes/B.png#center)
<div align="center">Figure 1b. $F_D(B)$</div>
![](/littlenotebook/assets/img/rna_polytopes/AoplusB_partial.png#center)
<div align="center">Figure 1c. $F_D(A \oplus B)$</div>
![](/littlenotebook/assets/img/rna_polytopes/AotimesB_partial.png#center)
<div align="center">Figure 1d. $F_D(A \otimes B)$</div>

This is stated formally as follow.

**Theorem 1.** Let $$D \subseteq \mathbb{R}^{d+1}$$, $$A$$ and $$B$$ be two polytopes in $$\mathbb{R}^{d+1}$$. For a given polytope $$P \in \mathbb{R}^{d+1}$$, we define 

$$
	F_D(P) = \{x \mid \exists d \in D, \langle x, d\rangle = h_P(d)\}.
$$

Then, one has

$$
    F_D(A \oplus B) \subseteq F_D(A) \oplus F_D(B), 
$$

and

$$
    F_D(A \otimes B) \subseteq F_D(A) \otimes F_D(B).
$$


_Proof._ Let $$d \in D$$, one has 

$$
\max_{x \in A \otimes B} \langle x, d \rangle = \max_{(a, b) \in A \times B} \langle a + b, d \rangle = \max_{(a, b) \in A \times B} \langle a, d \rangle + \langle b, d \rangle = \max_{a \in A} \langle a, d \rangle + \max_{b \in B} \langle b, d \rangle
$$

which implies that for $$x = a + b \in A \otimes B$$ for some $$a \in A$$ and $$b \in B$$, then if $$x \in F_D(A \otimes B)$$ then $$a \in F_D(A)$$ and $$b \in F_D(B)$$, implying $$x \in F_D(A) \otimes F_D(B)$$, as desired.

On the other hand, by definition, for any $$x \in A \oplus B$$, there exist $$\lambda, \mu \geq 0$$ and $$a \in A$$, $$b \in B$$, such that $$\lambda + \mu = 1$$ and $$x = \lambda a + \mu b$$. Now suppose $$\lambda > 0$$, and $$\langle x, d \rangle = h_{A \oplus B}(d)$$, if there exists $$a' \in A$$ such that $$\langle a', d \rangle > \langle a, d \rangle$$, then for $$x' = \lambda a' + \mu b \in A \oplus B$$, one has

$$
\langle x, d \rangle = \lambda \langle a, d \rangle + \mu \langle b, d \rangle < \lambda \langle a, d \rangle + \mu \langle b, d \rangle = \langle x', d \rangle, 
$$

a contradiction, thus $$a \in F_D(A)$$. Similarly, if $$\mu > 0$$ then $$b \in F_D(B)$$, which together implies $$x \in F_D(A) \oplus F_D(B)$$. <span style="float:right;">$$\square$$</span>

This results allows one to specify beforehand a quadrant $$D$$ of interest, and as they carry out the modified dynamic programming algorithm according to Pachter and Sturmfels {% cite Pachter2004 Pachter2005 -f rna_polytopes %}, only retain relevant parts of the polytope, thus significantly reduce the complexity of the polytopes and speed up the computation.

## 2. Toward robustness

Now suppose we know a pair $$(q, s)$$ is learnable for a parameter vector $$p$$. To determine the robustness of $$p$$, we shall consider $$D$$ to be of some special forms: in particular, we restrict $$D$$ to $$\mathbb{S}^d$$, and consider $$D$$ to be a closed ball $$B_{\mathbb{S}^d}(p, \theta) = \{x \in \mathbb{S}^d \mid d_{\mathbb{S}^d}(x, p) \leq \theta\}$$ for some $$\theta > 0$$. Let $$F$$ be a polytope of $$\mathcal{S}(\mathcal{P}(p))$$ of minimum dimension that contains $$p$$, then the robustness of $$p$$ is the minimum $$\theta$$ such that $$B_{\mathbb{S}^d}(p, \theta)$$ intersects another face $$F'$$ of $$\mathcal{S}(\mathcal{P}(p))$$, which one can limit to only the adjacent faces of $$F$$. Or, equivalently, it is the infimum of $$\theta$$ such that $$F$$ contains $$B_{\mathbb{S}^d}(p, \theta)$$.

Note that we know a priori that $$0 \leq \theta_P(p) \leq \pi$$, thus we can compute $$\theta_P(p)$$ via binary search. But it is also possible to compute the robustness with only one execution of the dynamic programming algorithm, and no binary search is necessary. Similar to the above, we shall modify the dynamic programming algorithm once more, as follow.

One crucial observation is that $$D$$ needs not be fixed throughout the dynamic programming process, but in fact can be decreased, i.e. the new set $$D'$$ for the next iteration of the dynamic programming needs only be a subset of $$D$$, since by definition, if $$D' \subseteq D$$, then for any polytope $$P$$, $$F_{D'}(P) \subseteq F_D(P)$$. This gives us the idea to set initially $$\theta = \pi$$ and $$D = B_{\mathbb{S}^d}(p, \theta) = \mathbb{S}^d$$, and as we carry out the dynamic programming with $$F_D$$, we may decrease $$\theta$$ until $$\theta$$ equals $$\theta_{\mathcal{P}(q)}(p)$$. 

There remains one last question: can we decrease $$\theta$$? Or, in other words, given $$F_{B_{\mathbb{S}^d}(p, \theta_A(p))}(A)$$ and $$F_{B_{\mathbb{S}^d}(p, \theta_B(p))}(B)$$, can we construct $$F_{B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))}(A \oplus B)$$ and $$F_{B_{\mathbb{S}^d}(p, \theta_{A \otimes B}(p))}(A \otimes B)$$? The answer is yes, and the following theorem provides rigorous reasoning.

**Theorem 2.** Let $$A, B \subseteq \mathbb{R}^{d+1}$$ be two polytopes, and $$p \in \mathbb{R}^{d+1}$$. With notations defined as in Theorem 1, one has 

$$
    F_{B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))}(A \oplus B) \subseteq F_{B_{\mathbb{S}^d}(p, \theta_A(p))}(A) \oplus F_{B_{\mathbb{S}^d}(p, \theta_B(p))}(B)
$$

and

$$
    F_{B_{\mathbb{S}^d}(p, \theta_{A \otimes B}(p))}(A \otimes B) \subseteq F_{B_{\mathbb{S}^d}(p, \theta_A(p))}(A) \otimes F_{B_{\mathbb{S}^d}(p, \theta_B(p))}(B).
$$

_Proof._ We recall Proposition 7.12 in Ziegler's _Lectures on Polytopes_ {% cite Ziegler1995 -f rna_polytopes %} stating that $$\mathcal{N}(A \otimes B) = \mathcal{N}(A) \wedge \mathcal{N}(B)$$, where $$\cdot \wedge \cdot$$ denotes the common refinement. Thus, it follows immediately by definition that $$\theta_{A \otimes B} (p) = \min (\theta_A(p), \theta_B(p))$$, whence together with Theorem 1, one deduces
<div align="center">
$$
\begin{split}
    F_{B_{\mathbb{S}^d}(p, \theta_{A \otimes B}(p))}(A \otimes B)
    & \subseteq F_{B_{\mathbb{S}^d}(p, \theta_{A \otimes B}(p))}(A) \otimes F_{B_{\mathbb{S}^d}(p, \theta_{A \otimes B}(p))}(B) \\
    & \subseteq F_{B_{\mathbb{S}^d}(p, \theta_A(p))}(A) \otimes F_{B_{\mathbb{S}^d}(p, \theta_B(p))}(B).
\end{split}
$$
</div>
	
Now let $$x \in F_{B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))}(A \oplus B)$$, and let $$q \in B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))$$ such that $$\langle x, q \rangle = h_{A \oplus B}(q)$$. We notice that for any $$r \in \mathbb{R}^{d+1}$$, one has $$h_{A \oplus B}(r) = \max(h_A(r), h_B(r))$$. Suppose without loss of generality that $$x \in A$$, then it follows that $$h_A(q) \geq h_B(q)$$, otherwise one would have $$\max(h_A(q), h_B(q)) = \langle x, q \rangle \leq h_A(q) < h_B(q)$$,
a contradiction.

Suppose by contradiction that $$h_A(p) < h_B(p)$$. Let $$y = \sigma_A(q)$$ then $$\langle y, q \rangle = h_{A \oplus B}(q)$$. Since $$q \in B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))$$, one has

$$
    h_A(p) \geq \langle y, p \rangle = h_{A \oplus B}(p) = h_B(p) > h_A(p),   
$$

a contradiction. Hence, $$h_A(p) \geq h_B(p)$$.

Then for any $$z \in A$$, if $$\langle z, p \rangle = h_A(p)$$, one has $$\langle z, p \rangle = h_{A \oplus B}(p)$$, and by assumption $$q \in B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))$$, it follows that $$\langle z, q \rangle = h_{A \oplus B}(q) = h_A(q)$$. Or, equivalently, $$q \in B_{\mathbb{S}^d}(p, \theta_A(p))$$, which implies $$x \in F_{B_{\mathbb{S}^d}(p, \theta_A(p))}(A)$$. Then, it follows that

$$
    F_{B_{\mathbb{S}^d}(p, \theta_{A \oplus B}(p))}(A \oplus B) \subseteq F_{B_{\mathbb{S}^d}(p, \theta_A(p))}(A) \oplus F_{B_{\mathbb{S}^d}(p, \theta_B(p))}(B), 
$$

as desired. <span style="float:right;">$$\square$$</span>

## 3. Computing robustness given a polytope

For a given polytope $$P$$ in the course of executing the dynamic programming, Theorem 2 allows one to consider only $$F_{B_{\mathbb{S}^d}(p, \theta_P(p))}(P)$$, but does not explicitly say how to compute it. Whilst it is possible via constructing the normal fan, we remind that translating from $$\mathcal{V}$$- to $$\mathcal{H}$$-representation introduces a bottleneck and is responsible for the complexity of the naive approach. For this reason, we introduce the following theorem which characterises $$D = B_{\mathbb{S}^d}(p, \theta)$$ when $$\theta = \theta_{\mathcal{P}(q)}(p)$$. 

First, we have a technical lemma.

**Lemma 3.** Let $$P \subseteq \mathbb{S}^d$$ be a polytope, $$x \in \mathring{P}$$, $$y \not\in P$$. Then $$[x, y]_{\mathbb{S}^d}$$ intersects $$\partial P$$.

_Proof._ Consider the map $$f : [0, 1] \mapsto [x, y]_{\mathbb{S}^d}$$, where $$f(t) = \frac{tx + (1-t)y}{\|tx + (1-t)y\|}$$. We denote $$T = \{t \mid f(t) \in P\}$$, and $$t^* = \sup T$$. A priori we know that $$0 < t^* < 1$$. Let $$z = f(t^*)$$, if $$z \not\in P$$, then there exists $$\theta > 0$$ such that $$B_{\mathbb{S}^d}(z, \theta) \cap P = \emptyset$$ and $$y \not\in B_{\mathbb{S}^d}(z, \theta)$$. Let $$\varepsilon = \frac{\theta}{d_{\mathbb{S}^d}(x, y)}$$ which is well-defined since $$x \neq y$$, one has $$[\theta - \varepsilon, \theta + \varepsilon] \cap T = \emptyset$$, a contradiction. Hence $$z \in P$$.

Now in similar fashion, if $$z \in \mathring{P}$$, then there exists $$\theta > 0$$ such that $$B_{\mathbb{S}^d}(z, \theta) \subsetneq P$$ and $$x \not\in B_{\mathbb{S}^d}(z, \theta)$$. Let $$\varepsilon = \frac{\theta}{d_{\mathbb{S}^d}(x, y)}$$, one has $$[\theta - \varepsilon, \theta + \varepsilon] \subsetneq T$$, a contradiction. Thus $$z \in \partial P$$, as desired. <span style="float:right;">$$\square$$</span>

**Theorem 4.** Let $$P \subseteq \mathbb{S}^d$$ be a polytope, $$p \in \mathring{P}$$, and $$p^* \in \partial P$$ such that 

$$
    d_{\mathbb{S}^d}(p, p^*) = d_{\mathbb{S}^d}(p, \partial P) = \min_{p' \in \partial P} d_{\mathbb{S}^d}(p, p').
$$

Then $$p^*$$ lies in the interior of a facet of $$P$$, i.e. no faces of $$P$$ with lower dimension contain $$p^*$$.

_Proof._ Suppose the contrary, there exist two different facets of $$P$$ containing $$p^*$$. Since $$p \neq p^*$$, the segment $$[p, p^*]$$ cannot be perpendicular to both of the facets, hence there must exist a facet $$F$$ of $$P$$ containing $$p^*$$ such that $$[p, p^*]$$ is not perpendicular to $$F$$. Now consider the projection $$p'$$ of $$p$$ onto $$F$$. If $$p' \in F$$, then since 

$$
    d_{\mathbb{S}^d}(p, \partial P) = d_{\mathbb{S}^d}(p, p^*) >  d_{\mathbb{S}^d}(p, p') \geq d_{\mathbb{S}^d}(p, \partial P),
$$

we have a contradiction. Otherwise, by Lemma 3, $$[p, p']$$ intersects $$\partial P$$ at a point $$p''$$, which necessarily differs from $$p'$$ and thus implies that 

$$
    d_{\mathbb{S}^d}(p, \partial P) = d_{\mathbb{S}^d}(p, p^*) > d_{\mathbb{S}^d}(p, p') \geq d_{\mathbb{S}^d}(p, p'') > d_{\mathbb{S}^d}(p, \partial P),
$$

a contradiction, as desired. <span style="float:right;">$$\square$$</span>

In simple terms, if the robustness of $$p$$ is positive, one knows that $$p$$ lies in the interior of some normal spherical polytope $$S \in \mathcal{S}(P)$$, which corresponds to one unique vertex $$x \in P$$ such that $$\langle x, p \rangle = h_P(p)$$, i.e. $$x = \sigma_P(p)$$. Theorem 4 tells us that the robustness $$\theta_P(p)$$ is obtained by an edge between $$x$$ and a neighbour vertex $$y$$ of $$P$$, or more precisely, if we denote $$H = \{z \in P | \langle z, p \rangle = 0\}$$ to be the hyperplane admitting $$p$$ as its normal vector, then

$$
    \theta_P(p) = \min_{y \in \vertex{P}\setminus\{x\}} \angle(y - x, H)= \min_{y \in \vertex{P}\setminus\{x\}} \sin^{-1} \frac{(h_P(p) - \langle y, p \rangle)\| p \|}{\| x - y\|}.
$$

Thus, one only needs to find _all_ points $$y \in \vertex{P}\setminus\{x\}$$ that minimise $$\angle(y - x, H)$$ and compute $$\theta_P(p)$$ without explicitly constructing $$\mathcal{N}(P)$$. This allows us to restrict ourselves to $$\mathcal{V}$$-representation.

### References 
{% bibliography --cited -f rna_polytopes %}