---
layout: post
title: "Private feature learning for single-index models (part 2)"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

This is part 2/3 of my M1 internship per requirement from ENS, done at Institute of Science and Technology Austria, under the supervision of Simone Bombari and Prof. Marco Mondelli, to be published.

## 3. First layer.

In this section, we only use the first dataset $$X^W$$ thus by abuse of notation, we write $$X$$, $$x_j$$, and $$y_j$$ to mean $$X^W$$, $$x^W_j$$, and $$y^W_j$$. First we state the main result of this section.

**Theorem 2.** After the first phase of Algorithm 1 (line 8), $$W_1$$ is $$(\varepsilon, \delta)$$-differentially private. Moreover, there exists a choice of second layer $$\overline{a}$$ such that for $$f^W(x) = \overline{a}^\intercal \sigma(W_1^\intercal x + b)$$, we have

$$
\newcommand{\P}{\mathbb{P}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\E}{\mathbb{E}}
\newcommand{\N}{\mathcal{N}}
\newcommand{\H}{\mathcal{H}}
\newcommand{\L}{\mathcal{L}}
\newcommand{\RR}{\text{RR}}
\newcommand{\he}{\text{He}}
\newcommand{\h}{\text{H}}
\newcommand{\sech}{\text{sech}}
\newcommand{\csch}{\text{csch}}
\newcommand{\id}{\text{Id}}
$$

$$
\begin{equation}
    \E_x \left[f^W(x) - f_*(x)\right]^2 + \lambda \|\overline{a}\|^2 = O_\P\left(\frac{1}{\log d}\right).
\end{equation}
$$

### 3.1. Rank-$$1$$ approximation of gradient.

Let $$f^0(x) = a_0^\intercal \sigma(W_0^\intercal x + b)$$ be the initial model. The gradient for $$w_i$$ is given by

$$
\begin{equation}
\begin{aligned}
    \nabla_{w_i} \L\left(f^0\right) & = \frac{a_{0, i}}{n} \sum_{j = 1}^n \left(f^0(x_j) - f_*(x_j)\right) \cdot \sigma'\left(\langle w^0_i, x_j \rangle + b_j\right)x_j \\
    & = a_{0, i} \left(\frac{1}{n}\sum_{j = 1}^n f^0(x_j) \sigma'\left(\langle w^0_i, x_j \rangle + b_j\right)x_j \right. \\
    & \phantom{=} \left. - \frac{1}{n}\sum_{k = 1}^n f_*(x_k) \cdot \sigma'\left(\langle w^0_i, x_k \rangle + b_k\right)x_j\right)
\end{aligned}
\end{equation}
$$

There are two terms, the former of which is negligible, and the latter can be approximated by its population counterpart. In particular, let

$$g(w^0_i) = \E_x\left[f_*(x) \sigma'\left(\langle w^0_i, x \rangle + b_i\right) x\right]$$

and

$$\hat{g}(w^0_i) = \frac{1}{n}\sum_{k = 1}^n f_*(x_k) \cdot \sigma'\left(\langle w^0_i, x_k \rangle + b_k\right)x_j,$$

we have the following lemma.

**Lemma 3.** We have

$$
\begin{equation}
\begin{gathered}
    \left\|\frac{1}{n}\sum_{j = 1}^n f^0(x_j) \sigma'\left(\langle w^0_i, x_j \rangle + b_j\right)\right\| \lesssim_\P \sqrt{\frac{d}{n}} \log d\\
    \left\|g(w^0_i) - \hat{g}(w^0_i)\right\| \lesssim_\P \sqrt{\frac{d}{n}} \log^{q+1} d.
\end{gathered}
\end{equation}
$$

_Proof._ Note that $$\sigma' = \sech^2 \leq 1$$, so

$$
\begin{equation}
\begin{gathered}
    \left\|\frac{1}{n}\sum_{j = 1}^n f^0(x_j) \sigma'\left(\langle w^0_i, x_j \rangle + b_j\right)\right\| \leq \left\| \frac{1}{n} \sum_{j = 1}^n f^0(x_j) x_j\right\| \\
    = \left\| \frac{1}{n} \sum_{j = 1}^n a_0^\intercal \sigma(W_0^\intercal x_j + b) x_j\right\|.
\end{gathered}
\end{equation}
$$

Let $$z_j = a_0^\intercal \sigma(W_0^\intercal x_j + b)$$. Since $\|\sigma\| \leq 1$, one has $\|a_{0, k} \sigma(\langle w^0_k, x_j \rangle + b)\| \leq \frac{1}{\sqrt{p}}$, implying $\|z_j\| \leq 1$, and so $$\|z_j x_j\| \leq_\P \sqrt{d}$$. Then, let $$v = \max\left(\|\sum_{j = 1}^n \E\left[z_j^2 x_j x_j^\intercal\right]\|, \|\sum_{j = 1}^n \E\left[z_j^2 x_j^\intercal x_j\right]\|\right)$$, we can bound

$$
\begin{equation}
    v \leq \sum_{j = 1}^n \E[z_j^2 \|x_j\|^2] \leq nd.
\end{equation}
$$

By matrix Bernstein's inequality {% cite Tropp2015 -f private_feature_learning_for_single_index_models.bib %}, we have

$$
\begin{equation}
    \P\left(\left\| \frac{1}{n} \sum_{j = 1}^n z_j x_j\right\| \geq t\right) \leq (d+1) \exp\left(\frac{-\frac{1}{2}n^2 t^2}{nd + \frac{1}{3}t \sqrt{d}}\right).
\end{equation}
$$

Letting $$t = \sqrt{\frac{d}{n}} \log d$$, the proof for the first inequality is complete.

For the second inequality, let $$u_j = x_j f_*(x_j) \sigma\left(\langle w^0_i, x_j \rangle  + b_i\right)$$, $$\overline{u_j} = u_j - \E[u_j]$$, and $$U = \sum_{j = 1}^n \overline{u_j}$$. Since $\|\sigma'\| \leq 1$, we have

$$
\begin{equation}
    \|\overline{u_j}\|  \leq \|x_j\| |\sigma_*\left(\langle \mu, x_j \rangle\right)| + \E\left[\|x_j\| |\sigma_*\left(\langle \mu, x_j \rangle\right)|\right]
\end{equation}
$$

Since $$x_j \sim \N(0, \id)$$, $$\langle \mu, x \rangle \sim \N(0, 1)$$, so by concentration of Gaussian variable and of Gaussian vector's norm {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %}, for all $$j$$ simultaneously, we have $$\|x_j\| \leq_\P \sqrt{d}$$ and $\|\langle \mu, x_j \rangle\| \leq_\P \log d$, implying $\|\sigma_*(\langle \mu, x_j \rangle)\| \lesssim_\P \log^q d$, and therefore, $$\|\overline{u_j}\| \lesssim_\P \sqrt{d} \log^q d$$. Arguing as done with $$z_j x_j$$ gives the second inequality. <span style="float:right;">$$\square$$</span>

By Stein's lemma, we write

$$
\begin{equation}
	g(w^0_i) = \mu \cdot \E_x \left[f'_*(x) \sigma' \left(\langle w^0_i, x \rangle + b\right)\right] + w^0_i \cdot \E_x \left[f_*(x) \sigma''\left(\langle w^0_i, x\rangle\right)\right].
\end{equation}
$$

Using Hermite coefficients, we have

$$
\begin{equation}
	\E_x \left[f_*(x) \sigma'' \left(\langle w^0_i, x \rangle + b\right)\right] = \sum_{j = 0}^\infty (j+1)(j+2) \alpha_{j+2}^{b_i} \alpha^*_j \langle \mu, w^0_i \rangle^j.
\end{equation}
$$

where Assumption 1 gives $$\alpha^*_0 = 0$$. On the other hand, initialisation gives $$w^0_i \sim \N\left(0, \frac{1}{d} \id\right)$$, so $$\langle \mu, w^0_i \rangle \sim \N\left(0, \frac{1}{d}\right)$$, and by tail bound of Gaussian variables {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %}, we have $\|\langle \mu, w^0_i \rangle\| \lesssim_\P \frac{\log d}{\sqrt{d}}$. Together this implies

$$
\begin{equation}
	\E_x \left[f_*(x) \sigma'' \left(\langle w^0_i, x \rangle + b\right)\right] = O_\P \left(\frac{\log d}{\sqrt{d}}\right).
\end{equation}
$$

Similarly, we write

$$
\begin{equation}
	\label{eq:first layer - second term of Stein lemma}
	\E_x \left[f'_*(x) \sigma' \left(\langle w^0_i, x \rangle + b\right)\right] = \sum_{j = 0}^\infty (j+1)^2 \alpha_{j+1}^{b_i} \alpha^*_{j+1} \langle \mu, w^0_i \rangle^j.
\end{equation}
$$

Now to handle $$\alpha^{b_i}_1$$, we have the following lemma.

**Lemma 4.** For $$\sigma = \tanh$$, we have

$$
\alpha^b_1 = \frac{1}{\sqrt{2\pi}} \int_\R \sigma(x+b)\cdot xe^{-\frac{x^2}{2}} \, dx = \Omega(e^{-2|b|}).
$$

_Proof._ By parity of $$\tanh$$, $$\alpha^b_1 = \alpha^{-b}_1$$, so we can assume without loss of generality that $$b \geq 0$$. By integration by part and using $$\tanh' = \sech^2 = \frac{4}{(e^x+e^{-x})^2}$$, we have

$$
\begin{equation}
    \begin{aligned}
        \alpha^b_1 & = \frac{1}{\sqrt{2\pi}}\int_\R \sigma(x + b)\cdot xe^{-\frac{x^2}{2}}\, dx \\
        & = \frac{1}{\sqrt{2\pi}}\left[\sigma(x+b) e^{-\frac{x^2}{2}} \Big\vert_{-\infty}^\infty - \int_\R \left[ -e^{-\frac{x^2}{2}}\right] \sech^2(x+b) \, dx\right] \\
        & = \frac{1}{\sqrt{2\pi}}\int_\R \frac{4}{(e^{x+b}+e^{-x-b})^2} e^{-\frac{x^2}{2}} \, dx. 
    \end{aligned}
\end{equation}
$$

Now note that for $$x \geq 0 \geq -b$$, we have $$\frac{2}{e^{x+b}+e^{-x-b}} \geq \frac{1}{e^{x+b}}$$, so noticing that the integrand is non-negative, we have the following lower bound.

$$
\begin{equation}
\begin{aligned}
    \alpha^b_1 
    & \geq \frac{1}{\sqrt{2\pi}}\int_0^\infty \frac{4}{(e^{x+b}+e^{-x-b})^2} e^{-\frac{x^2}{2}} \, dx \\
    & \geq \frac{1}{\sqrt{2\pi}}\int_0^\infty \frac{1}{e^{2(x+b)}} e^{-\frac{x^2}{2}} \, dx \\
    & = e^{-2b} \left(\frac{1}{\sqrt{2\pi}}\int_0^\infty e^{-\frac{x^2}{2} - 2x} \, dx\right)
\end{aligned}
\end{equation}
$$

as desired.<span style="float:right;">$$\square$$</span>

By initialisation $$b_i \sim \N(0, 1)$$ and concentration of Gaussian variable {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %}, for $$3\varepsilon_n - \frac{\varepsilon_p}{4} > \varepsilon_b > 0$$ fixed, we have $\|b_i\| \leq_\P \frac{\varepsilon_b}{4} \log d$ for all $$1 \leq i \leq p$$, implying $$\alpha^{b_i}_1 \geq_\P \frac{d^{-\varepsilon_b}{2}} > 0$$. On the other hand, Assumption 1 gives $$\alpha^*_1 \neq 0$$, thus in the right-hand side of $$\eqref{eq:first layer - second term of Stein lemma}$$, the term with $$j = 0$$ dominates the sum, and we have

$$
\begin{equation}
	g(w^0_i) = \alpha^{b_i}_1 \alpha^*_1 \mu + O_\P \left(\frac{\log d}{\sqrt{d}}\right).
\end{equation}
$$

Now note that by concentration of Gaussian's vector's norm {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %}, we have $$\|w_i^0\| \leq_\P 1$$. Due to our choice of hyper-parameters and Lemma 4, we have $$\eta_W \alpha_1^{b_i} \alpha^*_1 \gtrsim_\P \alpha^*_1 d^{\frac{3\varepsilon_n - \varepsilon_b}{2}} > 0$$. Therefore, the rank-$$1$$ term, $$a_{0, i} \alpha_1^{b_i} \alpha^*_1 \mu$$ is much larger than $$w_i^0$$. Combining with Lemma 3, we have 

$$
\begin{equation}
\begin{aligned}
    \nabla_{w_i} \L
    & = -a_{0, i} \left[\alpha_1^{b_i} \alpha^*_1 \mu + O_\P \left(\sqrt{\frac{d}{n}} \log^{q+1} d\right)\right] \\
    \Rightarrow w_i^0 - \eta_W \nabla_{w_i} \L 
    & = w_i^0 + \eta_W a_{0, 1}  \left[\alpha_1^{b_i} \alpha^*_1 \mu + O_\P \left(\sqrt{\frac{d}{n}} \log^{q+1} d\right)\right] \\
    & = \eta_W \alpha_1^{b_i} \alpha^*_1 \mu + \eta_W O_\P \left(\sqrt{\frac{d}{n}} \log^{q+1} d\right)
\end{aligned}
\end{equation}
$$

giving

$$
\begin{equation}
    \label{eq: first layer - rank-1 approximation}
    w_i^1 = \frac{w_i^0 - \eta_W \nabla_{w_i} \L}{\left\| w_i^0 - \eta_W \nabla_{w_i} \L\right\|} = \mu + O_\P \left(\sqrt{\frac{d}{n}} \log^{q+1} d\right).
\end{equation}
$$

### 3.2. Differential privacy of $W_1$.

Now recall the definition of $$\ell_2$$-sensitivity.

**Definition 2.** {% cite Dwork2013 -f private_feature_learning_for_single_index_models.bib %} For an algorithm $$\mathcal{A}$$ returning values in $$\R^d$$, its $$\ell_2$$-sensitivity is given by

$$
\begin{equation}
    \Delta_2(\mathcal{A}) = \sup_{D, D' \text{ adjacent}} \left\| \mathcal{A}(D) - \mathcal{A}(D') \right\|_2.
\end{equation}
$$

This in turn dictates the level of noise one must add to ensure differential privacy.

**Theorem 5.** {% cite Dwork2013 -f private_feature_learning_for_single_index_models.bib %} Let $$\varepsilon \in (0, 1)$$ be arbitrary. For $$c^2 > 2 \log(1.25/\delta)$$, the Gaussian mechanism with parameter $$\sigma \geq c \Delta_2(\mathcal{A})/\varepsilon$$ is $$(\varepsilon, \delta)$$-differentially private.

In our case, $$\eqref{eq: first layer - rank-1 approximation}$$ implies that we have $$\Delta_2 \lesssim_\P \sqrt{\frac{d^{1 + \varepsilon_b}}{n}} \log^{q+1} d < \sqrt{\frac{d^{1 + \varepsilon_b}}{n}} \log^{q+2} d$$ for $$d$$ large enough, thus it suffices to add the noise of level

$$
\begin{equation}
	\sigma_W = \sqrt{\frac{d^{1 + \varepsilon_b}}{n}} \log^{q+2} d \frac{\sqrt{2\log(1.25/\delta)}}{\varepsilon} < \sqrt{\frac{d^{1 + \varepsilon_b}}{n}} \log^{q+3} d
\end{equation}
$$

where the extra $$\log d$$ factor (notice the exponent $$\log^{q+2} d$$) guarantees that Theorem 5 applies, consistent with our choice of hyper-parameters (Assumption 3). Then, standard concentration of Gaussian variables {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %} adds another $$\log d$$ factor, and we have

$$
\begin{equation}
	w_i^1 = \mu + O_\P \left(\sqrt{\frac{d^{1 + \varepsilon_b}}{n}} \log^{q+1} d \left[1 + \log^2 d \frac{\sqrt{\log(1.25/\delta)}}{\varepsilon}\right]\right).
\end{equation}
$$

### 3.3. Universal approximation theorem for $$\tanh$$.
As analogue of universal approximation theorem for ReLU {% cite Damian2022 Ba2023 Oko2024 -f private_feature_learning_for_single_index_models.bib %}, we have the following result for $$\tanh$$.

**Theorem 6.** For $$f$$ a degree-$$q$$ polynomial, $$\eta, \gamma > 0$$, there exists a function $$a$$ such that for $$f_\gamma(x) = \E_{b \sim \N(0, 1)} [a(b)\sigma(\eta \cdot z + b)]$$ and $$\phi(b) = \frac{1}{\sqrt{2\pi}}e^{-\frac{b^2}{2}}$$ being the probability mass function for $$b \sim \N(0, 1)$$, we have

$$
\begin{equation}
\begin{gathered}
    \E_{z \sim \N(0, 1)}[(f(z) - f_\gamma(z))^2] < 2\gamma, \\
    \sup_b \left| a(b) \phi(b)\right| \lesssim \gamma^{-\frac{q-2}{2}} \eta \cdot e^{\frac{\pi^2 \gamma}{4}} \cdot \max\left[1, \left(\frac{\gamma}{\eta}\right)^{q-1}\right].
\end{gathered}
\end{equation}
$$

_Proof._ The idea is to interpret $$\E_{b \sim \N(0, 1)}[a(b)\sigma(\eta \cdot z + b)]$$ as a convolution, then as we pass to and from Fourier transform, we can construct $$a$$. This is done in four steps.
	
**Step 1: Convolution theorem.**
	
Let $$f = \sum_n c_n x^n$$. For $$\gamma > 0$$, let $$f_\gamma = f \cdot e^{-\gamma x^2}$$, then

$$
\begin{equation}
\begin{aligned}
    \mathbb{E}_{x} [(f - f_a)^2] 
    & = \mathbb{E}_{x} [[f(1 - e^{-ax^2})]^2] \\
    & \leq \left(\sup_{x} \frac{1}{\sqrt{2\pi}} |f(x)|^2 e^{-\frac{x^2}{2}}\right) \mathbb{E}_x [(1 - e^{-ax^2})^2] \\
    & \lesssim \mathbb{E}_x [(1 - e^{-ax^2})^2].
\end{aligned}
\end{equation}
$$
    
Now,

$$
\begin{equation}
    \begin{aligned}
        \mathbb{E}_x [(1 - e^{-\gamma x^2})^2]
        & = \frac{1}{\sqrt{2\pi}} \int_\mathbb{R} \left(1 - 2e^{-\gamma x^2} + e^{-2 \gamma x^2}\right) e^{-\frac{x^2}{2}} \,dx \\
        & = 1 - \frac{2}{\sqrt{2\pi}} \int_\mathbb{R} e^{-\gamma x^2} e^{-\frac{x^2}{2}} \,dx \\
        & \phantom{=} + \frac{1}{\sqrt{2\pi}} \int_\mathbb{R}e^{-2\gamma x^2} e^{-\frac{x^2}{2}} \,dx \\
        & = 1 - \frac{2}{\sqrt{2\pi}} \sqrt{2\pi \frac{1}{1 + 2\gamma}} + \frac{1}{\sqrt{2\pi}} \sqrt{2\pi \frac{1}{1+4\gamma}} \\
        & = 1 - \frac{2}{\sqrt{1+2\gamma}} + \frac{1}{\sqrt{1 + 4\gamma}} \\
        & < 1 - \frac{1}{\sqrt{1 + 2\gamma}} < 2\gamma
    \end{aligned}
\end{equation}
$$

Note that 

$$
\begin{equation}
    \mathbb{E}_{b \sim \mathcal{N}(0, 1)} [a(b) \sigma(\eta z - b)] = -\mathbb{E}_{b \sim \mathcal{N}(0, 1)} [a(b) \sigma(-\eta z + b)]
\end{equation}
$$
    
so we can work the the former, merely by changing $$f(x)$$ with $$-f(-x)$$. Now we construct $$a$$ such that $$\mathbb{E}_{b \sim \mathcal{N}(0, 1)} [a(b) \sigma(z - b)] = f_\gamma(z)$$. Due to the factor $$e^{-ax^2}$$, $$f_\gamma$$ is a Schwartz function; suppose $$a$$ is also Schwartz, we can use convolution theorem. Note that by time scaling property, if we denote $$\sigma_\eta(x) = \sigma(\eta x)$$ then $$\widehat{\sigma_\eta}(t) = \frac{1}{\eta} \widehat{\sigma}\left(\frac{t}{\eta}\right)$$.

$$
\begin{equation}
    \begin{aligned}
        \widehat{a \cdot \phi}(t) \cdot \frac{1}{\eta}\widehat{\sigma}\left(\frac{t}{\eta}\right) 
        & = \widehat{a\cdot\phi}(t) \cdot \widehat{\sigma_\eta}(t) \\
        & = \widehat{f_\gamma}(t) = \sum_n c_n \widehat{x^n e^{-\gamma x^2}} \\
        & = \sqrt{\frac{\pi}{\gamma}} \sum_n c_n \left(\frac{i}{2\pi}\right)^n \frac{d^n}{dt^n}e^{-\pi^2 \frac{t^2}{\gamma}} \Big|_t \\
        & = \sqrt{\frac{\pi}{\gamma}} \sum_n c_n \left(-\frac{i}{2\pi}\right)^n  \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right) e^{-\pi^2 \frac{t^2}{\gamma}}
    \end{aligned}
\end{equation}
$$

**Step 2: Calculating $$\widehat{\sigma}(t)$$.**

Unfortunately, $$\tanh$$ is not in $$L^1 (\mathbb{R})$$ so we cannot calculate its Fourier transform via definition. That being said, it is possible to calculate it as a distribution.

First, we calculate the Fourier transform of $$\sech$$, which, as $$\sech$$ is in $$L^1(\mathbb{R})$$, is well-defined. In fact, $$\sech$$ is even, and

$$
\begin{equation}
    |\sech(x)| = \frac{2e^{-|x|}}{1+e^{-2|x|}} \leq 2e^{-|x|} \in L^1(\mathbb{R}), 
\end{equation}
$$

whence follows

$$
\begin{equation}
    \int_\mathbb{R} \sech(x) e^{-i2\pi tx} \, dx = \int_\mathbb{R} \sech(x) \cos(2\pi tx)\, dx. 
\end{equation}
$$

Consider the rectangular contour $$C_R$$ connecting, in order, $$-R$$, $$R$$, $$R+i\pi$$, $$-R + i\pi$$, i.e. in counterclockwise fashion. The function $$g(z) = \frac{2\cos(2\pi tz)}{e^z + e^{-z}}$$ has only one simple pole in the region bounded by the contour at $$z = \frac{\pi}{2}i$$. Residue theorem gives

$$
    \int_{C_R} \frac{2\cos(2\pi tz)}{e^z + e^{-z}} \, dz = 2\pi i \Re_{\frac{\pi}{2}i}(g) = 2\pi \lim_{z \to \frac{\pi}{2} i} \left(z - \frac{\pi}{2}i\right) g(z) = \pi \left(e^{\pi^2t} + e^{-\pi^2 t}\right).
$$

Since $$g$$ is even, the contour integral can be written as

$$
\begin{equation}
    \begin{aligned}
        \int_{C_R} g(z) \, dz
        & = \int_{-R}^R g(x) \, dx + i \int_0^\pi g(R + ix) \, dx - \int_{-R}^R g(x + i\pi) \, dx \\
        & = \frac{2 + e^{2\pi^2 t} + e^{-2\pi^2 t}}{2} \int_{-R}^R g(x) \, dx \\
        & \phantom{=} + i \int_0^\pi \frac{e^{2\pi it} R e^{-2\pi tx} + e^{-2\pi itR} e^{2\pi tx}}{e^{-R}e^{itx} + e^R e^{-itx}} \, dx.
    \end{aligned}
\end{equation}
$$

Now letting $$R$$ tend to infinity and using dominated convergence theorem, we have

$$
\begin{equation}
    \frac{1}{2} \left(e^{\pi^2 t} + e^{-\pi^2 t}\right)^2 \widehat{\sech}(t) = \pi \left(e^{\pi^2 t} + e^{-\pi^2 t}\right),
\end{equation}
$$

or, $$\widehat{\sech}(t) = \pi \sech(\pi^2 t)$$. By convolution theorem, we have $$\widehat{\sech^2} = \widehat{\sech} \star \widehat{\sech}$$, meaning we have to calculate the integral

$$
\begin{equation}
    I = \int_\mathbb{R} \sech\left(\pi^2(t - y)\right) \sech\left(\pi^2 y\right) \, dy.
\end{equation}
$$

Its Cauchy principal value reads

$$
\begin{equation}
    I = \left(\frac{1}{\pi^2}\csch(\pi^2 t) \ln \frac{\csch(\pi^2 y)}{\csch(\pi^2 (y-t))}\right) \Big|_{y=-\infty}^\infty = 4t \csch(\pi^2 t),
\end{equation}
$$

Dividing by $$it$$, we have $$\widehat{\tanh}(t) = -4i \csch(\pi^2 t)$$.

**Step 3: Reconstructing $$a(t)$$.**

We now have

$$
\begin{equation}
    \widehat{a\cdot\phi}(t) = \eta\frac{i}{8} \sqrt{\frac{\pi}{\gamma}} \sum_n c_n \left(-\frac{i}{2\pi}\right)^n \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right) e^{-\pi^2 \frac{t^2}{\gamma}} \left(e^{\pi^2 \frac{t}{\eta}} - e^{-\pi^2 \frac{t}{\eta}}\right).
\end{equation}
$$

Now it is clear why we multiple $$f$$ with $$e^{-\gamma x^2}$$: being Gaussian up to a multiplicative function, we now have $$e^{-\pi^2 \frac{t^2}{\gamma}}$$ in the expression of $$\widehat{a}$$, which dominates the term $$\left(e^{\pi^2 \frac{t}{\eta}} - e^{-\pi^2 \frac{t}{\eta}}\right)$$ for all choice of $$\gamma$$ and $$\eta$$. This makes $$\widehat{a}$$, and in turn $$a$$, a Schwartz function, so that all calculation in this proof is justified in classical sense, and not only in distribution sense.

With some more pain, let us calculate $$a$$ explicitly. First, let us consider

$$
\begin{equation}
\begin{aligned}
    \Delta^+_n
    & = \frac{\pi}{\sqrt{\gamma}} \left(-\frac{i}{2\pi}\right)^n \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right) e^{-\pi^2 \frac{t^2}{\gamma}} e^{\pi^2 \frac{t}{\eta}} \\
    & = \frac{\pi}{\sqrt{\gamma}} \left(-\frac{i}{2\pi}\right)^n \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right) e^{-\pi^2 \frac{t^2}{\gamma} + \pi^2 \frac{t}{\eta}}
\end{aligned}
\end{equation}
$$

Completing the square

$$
\begin{equation}
    \pi^2 \frac{t^2}{\gamma} - \pi^2 \frac{t}{\eta} = \left(\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right)^2 - \frac{\pi^2 \gamma}{4\eta^2}
\end{equation}
$$

thus

$$
\begin{equation}
    \label{eq: first layer - expression of Delta_n^+}
    \Delta^+_n = \frac{\pi}{\sqrt{\gamma}} e^{\frac{\pi^2 \gamma}{4}} \left(-\frac{i}{2\pi}\right)^n \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right) e^{-\left[\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right]^2}.
\end{equation}
$$

Now expanding $$\h_n$$ in order to match the exponent in the exponentiation, 

$$
\begin{equation}
    \begin{aligned}
        (-1)^n \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right)
        & = (-1)^n \sum_{k = 0}^n \binom{n}{k} \h_k \left(\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right) \left[-\frac{\pi \sqrt{\gamma}}{\eta}\right]^{n-k} \\
        & = \sum_{k = 0}^n \binom{n}{k}(-1)^k \h_k \left(\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right) \left[\frac{\pi \sqrt{\gamma}}{\eta}\right]^{n-k} \\
    \end{aligned}
\end{equation}
$$

which, combining with the definition of $$\Delta^+_n$$, gives

$$
\begin{equation}
    \begin{aligned}
        \Delta^+_n
        = & \frac{\pi}{\sqrt{\gamma}} e^{\frac{\pi^2 \gamma}{4}} e^{-\left[\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right]^2} \left(\frac{i}{2\pi}\right)^n \\
        & \phantom{=} \sum_{k = 0}^n \binom{n}{k}(-1)^k \h_k \left(\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right) \left[\frac{\pi \sqrt{\gamma}}{\eta}\right]^{n-k} \\
        = &  e^{\frac{\pi^2 \gamma}{4}} \sum_{k = 0}^n \binom{n}{k} \left[\frac{i}{2\pi}\frac{\pi \sqrt{\gamma}}{\eta}\right]^{n-k} \cdot \left(\frac{i}{2\pi}\right)^k \frac{\pi}{\sqrt{\gamma}} (-1)^k \\
        & \phantom{=} \h_k \left(\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right) e^{-\left[\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right]^2} \\
        = & e^{\frac{\pi^2 \gamma}{4}} \sum_{k = 0}^n \binom{n}{k} \left[\frac{i\sqrt{\gamma}}{2\eta}\right]^{n-k} \cdot \left(\frac{i}{2\pi}\right)^k \frac{d^k}{dt^k} e^{-\left[\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right]^2}.
    \end{aligned}
\end{equation}
$$

Then we calculate the inverse Fourier transform of $$e^{-\left[\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right]^2}$$. For $$e^{-(at-b)^2}$$, its inverse Fourier transform is 

$$
\begin{equation}
    \frac{\sqrt{\pi}}{a} e^{-\frac{\pi^2}{a^2} x^2 + i2\pi bx},
\end{equation}
$$

and in particular, for $$e^{-\left[\frac{\pi}{\sqrt{\gamma}} t - \frac{\pi \sqrt{\gamma}}{2\eta}\right]^2}$$, we have its inverse Fourier transform to be

$$
\begin{equation}
    \frac{\sqrt{\pi}}{\frac{\pi}{\sqrt{\gamma}}} e^{-\frac{\pi^2}{\left(\frac{\pi}{\sqrt{\gamma}}\right)^2} x^2 + i2\pi \frac{\pi \sqrt{\gamma}}{2\eta}x} = \sqrt{\frac{\gamma}{\pi}}e^{-\gamma x^2 + i\pi^2\frac{\sqrt{\gamma}}{\eta}}
\end{equation}
$$

and so

$$
\begin{equation}
\begin{aligned}
    \check{\Delta^+_n} 
    & = e^{\frac{\pi^2 \gamma}{4}} \sum_{k = 0}^n \binom{n}{k} \left[\frac{i \sqrt{\gamma}}{2\eta}\right]^{n-k} \cdot x^k \sqrt{\frac{\gamma}{\pi}}e^{-\gamma x^2 + i\pi^2\frac{\sqrt{\gamma}}{\eta}} \\
    & = \sqrt{\frac{\gamma}{\pi}} e^{\frac{\pi^2 \gamma}{4}-\gamma x^2 + i\pi^2\frac{\sqrt{\gamma}}{\eta}} \left[x+\frac{i \sqrt{\gamma}}{2\eta}\right]^n
\end{aligned}
\end{equation}
$$

Similar calculation for

$$
\begin{equation}
    \Delta^-_n = \frac{\pi}{\sqrt{\gamma}} \left(-\frac{i}{2\pi}\right)^n \h_n \left(\frac{\pi}{\sqrt{\gamma}} t\right) e^{-\pi^2 \frac{t^2}{\gamma} - \pi^2 \frac{t}{\eta}}
\end{equation}
$$

yields

$$
\begin{equation}
    \check{\Delta^-_n} = \sqrt{\frac{\gamma}{\pi}} e^{\frac{\pi^2 \gamma}{4}-\gamma x^2 - i\pi^2\frac{\sqrt{\gamma}}{\eta}} \left[x-\frac{i \sqrt{\gamma}}{2\eta}\right]^n
\end{equation}
$$

Now plugging everything into $$ \widehat{a\cdot\phi}(t) = \frac{i\eta}{8\sqrt{\pi}} \sum_n c_n (\Delta^+_n - \Delta^-_n)$$, 

$$
\begin{equation}
\begin{aligned}
    a(x) \phi(x)
    & = \frac{i\eta\sqrt{\gamma}}{8\pi} e^{\frac{\pi^2 \gamma}{4}-\gamma x^2}\sum_n c_n \left[e^{i\pi^2 \frac{\sqrt{\gamma}}{\eta}}\left(x+\frac{i \sqrt{\gamma}}{2\eta}\right)^n \right.\\
    & \phantom{=} \left. - e^{-i\pi^2 \frac{\sqrt{\gamma}}{\eta}}\left(x-\frac{i \sqrt{\gamma}}{2\eta}\right)^n \right].
\end{aligned}
\end{equation}
$$

**Step 4. Bounding $\sup_{b \in \mathbb{R}} \|a(b) \phi(b)\|$.**

One has a violent bound

$$
\begin{equation}
    \begin{aligned}
        |a(x)\phi(x)|
        & \leq \frac{\eta\sqrt{\gamma}}{8\pi} e^{\frac{\pi^2 \gamma}{4}-\gamma x^2}\sum_{n = 0}^q |c_n| \left|\left|x+\frac{\sqrt{\gamma}}{2\eta}\right|^n - \left|x-\frac{\sqrt{\gamma}}{2\eta}\right|^n\right| \\
        & = \frac{\eta\sqrt{\gamma}}{8\pi} e^{\frac{\pi^2 \gamma}{4}-\gamma x^2}\sum_{n = 0}^q |c_n| \left|\sum_{i = 0}^{n-1} x^{n-1-i} \left(\frac{\sqrt{\gamma}}{2\eta}\right)^i \left[1 - (-1)^i\right] \right| \\
        & \leq \frac{\eta\sqrt{\gamma}}{4\pi} e^{\frac{\pi^2 \gamma}{4}-\gamma x^2}\sum_{n = 0}^q |c_n| \left|\sum_{i = 0}^{n-1} x^{n-1-i} \left(\frac{\sqrt{\gamma}}{2\eta}\right)^i \right|\\
    \end{aligned}
\end{equation}
$$

Consider $$g(x) = e^{-\gamma x^2} x^k$$. We have $$g'(x) = kx^{k-1}e^{-\gamma x^2} - 2\gamma x e^{-\gamma x^2} x^k$$. Solving for $$g'(x) = 0$$, either $$x = 0$$, which gives $$g(x) = 0$$, so eliminated, or $$2\gamma x^2 = k$$. The equation has two roots, at the positive of which $$g$$ attains its maximum, meaning at $$x = \sqrt{\frac{k}{2\gamma}}$$. So, 

$$
\begin{equation}
\begin{aligned}
    g(x) 
    & \leq \left(\sqrt{\frac{k}{2\gamma}}\right)^n e^{-\gamma \frac{k}{2\gamma}} = \left(\frac{k}{2e\gamma}\right)^{\frac{k}{2}} \\
    \Rightarrow e^{-\gamma x^2} x^{n-1-i} \left(\frac{\sqrt{\gamma}}{2\eta}\right)^i 
    & \leq \left(\frac{n-1-i}{2e\gamma}\right)^{\frac{n-1-i}{2}} \left(\frac{\sqrt{\gamma}}{2\eta}\right)^i \\
    & \lesssim \frac{\gamma^{i - \frac{n-1}{2}}}{\eta^i}.
\end{aligned}
\end{equation}
$$

Thus to conclude

$$
\begin{equation}
    |a(x)\phi(x)| \lesssim \gamma^{-\frac{q-2}{2}} \eta \cdot e^{\frac{\pi^2 \gamma}{4}} \cdot \max\left[1, \left(\frac{\gamma}{\eta}\right)^{q-1}\right].
\end{equation}
$$

### 3.4. Certificate construction.

Now the construction of our certificate is done in three steps.

**Step 1.** First, we consider the perfect infinite model, $$f_\infty(x) = \E_{b \sim \N(0, 1)}[a(b)\sigma(\langle \mu, x \rangle + b)]$$. Let $$\gamma = d^{-\frac{\varepsilon_p}{4(q+1)}}$$, by Theorem 6 for $$\eta = 1$$, note that as $$\langle \mu, x \rangle \sim \N(0, 1)$$, there exists a function $$a$$ such that

$$
\begin{equation}
	\label{eq: first layer - step 1}
	\E_{x \sim \N(0, \id)} \left[(f_*(x) - f_\infty(x))^2\right] < 2d^{-\frac{\varepsilon_p}{4(q+1)}}.
\end{equation}
$$

**Step 2.** Next, we consider the perfect finite model, $$f_p(x) = \sum_{i = 1}^p \overline{a}_i \sigma(\langle \mu, x \rangle + b_i)$$. Let $$C_b > 0$$ be some _fixed_ constant, and denote $\mathcal{I} = \\{i \mid \|b_i\| \leq C_b\\}$. Due to the standard normal initialisation of $$b_i$$, there exists a constant $$c_b > 0$$ depending only on $$C_b$$ such that $\|\mathcal{I}\| \geq_\P c_b \cdot p$.

First let us consider the case where $\|\langle \mu, x \rangle\| \leq \log d$. Denote $\H_{\log d} = \\{x \mid \|\langle \mu, x\rangle\| \leq \log d\\}$, and let $\overline{a}_i = \frac{1\_\{i \in \mathcal{I}\}}{\|\mathcal{I}\|} a(b_i)$, then

$$
\begin{equation}
	f_p(x) = \frac{1}{|\mathcal{I}|} \sum_{i \in \mathcal{I}} a(b_i) \sigma(\langle \mu, x \rangle + b_i).
\end{equation}
$$

Note that by Theorem 6 and the definition of $$\overline{a}_i$$ and $$\mathcal{I}$$, we have 

$$
\begin{equation}
	\label{eq: first layer - bound on a}
	\left| a(b_i) \sigma(\langle \mu, x \rangle + b_i)\right| \leq |a(b_i)| \lesssim \left(d^\frac{\varepsilon_p}{4(q+1)}\right)^{q-1-\frac{q-2}{2}} = \left(d^\frac{\varepsilon_p}{4(q+1)}\right)^{\frac{q}{2}} \ll d^\frac{\varepsilon_p}{8}.
\end{equation}
$$

which means if we replace one $$b_i$$ by another $$b'_i$$, we will have that the resulting model does not differ too much. In particular, let $$f'_p$$ be $$f_p$$ but with $$b'_i$$ replacing $$b_i$$, $F = \sup_x \|f_p(x) - f_\infty(x)\|$, and $F' = \sup_x \|f'_p(x) - f_\infty(x)\|$, then 

$$
\begin{equation}
	|F - F'| \lesssim \frac{2}{|\mathcal{I}|} d^\frac{\varepsilon_p}{8} \leq_\P \frac{2}{c_b \cdot p} d^\frac{\varepsilon_p}{8} \asymp \frac{1}{c_b \cdot p} d^\frac{\varepsilon_p}{8}.
\end{equation}
$$

Thus by McDiarmid's inequality {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %}, Assumption 2, and letting $$t = d^{-\frac{\varepsilon_p}{4q}}$$, 

$$
\begin{equation}
	\label{eq: first layer - bound of F-EF}
	\P \left[|F - \E[F]| \geq t \right] \leq 2 \exp\left(-\frac{2t^2}{p \frac{4}{c^2_b \cdot p^2} d^\frac{\varepsilon_p}{4}}\right) = 2\exp\left(-\frac{c^2_b \cdot pt^2}{2d^\frac{\varepsilon_p}{4}}\right)
\end{equation}
$$

which implies $\|F - \E[F]\| \leq_\P d^{-\frac{\varepsilon_p}{4(q+1)}}$.

To bound $$\E[F]$$, consider the function $g_x(b) = a(b) \mathbb{1}\_{\|b\| \leq C_b} \cdot \sigma(\langle \mu, x\rangle + b)$ which defines a function class $$\mathcal{G}$$. Let $$\varepsilon_i \sim \mathcal{U}(\{1, -1\})$$, Rademacher's inequality and the definition of $$\mathcal{I}$$ give

$$
\begin{equation}
	\begin{aligned}
		\E[F]
		& = \E_b \left[\sup_{x \in \H_{\log d}} \left| \E_b(g_x(b)) - \frac{1}{|\mathcal{I}|} \sum_{i \in \mathcal{I}} g_x(b_i)\right|\right] \\
        & = \E_b \left[\sup_{g \in \mathcal{G}} \left| \E_b(g(b)) - \hat{E}_b g(b)\right|\right] \\
        & \leq  2\E_{b, \varepsilon_1, \ldots, \varepsilon_p} \left[\sup_{g \in \mathcal{G}} \left| \frac{1}{|\mathcal{I}|} \sum_{i \in \mathcal{I}} \varepsilon_i g(b_i)\right|\right].
	\end{aligned}
\end{equation}
$$

Then by Lipschitz property of $$\sigma$$, 

$$
\begin{equation}
	\begin{aligned}
		\mathbb{E}[F]
		& \leq 2 \mathbb{E}_{b_i, \varepsilon_i} \left[\left| \frac{1}{|\mathcal{I}|}\sum_{i \in \mathcal{I}} \varepsilon_i a(b_i) \sigma(b_i)\right| + \sup_{x \in \H_{\log d}} \left| \frac{1}{|\mathcal{I}|}\sum_{i \in \mathcal{I}} \varepsilon_i a(b_i) \eta \langle \mu, x\rangle\right|\right] \\
		& = 2 \mathbb{E}_{b_i, \varepsilon_i} \left[\left|\frac{1}{|\mathcal{I}|} \sum_{i \in \mathcal{I}} \varepsilon_i a(b_i) \sigma(b_i)\right| + \eta \log d \left| \frac{1}{|\mathcal{I}|}\sum_{i \in \mathcal{I}} \varepsilon_i a(b_i) \right|\right] \\
	\end{aligned}
\end{equation}
$$

Now

$$
\E_{b_i, \varepsilon_i} \left[\varepsilon_i a(b_i) \mathbb{1}_{|b_i| \leq C_b} \sigma(b_i)\right] = \E_{b_i, \varepsilon_i} [\varepsilon_i a(b_i) \mathbb{1}_{|b_i| \leq C_b}] = 0,
$$

and

$$
    |\varepsilon_i a(b_i) \mathbb{1}_{|b_i| \leq C_b} \sigma(b_i)| \leq |\varepsilon_i a(b_i) \mathbb{1}_{|b_i| \leq C_b}| \lesssim d^\frac{\varepsilon_p}{8}
$$


so by standard concentration inequality and $$\eqref{eq: first layer - bound of F-EF}$$.

$$
\begin{equation}
	\E[F] \lesssim \frac{d^\frac{\varepsilon_p}{8} \log d}{\sqrt{\mathcal{I}}} \lesssim \frac{d^\frac{\varepsilon_p}{8} \log d}{\sqrt{p}}  \Rightarrow |F| \leq_\P d^{-\frac{\varepsilon_p}{4(q+1)}}
\end{equation}
$$

which implies

$$
\begin{equation}
    \E_x[(f_\infty(x) - f_p(x))^2 \mathbb{1}_{|\langle \mu, x \rangle \leq \log d}] \leq |F|^2 \leq_\P d^{-\frac{\varepsilon_p}{2(q+1)}}.
\end{equation}
$$

As for the case $$|\langle \mu, x\rangle| > \log d$$, since $$f_\gamma(x) = f(x)\cdot e^{-\gamma \langle \mu, x \rangle^2} = \sigma_*(\langle \mu, x \rangle) e^{-\gamma \langle \mu, x \rangle^2}$$, we have $$f^2_\gamma(x) \lesssim \langle \mu, x \rangle^{2q} e^{-2\gamma \langle \mu, x \rangle^2}$$. For function $$f(z) = z^q e^{-2\gamma z}$$, $$f'(z) = (q z^{q-1} - 2\gamma z^q) e^{-2\gamma z}$$, so $$f'(z) = 0$$ when $$z = 0$$ or $$z = \frac{q}{2\gamma}$$, implying $$f(z) \leq \left(\frac{q}{2\gamma}\right)^q e^{-\frac{q}{2}} = \left(\frac{q}{2\sqrt{e}\gamma}\right)^q$$.
Thus

$$
\begin{equation}
\begin{aligned}
    \E\left[f_\gamma (x)^2 \mathbb{1}_{|\langle \mu, x \rangle| > \log d}\right] 
    & \leq \sup_x f_\gamma(x)^2 \P(|\langle \mu, x \rangle| > \log d) \\
    & \lesssim \left(\frac{q}{2\sqrt{e}\gamma}\right)^q e^{-\frac{\log^2 d}{2}} \asymp d^{\frac{\varepsilon_p}{4}} e^{-\frac{\log^2 d}{2}}.
\end{aligned}
\end{equation}
$$

And similarly, 

$$
\begin{equation}
\begin{aligned}
    \E\left[f_p (x)^2 \mathbb{1}_{|\langle \mu, x \rangle| > \log d}\right] 
    & \lesssim \P(|\langle \mu, x \rangle| > \log d) \sup_{b, x} \left| a(b) \mathbb{1}_{|b| \leq C_b} \sigma(\langle \mu, x \rangle + b) \right| \\
    & \lesssim e^{-\frac{\log^2 d}{2}} d^{\frac{\varepsilon_p}{8}}.
\end{aligned}
\end{equation}
$$

All in all, we have

$$
\begin{equation}
	\begin{aligned}
		\E_x[(f_\infty(x) - f_p(x))^2] 
		& = \E_x[(f_\infty(x) - f_p(x))^2 \mathbb{1}_{|\langle \mu, x \rangle \leq \log d}] \\
        & \phantom{=} + \E_x[(f_\infty(x) - f_p(x))^2 \mathbb{1}_{|\langle \mu, x \rangle \leq \log d}] \\
		& \lesssim_\P d^{-\frac{\varepsilon_p}{4(q+1)}} + \E\left[f_\infty(x)^2 \mathbb{1}_{|\langle \mu, x \rangle| > \log d}\right] \\
        & \phantom{=} + \E\left[f_p (x)^2 \mathbb{1}_{|\langle \mu, x \rangle| > \log d}\right] \\
		& \lesssim_\P d^{-\frac{\varepsilon_p}{4(q+1)}} + d^{\frac{\varepsilon_p}{4}} e^{-\frac{\log^2 d}{2}} + e^{-\frac{\log^2 d}{2}} d^{\frac{\varepsilon_p}{8}} \\
        & = O_\P\left(d^{-\frac{\varepsilon_p}{4(q+1)}}\right).
	\end{aligned}
\end{equation}
$$

or, to summarise,
 
$$
\begin{equation}
	\label{eq: first layer - step 2}
    \E_x[(f_\infty(x) - f_p(x))^2] = O_\P\left(d^{-\frac{\varepsilon_p}{4(q+1)}}\right).
\end{equation}
$$


**Step 3.** Finally, we consider the initialisation and the noise. Let $$f^W(x) = \sum_{i = 1}^p \overline{a}_i \sigma (\langle w^1_i, x \rangle + b)$$, then by $$\eqref{eq: first layer - rank-1 approximation}$$.

$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} \E_x\left[\left(f^W(x) - f_p(x)\right)^2\right] \\
		& = \mathbb{E}_x \left[\left(\frac{1}{|\mathcal{I}|}\sum_{i \in \mathcal{I}} a(b_i) \left[ \sigma (\langle w_i^1 , x \rangle + b_i) - \sigma(\langle \mu, x\rangle + b_i)\right] \right)^2 \right] \\
		& \leq \frac{1}{|\mathcal{I}|^2} \|a(b)\|^2 \cdot \sum_{i \in \mathcal{I}}  \mathbb{E}_x \left[ \left[ \sigma (\langle w_i^1 , x \rangle + b_i) - \sigma(\langle \mu, x\rangle + b_i)\right]^2 \right] \\
		& \lesssim \frac{d^{\frac{\varepsilon_p}{4}}}{p} \sum_{i \in \mathcal{I}}  \mathbb{E}_x \left[ \left\langle O_\P \left( \sqrt{\frac{d^{1 + \varepsilon_b}}{n}} \log^{q+1} d \left[1  \right. \right. \right. \right. \\
        & \phantom{=} \left. \left. \left. \left. + \log^2 d \frac{\sqrt{\log(1.25/\delta)}}{\varepsilon}\right] \right), x\right\rangle^2 \right]\\
		& = O_\mathbb{P} \left(\frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+2} d + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right),
	\end{aligned}
\end{equation}
$$

or, to summarise, 

$$
\begin{equation}
    \label{eq: first layer - step 3}
    \E_x\left[\left(f^W(x) - f_p(x)\right)^2\right] \\ = O_\mathbb{P} \left(\frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+2} d + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right)
\end{equation}
$$

By our choice of $$\varepsilon_b$$, Assumption 1 and Assumption 3, combining $$\eqref{eq: first layer - step 1}$$, $$\eqref{eq: first layer - step 2}$$, and $$\eqref{eq: first layer - step 3}$$, we have

$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} \E_x\left[\left(f^W(x) - f_*(x)\right)^2\right] \\
		& \lesssim \E_x\left[\left(f^W(x) - f_p(x)\right)^2\right] + \E_x\left[\left(f_p(x) - f_\infty(x)\right)^2\right] \\
        & \phantom{=} + \E_x\left[\left(f_\infty(x) - f_*(x)\right)^2\right] \\
		& = O_\P \left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right) \\
        & =  O_\P \left(d^{-\frac{\varepsilon_p}{8(q+1)}}\right)
	\end{aligned}
\end{equation}
$$


Moreover, by our construction of $\overline{a}\_i = \frac{\mathbb{1}\_{\|b_i\| \leq C_b}}{\|\mathcal{I}\|} a(b_i)$ and bound on $$a$$ via Theorem 6 (cf. \eqref{eq: first layer - bound on a}), we also have $$\|\overline{a}\|^2 \lesssim \frac{d^\frac{\varepsilon_p}{4}}{p}$$, which implies

$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} \E_x\left[\left(f^W(x) - f_*(x)\right)^2\right] + \lambda \|\overline{a}\|^2 \\
		& = O_\P \left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2} + d^{\frac{\varepsilon_p}{4} - 2\varepsilon_n}\right) \\
		& = O_\P\left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right).
	\end{aligned}
\end{equation}
$$

as desired.

## References

{% bibliography --cited -f private_feature_learning_for_single_index_models.bib %}