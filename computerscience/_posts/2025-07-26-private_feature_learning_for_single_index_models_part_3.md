---
layout: post
title: "Private feature learning for single-index models (part 3)"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

This is part 3/3 of my M1 internship per requirement from ENS, done at Institute of Science and Technology Austria, under the supervision of Simone Bombari and Prof. Marco Mondelli, to be published.

## 4. Second layer.

In this section, we only use the dataset $$X^a$$, so by abuse of notation, we write $$X$$, $$x_j$$, and $$y_j$$ to mean $$X^a$$, $$x^a_j$$, and $$y^a_j$$.

### 4.1. Differential privacy of $$a_T$$.

This follows from our choice of hyper-parameters (Assumption 4) and the design of gradient descent algorithm {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %}.

### 4.2. Ridge regression.

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

Consider the Hilbert space $$\H = \{f_a : x \mapsto a^\intercal \sigma(W_1^\intercal x + b) \mid a \in \R^d\}$$ corresponding to the second layer. Let $$f_\H = \arg \min_{f_a \in \H} \E_x\left[\left(f_a(x) - f_*(x)\right)\right]^2 + \lambda\|a\|^2$$, and $$a_\H$$ be the corresponding parameter generating $$f_\H$$, which must exist by Riesz representation theorem.

Now denote the feature matrix $$\Phi = \sigma\left(W_1^\intercal X + b \cdot \mathbb{1}^\intercal\right) \in \R^{p \times n}$$ and the associated kernel $$K_n = \frac{1}{n} \Phi \Phi^\intercal \in \R^{p \times p}$$. For ridge regression

$$
\begin{equation}
	f_\RR = \arg \min_{f_a \in \H} \L (f_a) = \arg \min_{f_a \in \H} \frac{1}{n} \sum_{i = 1}^n \left(f_a(x_i) - y_i\right)^2 + \lambda \| a\|^2,
\end{equation}
$$

we have $$f_\RR(x) = a_\RR^\intercal \sigma(W_1^\intercal x + b)$$, where $$a_\RR = \frac{1}{n} \left(K_n + \lambda \id\right)^{-1} \Phi Y$$ and $$Y = (y_1, \ldots, y_n)^\intercal$$.

First, we write

$$
\begin{equation}
	\label{eq:second layer - ridge regression - inequality}
	\begin{aligned}
        \mathcal{R}(f_\RR) 
        & = \E_x\left[\left(f_\RR(x) - f_*(x)\right)^2\right] \\
        & \lesssim \E_x\left[\left(f_\H(x) - f_*(x)\right)^2\right] + \E_x\left[\left(f_\RR(x) - f_\H(x)\right)^2\right].
    \end{aligned}
\end{equation}
$$


By Theorem 2, we have that

$$
\begin{equation}
    \label{eq:second layer - ridge regression - first term}
	\begin{aligned}
		& \phantom{=} \E_x\left[\left(f_\H(x) - f_*(x)\right)^2\right] \\
		& \leq \E_x\left[\left(f_\H(x) - f_*(x)\right)^2\right] + \lambda \|a_\H\|^2 \\
        & = \min_{f_a \in \H} \E_x\left[\left(f_a(x) - f_*(x)\right)\right]^2 + \lambda\|a\|^2 \\
		& \leq \E_x\left[\left(f^W(x) - f_*(x)\right)\right]^2 + \lambda\|\overline{a}\|^2 \\
        & = O_\P\left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right).
	\end{aligned}
\end{equation}
$$

For the second term, note that the function $$\L(f_\RR)$$ is in $$n$$ variables $$x_1, \ldots, x_n$$. Let $$x' \sim \N(0, 1)$$, $$X^{(i)} = (x_1^{(i)}, \ldots, x_n^{(i)})$$ by $$X$$ but with $$x_i$$ replaced by $$x'$$, $$a_\RR^{(i)}$$ be the result of ridge regression on $$X^{(i)}$$ and its corresponding labels, and $$f_\RR^{(i)} (x) = a_\RR^{(i)} \sigma(W_1^\intercal x + b)$$. By definition of ridge regression, we have


$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} \mathcal{L}\left(f^{(i)}_\text{RR}\right) - \mathcal{L}(f_\text{RR}) \\
		& = \underbrace{\frac{1}{n} \sum_{j = 1}^n \left(f^{(i)}_\text{RR}\left(x^{(i)}_j\right) - y_i\right)^2 + \lambda \| a^{(i)}_\text{RR} \|^2 - \frac{1}{n} \sum_{j = 1}^n \left(f_\text{RR}\left(x^{(i)}_j\right) - y_i\right)^2 - \lambda \| a_\text{RR} \|^2}_{\leq 0} \\
		& \phantom{=} + \frac{1}{n} \left[f_\text{RR}(x') - f^{(i)}_\text{RR}(x') + f^{(i)}_\text{RR}(x_i) - f_\text{RR}(x_i)\right].
	\end{aligned}
\end{equation}
$$

Note that $$\left|f_\text{RR}(x') - f^{(i)}_\text{RR}(x')\right| = \left|\left(a_\text{RR} - a^{(i)}_\text{RR}\right) \sigma(W_1^\intercal x' + b)\right| \leq \left\|a_\text{RR} - a^{(i)}_\text{RR}\right\|\sqrt{p}$$,
and similarly for $$f^{(i)}_\text{RR}(x_i) - f_\text{RR}(x_i)$$. Triangle inequality gives 

$$
\begin{equation}
\begin{aligned}
    \mathcal{L}\left(f^{(i)}_\text{RR}\right) - \mathcal{L}(f_\text{RR}) 
    & \leq \frac{1}{n} \left[\left|f_\text{RR}(x') - f^{(i)}_\text{RR}(x')\right| + \left|f^{(i)}_\text{RR}(x_i) - f_\text{RR}(x_i)\right|\right] \\
    & \leq 2 \frac{\sqrt{p}}{n}\left\|a_\text{RR} - a^{(i)}_\text{RR}\right\|.
\end{aligned}
\end{equation}
$$


Since $$\mathcal{L}$$ is $$2\lambda$$-strongly convex, we have

$$
\begin{equation}
	\label{eq:second layer - bound on a_RR^i - a_RR}
	\mathcal{L}(a^{(i)}_\text{RR}) - \mathcal{L}(a_{\text{RR}})\geq \lambda \|a^{(i)}_\text{RR} - a_{\text{RR}}\|^2 \Rightarrow \|a^{(i)}_\text{RR} - a_{\text{RR}}\| \leq \frac{2\frac{\sqrt{p}}{n}}{\frac{p}{d^{2\varepsilon_n}}} = \frac{2d^{2\varepsilon_n}}{n\sqrt{p}} .
\end{equation}
$$

McDiarmid's inequality {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %} then implies that

$$
\begin{equation}
	\mathbb{P}\left(\left\|a_\text{RR} - a_\mathcal{H}\right\| \geq t\right) = \mathbb{P}\left(\left\|a_\text{RR} - \mathbb{E}_X\left[a_\text{RR}\right]\right\| \geq t\right)
	\leq 2\exp\left(-\frac{2t^2}{n \frac{4d^{4\varepsilon_n}}{pn^2}}\right)
\end{equation}
$$

thus implies $$\left\|a_\text{RR} - a_\mathcal{H}\right\| \leq_\mathbb{P} \frac{d^{2\varepsilon_n}}{\sqrt{pn}} \log d$$, whence

$$
\begin{equation}
	\label{eq:second layer - ridge regression - second term}
	\begin{aligned}
        \mathbb{E}_x\left[f_\text{RR}(x) - f_\mathcal{H}(x)\right]^2
        & \leq \left\|a_\text{RR} - a_\mathcal{H}\right\|^2 \mathbb{E}_x\left[\left\|\sigma(W_1^\intercal x + b)\right\|^2\right] \\
	    & \leq_\mathbb{P} \frac{d^{4\varepsilon_n}}{pn} p \log^2 d = O\left(\frac{d^{4\varepsilon_n}}{n} \log^2 d\right).
    \end{aligned}
\end{equation}
$$

Combining $$\eqref{eq:second layer - ridge regression - inequality}$$, $$\eqref{eq:second layer - ridge regression - first term}$$, and $$\eqref{eq:second layer - ridge regression - second term}$$, we have

$$
\begin{equation}
	\label{eq:second layer - ridge regression test loss}
	\E_x\left[\left(f_\RR(x) - f_*(x)\right)^2\right] = O_\P\left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right).
\end{equation}
$$

### 4.3. Clipping analysis.

Now we will show that with high probability, clipping does not occur throughout the training of second layer, proceeding as did Brown et al. {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %}. Let $$\left(\tilde{a}_i\right)_{i = 1}^T$$ be the iterations of gradient descent with noise and clipping and $$f_i = \tilde{a}_i^\intercal \sigma(W_1^\intercal x + b)$$, that is

$$
\begin{equation}
	\tilde{a}_{i+1} = \tilde{a}_i - \eta_a \nabla_a \L(f_i) + \eta_a z_i.
\end{equation}
$$

where $$z_i \sim \N(0, \Sigma^2 \id)$$ and $$\Sigma = \frac{2C_a \sqrt{T}}{n} \frac{\sqrt{8 \log(1/\delta)}}{\varepsilon}$$. Since $$a_\RR = \frac{1}{n} (K_n + \lambda \id)^{-1} \Phi Y$$, letting $$L = 2K_n + 2\lambda\id$$ and expanding the definition, we have

$$
\begin{equation}
	\label{eq:second layer - recursion for tilde a}
	\begin{aligned}
		\tilde{a}_{i+1} 
		& = \tilde{a}_i -  \eta_a \left[L \tilde{a}_i - 2(K_n + \lambda \id)\frac{1}{n} (K_n + \lambda\id)^{-1} \Phi Y\right] + \eta_a z_i \\
		& = \tilde{a}_i - \eta_a [L \tilde{a}_i - L a_\RR] + \eta_a z_i \\
        & = (\id - \eta_a L)\tilde{a}_i + \eta_a L a_\RR + \eta_a z_i\\
		\Rightarrow \tilde{a}_{i+1} - \tilde{a}_i 
        & = (\id - \eta_a L)(\tilde{a}_i - a_\RR) + \eta_a z_i \\
        & = (\id - \eta_a L)^{i+1} (a_0 - a_\RR) + \eta_a z'_{i+1}
	\end{aligned}
\end{equation}
$$

where $$z'_i = \sum_{j = 1}^i z_i$$. Let $$\varphi_j = \sigma(W_1^\intercal x + b)$$ and $$g_{i, j} = (\varphi_j^\intercal \tilde{a}_i - y_j) \varphi_j + 2\lambda \tilde{a}_i$$ be the gradient of the loss on $$x_j$$ with respect to $$\tilde{a}_i$$, we show that simultaneously for all $$i$$ and $$j$$, we have that $$\|g_{i, j}\| \leq_\P C_a$$, meaning with high probability the clipping does not happen, by showing that $$\left\| 2\lambda \tilde{a}_i\right\| \leq_\P \frac{C_a}{2}$$ and $$\left\| (\varphi_j^\intercal \tilde{a}_i - y_j) \varphi_j\| \right\|\leq \frac{C_a}{2}$$.

#### 4.3.1 Bounding $$2\lambda \tilde{a}_i$$

Writing the singular value decomposition of $$\Phi = U\Lambda V^\intercal$$, we have

$$
\begin{equation}
    L = 2K_n + 2\lambda \id = 2 U \left(\frac{1}{n}\Lambda^2 + \lambda\id\right) U^\intercal \\
    \Rightarrow \id - \eta_a L = U \left[\id - 2\eta_a\left(\frac{1}{n}\Lambda^2 + \lambda\id\right)\right] U^\intercal
\end{equation}
$$

giving

$$
\begin{equation}
	\label{eq:second layer - bound on id - eta_a L}
	\lambda_i(\id - \eta_a L) = 1 - 2\eta_a\left(\frac{\lambda_i^2}{n} + \lambda\right) \\
    \Rightarrow \left\|\id - \eta_a L \right\| = 1 - 2\eta_a\left(\frac{\lambda_{\min}(\Phi)^2}{n} + \lambda\right) < 1 - 2\frac{d^{\varepsilon_n}\log^2 d}{p} \frac{p}{d^{2\varepsilon_n}} < 1
\end{equation}
$$

for $$d$$ sufficiently large. Using $$\eqref{eq:second layer - recursion for tilde a}$$ and triangle inequality, we have that

$$
\begin{equation}
\begin{aligned}
    \|\tilde{a}_i\| 
    & \leq \| a_\RR\| + \|\tilde{a}_i - a_\RR\| \\
    & \leq \| a_\RR\| +  \left\| \id - \eta_a L\right\|^i \|a_0 - a_\RR\| + \eta_a \|z'_i\| \\
    & \leq \|a_0\| + 2\|a_\RR\| + \eta_a\|z'_i\|.
\end{aligned}
\end{equation}
$$

Initialisation of $$a$$ gives

$$
\begin{equation}
	\label{eq:second layer - bound on a_0}
	\|a_0\| = 1 \ll \frac{d^{\varepsilon_n}\sqrt{p} \log^{q+1} d}{8}.
\end{equation}
$$

Next, using again the singular decomposition and Cauchy-Schwarz inequality give

$$
\begin{equation}
	(K_n + \lambda \id)^{-1} \Phi = n \cdot U(\Lambda^2 + \lambda n \id)^{-1} \Lambda V^\intercal \\
    \Rightarrow \left\|(K_n + \lambda \id)^{-1} \Phi\right\| = \frac{n\|\Phi\|}{\|\Phi\|^2 + \lambda n} \leq \frac{n \|\Phi\|}{2\|\Phi\|\sqrt{\lambda n}} = \frac{d^{\varepsilon_n}}{2}\sqrt{\frac{n}{p}}. 
\end{equation}
$$

On the other hand, since $$x_i \sim \N(0, \id)$$, $$\langle \mu, x_i \rangle \sim \N(0, 1)$$, and concentration of Gaussian variables {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %} gives $\|\langle \mu, x_i \rangle\| \leq_\P \log d$ for all $$i$$ simultaneously. Since $$\sigma_*$$ is a degree-$$q$$ polynomial, we thus have $$y_i = \sigma_* (\langle \mu, x_i \rangle) \lesssim \langle \mu, x_i \rangle^q \leq \log^q d$$, which implies $$\|Y\| \lesssim_\P \sqrt{n} \log^q d$$. Plugging these two observations back to the definition of $$a_\RR$$ gives

$$
\begin{equation}
\label{eq:second layer - bound on a_RR}
\begin{aligned}
    \|a_\RR\| 
    & \leq \frac{1}{n} \left\|(K_n + \lambda \id)^{-1} \Phi\right\| \|Y\| \\
    & \lesssim_\P \frac{d^{\varepsilon_n}}{n} \sqrt{\frac{n}{p}} \sqrt{n} \log^q d \\
    & = \frac{d^{\varepsilon_n}}{\sqrt{p}} \log^q d \\
    & \ll \frac{d^{\varepsilon_n}\sqrt{p} \log^{q+1} d}{8}.
\end{aligned}
\end{equation}
$$

Finally, to conclude, expanding the definition of $$z'_i$$, using triangle inequality, concentration bound on Gaussian variables {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %} and (\ref{eq:second layer - bound on id - eta_a L}), we have

$$
\begin{equation}
	\begin{aligned}
		\|z'_i\| 
        & \leq \sum_{k = 1}^i \left\|(\id - \eta_a L)^{i-1} z_{i - k}\right\| \\
        & \leq \sum_{k = 1}^i \left\|\id - \eta_a L \right\|^{i-1} \left\| z_{i - k}\right\| \\
        & \leq_\P \Sigma \sum_{k = 1}^i \left\|\id - \eta_a L \right\|^{i-1} \\
        & \leq \frac{\Sigma}{1 - \|\id - \eta_a L\|} \\
        & \lesssim \frac{\Sigma}{2\eta_a \lambda}
	\end{aligned}
\end{equation}
$$

implying with Assumption 3 and Assumption 4

$$
\begin{equation}
	\begin{aligned}
		\label{eq:second layer - bound on eta_a z'_i}
		\eta_a \|z'_i\| 
        & \lesssim \frac{\Sigma}{2\lambda} \\
		& = \frac{1}{2\frac{p}{d^{2\varepsilon_n}}}\frac{2C_a\sqrt{T}}{n} \frac{\sqrt{8\log(1/\delta)}}{\varepsilon} \\
        & = \frac{C_a d^\frac{\varepsilon_n}{2}}{n}\frac{d^{2\varepsilon_n}}{p} \frac{\sqrt{8\log(1/\delta)}}{\varepsilon} \\
		& = C_a\frac{d^{\frac{5}{2}\varepsilon_n}}{np} \frac{\sqrt{8\log(1/\delta)}}{\varepsilon} \\
        & = \frac{d^{\frac{5}{2}\varepsilon_n}\log^{q+1} d}{n \sqrt{p}} \frac{\sqrt{8\log(1/\delta)}}{\varepsilon} \\
        & \lesssim \frac{d^{\frac{5}{2}\varepsilon_n}\log^{q+2} d}{n \sqrt{p}} \\
        & \ll \frac{C_a}{8}.
	\end{aligned}
\end{equation}
$$


#### 4.3.2 Bounding $$(\varphi_j^\intercal \tilde{a}_i - y_j) \varphi_j$$.

Since $\|\sigma\| \leq 1$, one trivially has $$\|\varphi_j\| \leq \sqrt{p}$$, and it suffices to show that $\|\varphi_j^\intercal \tilde{a}_i - y_j\| \leq \frac{1}{2} d^{\varepsilon_n}\log^{q+1} d$$. We proceed by induction.

For the base case $$i = 0$$, we have $$\tilde{a}_0 = a_0 = \frac{1}{\sqrt{p}}(1, \ldots, 1)^\intercal$$, so $\|\varphi_j^\intercal a_0\| \leq 1$. On the other hand, we have shown above that $\|y_j\| \lesssim_\P \log^q d$, so overall

$$
\begin{equation}
	|\varphi_j^\intercal \tilde{a}_i - y_j| \leq |\varphi_j^\intercal \tilde{a}_i| + |y_j| \lesssim_\P \log^q d \ll \frac{1}{2} d^{\varepsilon_n}\log^{q+1} d.
\end{equation}
$$

For the induction step, we write

$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} |\varphi_j^\intercal \tilde{a}_i - y_j| \\
		& = |\varphi_j^\intercal \tilde{a}_i - \varphi_j^\intercal a_\text{RR} + \varphi_j^\intercal a_\text{RR} - y_j| \\
        & \leq |\varphi_j^\intercal \tilde{a}_i - \varphi_j^\intercal a_\text{RR}| + |\varphi_j^\intercal a_\text{RR} - y_j| \\
        & = |\varphi_j^\intercal \left(\tilde{a}_i -  a_\text{RR}\right)| + |\varphi_j^\intercal a_\text{RR} - y_j| \\
		& = |\varphi_j^\intercal \left[(\id - \eta L)^i \left(\tilde{a}_0 -  a_\text{RR}\right) + \eta_a z'_i\right]| + |\varphi_j^\intercal a_\text{RR} - y_j| \\
        & \leq |\varphi_j^\intercal (\id - \eta L)^i \left(a_0 -  a_\text{RR}\right)| + |\eta_a \varphi_j^\intercal z'_i| + |\varphi_j^\intercal a_\text{RR} - y_j| \\
		& \leq \underbrace{|\varphi_j^\intercal (\id - \eta L)^i \left(a_0 - a^{(i)}_\text{RR}\right)|}_{S_1} + \underbrace{|\varphi_j^\intercal (\id - \eta L)^i \left(a^{(i)}_\text{RR} - a_\text{RR}\right)|}_{S_2} \\
        & \phantom{=} + \underbrace{|\eta_a \varphi_j^\intercal z'_i|}_{S_3} + \underbrace{|\varphi_j^\intercal a_\text{RR} - y_j|}_{S_4}. \\
	\end{aligned}
\end{equation}
$$

For $$S_4$$, we combine the trivial bound $\|y_j\| \lesssim_\P \log^q d$ and $$\|\varphi_j\| \leq \sqrt{p}$$ from above and $$\eqref{eq:second layer - bound on a_RR}$$ to obtain

$$
\begin{equation}
\begin{aligned}
    |\varphi_j^\intercal a_\RR - y_j| 
    & \leq |\varphi_j^\intercal a_\RR| + |y_j| \\
    & \leq \|\varphi_j^\intercal\|  \|a_\RR\| + |y_j| \\
    & \lesssim_\P \sqrt{p} \frac{d^{\varepsilon_n}}{\sqrt{p}} \log^q d + \log^q d \\
    & \ll \frac{d^{\varepsilon_n} \log^{q+1} d}{8}.
\end{aligned}
\end{equation}
$$


For $$S_1$$, write $$a_0 - a_\RR^{(i)} = \psi v$$ and $$\varphi_j = \nu u$$ where $$u$$ and $$v$$ are unit vectors, note that $$v$$ is fixed after initialisation and $$u$$ is uniformly distributed on $$\mathbb{S}^{d-1}$$. Standard bound on inner product {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %} gives $\|\langle u, v\rangle\| \leq_\P \frac{\log^2 d}{\sqrt{p}}$. $$\eqref{eq:second layer - bound on id - eta_a L}$$, $$\eqref{eq:second layer - bound on a_0}$$, $$\eqref{eq:second layer - bound on a_RR}$$, and the trivial bound $$\nu \leq \sqrt{p}$$ give

$$
\begin{equation}
	|S_1| \leq_\P \frac{\log^2 d}{\sqrt{p}}\sqrt{p} = \log^2 d \ll \frac{d^{\varepsilon_n} \log^{q+1} d}{8}.  
\end{equation}
$$

For $$S_2$$, (\ref{eq:second layer - bound on a_RR^i - a_RR}) and (\ref{eq:second layer - bound on id - eta_a L}) implies

$$
\begin{equation}
\begin{aligned}
    |S_2| 
    & \leq \|\phi_j\| \|(\id - \eta_a L)^i\| \|a^{(i)}_\RR - a_\RR\| \\
    & \lesssim_\P \sqrt{p} \log^q d \frac{d^{2\varepsilon_n}}{n\sqrt{p}} \\
    & = \frac{d^{2\varepsilon_n}}{n} \log^q d \\
    & \ll \frac{d^{\varepsilon_n} \log^{q+1} d}{8}.
\end{aligned}
\end{equation}
$$

To conclude, for $$S_3$$, we write

$$
\begin{equation}
	|\eta_a \varphi_j^\intercal z'_i|
	= \eta_a \left| \sum_{k = 1}^i \varphi_j^\intercal (\id - \eta L)^{k-1} z_{i - k}\right| \leq \eta_a \sum_{k = 1}^i \left|\varphi_j^\intercal (\id - \eta L)^{k-1} z_{i - k}\right|. \\
\end{equation}
$$

Due to independence between $$\varphi_j$$ and $$z_{i-k}$$, we can fix $$\varphi_j$$ and have $$\varphi_j^\intercal (\id - \eta L)^{k-1} z_{i - k} \sim \|\varphi_j^\intercal (\id - \eta L)^{k-1}\| \Sigma \cdot \mathcal{N}(0, \id)$$
which gives the high probability bound

$$
\begin{equation}
\begin{aligned}
    \left|\varphi_j^\intercal (\id - \eta L)^{k-1} z_{i - k}\right| 
    & \leq_\mathbb{P} \|\varphi_j^\intercal (\id - \eta L)^{k-1}\| \Sigma \log d \\
    & \leq \|\id - \eta L\|^{k-1} \| \|\varphi_j\| \Sigma \log d.
\end{aligned}
\end{equation}
$$

together with the trivial bound $$\|\varphi_j\| \leq \sqrt{p}$$ implying

$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} \sum_{k = 1}^i \left|\varphi_j^\intercal (\id - \eta L)^{k-1} z_{i - k}\right| \\\
		& \leq_\P \|\varphi_j\| \Sigma \log d  \sum_{k = 1}^i  (\id - \eta L)^{k-1} \\
        & \leq \sqrt{p} \log d \frac{2C_a\sqrt{T}}{n} \frac{\sqrt{8\log(1/\delta)}}{\varepsilon} \frac{1}{1 - \|\id - \eta_a L\|}\\
		& \leq \sqrt{p} \log d \frac{2d^{\frac{\varepsilon_n}{2}}d^{\varepsilon_n}\sqrt{p}\log^{q+1} d}{n} \log d \frac{1}{2\eta_a \lambda} \\
		& \leq \frac{d^{\frac{3}{2}\varepsilon_n}\log^{q+3} d}{n} \frac{p}{d^{\varepsilon_n} \log^2 d} \frac{d^{2\varepsilon_n}}{p} \\
        & = \frac{d^{\frac{5}{2}\varepsilon_n}\log^{q+1} d}{n} \\
        & \ll \frac{d^{\varepsilon_n} \log^{q+1} d}{8}.
	\end{aligned}
\end{equation}
$$


### 4.5. Noise analysis and early stopping.

Knowing that the clipping does not occur, the final task is to show that the noise added and the fact that we stop the gradient descent after $$T$$ iteration do not hurt the performance, meaning that after $$T$$ iterations, we will have gotten sufficiently close to $$a_\RR$$. Using $$\eqref{eq:second layer - recursion for tilde a}$$, $$\eqref{eq:second layer - bound on id - eta_a L}$$, and $$\eqref{eq:second layer - bound on eta_a z'_i}$$ we have that

$$
\begin{equation}
\begin{aligned}
    \|\tilde{a}_T - a_\RR\| 
    & \leq \left\| \id - \eta_a L\right\|^T \|a_0 - a_\RR\| + \eta_a \|z'_T\| \\
    & \lesssim_\P \frac{d^{\frac{5}{2}\varepsilon_n}\log^{q+1} d}{n \sqrt{p}} \frac{\sqrt{\log(1/\delta)}}{\varepsilon} + \|a_0 - a_\RR\| \left(1 - 2\frac{\log^2 d}{d^{\varepsilon_n}}\right)^\intercal
\end{aligned}
\end{equation}
$$

where $$\left(1 - 2\frac{\log^2 d}{d^{\varepsilon_n}}\right)^T = \left(1 - 2\frac{\log^2 d}{d^{\varepsilon_n}}\right)^{d^{\varepsilon_n}} \leq \exp(-2\log^2 d)$$. Thus combining with $$\eqref{eq:second layer - ridge regression test loss}$$, we reach the conclusion of Theoreom 1

$$
\begin{equation}
	\begin{aligned}
		& \phantom{=} \mathcal{R}(f_\text{NN}) = \mathcal{R}(f_T) \\
		& = \E_x \left[\left(f_T(x) - f_*(x)\right)^2\right] \\
        & \lesssim \E_x \left[\left(f_T(x) - f_\RR(x)\right)^2\right] + \E_x \left[\left(f_\RR(x) - f_*(x)\right)^2\right] \\
		& \lesssim \|a_T - a_\RR\|^2 \E_x \left[\|\sigma(W_1^\intercal x + b)\|^2\right] + O_\P\left(\frac{1}{\log d}\right) \\
		& = \frac{d^{5\varepsilon_n}\log^{2q+2} d}{n^2 p} \frac{\log(1/\delta)}{\varepsilon^2} p \\
        & \phantom{=} + O_\P\left(d^{-\frac{\varepsilon_p}{8(q+1)}}  \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right)\\
		& = O_\P\left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_b + \frac{\varepsilon_p}{4}}}{n} \log^{2q+6} d \frac{\log(1.25/\delta)}{\varepsilon^2}\right) \\
        & = O_\P\left(d^{-\frac{\varepsilon_p}{8(q+1)}} + \frac{d^{1 + \varepsilon_n}}{n}\frac{\log(1.25/\delta)}{\varepsilon^2}\right).
	\end{aligned}
\end{equation}
$$

## 5. Conclusion.

### 5.1. Summary.
In this work, we have shown that with a sufficiently wide two-layer neural network using $$\tanh$$ as activation function, one can learn a single-index model whilst ensuring strong differential privacy. Our approach largely follow Damian et al. {% cite Damian2022 -f private_feature_learning_for_single_index_models.bib %} and Ba et al. {% cite Ba2023 -f private_feature_learning_for_single_index_models.bib %}, based on characterising the first layer after one gradient step and, using a construction of second layer with small regularised test loss, one can, via ridge regression, achieve the same. Then we follow strategy of Brown et al. {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %} to show that training second layer via differentially private gradient descent does not incur any performance loss, where to make up for the lack of an underlying linear data model, we adapt leave-one-out argument similar to that of Bombari and Mondelli {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %}.

It should be noted that the choice of parameters are not tight, and our proofs leave some room. For example, in Assumption 3, one can loosen the third condition to $$\frac{\varepsilon}{\sqrt{\log(1/\delta)}} = \omega\left(\frac{1}{\text{polylog} d}\right)$$ or even $$\omega\left(\frac{\log d}{d^{\varepsilon_n}}\right)$$.

Theorem 6 mirrors the universal approximation result on ReLU {% cite Damian2022 Ba2023 Oko2024 -f private_feature_learning_for_single_index_models.bib %}. For this reasons, we expect that the same result as ours holds for ReLU. Conversely, a number of other works based on this strategy have adapted ReLU as activation function due to the crucial reliance on the existence of such universal approximation result, as well as the fact that the proof for ReLU is specific and not easy to adapt to other activations. As $$\tanh$$ is another popular activation function due to $$\alpha_0 = 0$$ and $$\alpha_1 \neq 0$$ {% cite Bombari2025 Moniri2024  -f private_feature_learning_for_single_index_models.bib %}, we hope that Theorem 6 will help analysis along this line of argument.

### 5.2. Limitations and future works.

The proof of Theorem 6 remains specific to polynomial targets, whilst other works studying one step of gradient descent on the first layer consider more general targets, such as $$\Theta(1)$$-Lipschitz {% cite Ba2022 Moniri2024 -f private_feature_learning_for_single_index_models.bib %}. It is thus an open question to generalise this result for other activation functions.

One point that we find unsatisfactory is the different nature of gradient descent on two layers: only the first layer is normalised. This in fact can be remedied, but at the cost of $$\eta_W$$ and $$\sigma_W$$ depend on the $a_0$, $$\alpha^{b_i}_1$$, and most importantly, $$\alpha^*_1$$, which is unknown a priori and thus requires a bound on $$\alpha^*_1$$.

Another point is that the $$a_0$$ is initialised to be $$\frac{1}{\sqrt{p}}(1, \ldots, 1)^\intercal$$. If we instead sampled $$a_0$$ from $$\N\left(0, \frac{1}{p} \id\right)$$ to resemble initialisation of other parameters, the step size for the first layer, $$\eta_W$$ must be set to $\frac{d^{\frac{\varepsilon_n}{2}}}{\min\_{1 \leq i \leq p} \|a_{0, i}\|}$, so that $$w_i^0$$ is hidden in the error term of the rank-$1$ approximation. Since there is no effective concentration bound on Cauchy variables, $$\eta_W$$ must be set _a posteriori_.

One can also remedy this by first normalising the gradient $$\nabla_{w_i} \mathcal{L}$$, then, with $$\eta_W = \sqrt{\frac{n}{d^{1 + \varepsilon_b}}}$$, updating the first layer, then normalising _again_ before adding Gaussian noise. That way, our rank-$1$ approximation still holds, and the rest of the analysis remains the same. We chose not to do this due to the artificiality of training procedure at that point.

Our strategy to ensure differential privacy is similar to that of Bombari and Mondelli {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %} and Brown et al. {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %}, meaning we show that one can set the clipping constant large enough for the clipping not to happen at all with high probability, but also small enough so that the necessary noise added does not degrade the final performance. Interestingly, Bombari et al. {% cite Bombari2025a -f private_feature_learning_for_single_index_models.bib %} recently showed that for linear regression, it is more beneficial to allow more aggressive clipping, that is, the clipping constant to be of the same order as, if not smaller than, the expected norm of per-sample gradient.

Beyond single-index models, one may ask if privacy for free can be achieved for multi-index models. Unfortunately, Dandi et al. {% cite Dandi2024a -f private_feature_learning_for_single_index_models.bib %} showed that for shallow neural networks without unit bias, the landscape of learnable targets is more complicated, and even without differential privacy, it is not always learnable. In this case, one will have to take these constraints into account.

## References

{% bibliography --cited -f private_feature_learning_for_single_index_models.bib %}