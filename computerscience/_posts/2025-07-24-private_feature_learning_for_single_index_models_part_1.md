---
layout: post
title: "Private feature learning for single-index models (part 1)"
author: "Nguyen Doan Dai"
categories: computerscience
tags: math, computerscience
---

This is part 1/3 of my M1 internship per requirement from ENS, done at Institute of Science and Technology Austria, under the supervision of Simone Bombari and Prof. Marco Mondelli, to be published.

### Abstract

$$(\varepsilon, \delta)$$-differential privacy has been gold standard in gauging protection of training data. It was commonly believed that there exists some privacy-performance trade-off, due to more parameters to be kept private. But recently, it was shown that one needs not to sacrifice performance for privacy when working with sufficiently over-parameterised linear regression or random feature models. Nevertheless, neither of the two models are capable of feature learning, and it has been open whether one could achieve feature learning whilst guaranteeing privacy.
		
In this work, we show that with sufficiently wide shallow neural networks trained with quadratic loss, one can learn single-index models, i.e. $$\mathcal{R} = o(1)$$, whilst maintaining strong privacy, i.e. $$\varepsilon, \delta = o(1)$$.

### Acknowledgement

I thank Simone Bombari and Marco Mondelli for suggesting this problem, countless helpful discussions, warmest hospitality at Institute of Science and Technology Austria, and support throughout this project. I acknowledge financial help from Department of Computer Science, École Normale Supérieure - PSL which have made my trip to Vienna possible.

Last but not least, I have my utmost gratitude to my dear friend Anna-Maria Edlinger for her unwavering support through difficulties during my staying in Vienna, without which I would have not been able to finish this project.

As such, this work is dedicated to her.

## 1. Introduction.

### 1.1. Privacy.

As machine learning finds more applications in privacy-sensitive setting, such as in medicine {% cite Shehab2022 -f private_feature_learning_for_single_index_models.bib %}, communication {% cite Dada2019 -f private_feature_learning_for_single_index_models.bib %}, finance {% cite HernandezAros2024 -f private_feature_learning_for_single_index_models.bib %}, and personalised services {% cite Weiner2024 -f private_feature_learning_for_single_index_models.bib %}, arises the question of protecting data. Indeed, from a publicly available model, one might infer unexpected properties {% cite Ganju2018 -f private_feature_learning_for_single_index_models.bib %}, extracting {% cite Carlini2021 -f private_feature_learning_for_single_index_models.bib %} or reconstruct {% cite Fredrikson2015 -f private_feature_learning_for_single_index_models.bib %} training data, or test whether a sample is in the dataset {% cite Shokri2017 -f private_feature_learning_for_single_index_models.bib %}.

One possible approach to mitigate these attacks is to ensure $$(\varepsilon, \delta)$$-differential privacy {% cite Dwork2013 -f private_feature_learning_for_single_index_models.bib %}, for example by using differentially private variants of stochastic gradient descent {% cite Abadi2016}, but this training scheme might hurt the performance {% cite Tramer2021 De2022 -f private_feature_learning_for_single_index_models.bib %} and/or incur computational cost {% cite Kurakin2022 -f private_feature_learning_for_single_index_models.bib %}. Deferring the precise definition to later sections, $$\varepsilon$$ and $$\delta$$ represent the privacy budget: the smaller they are, the stronger the privacy guarantee. As no lunch should come free, one might expect some trade-off between performance and privacy {% cite Chaudhuri2009 Chaudhuri2011 -f private_feature_learning_for_single_index_models.bib %}, similar to that between bias and variance in classical theory. The more parameters to train, the harder it is to have both privacy and performance due to the noise level added {% cite Jayaraman2019 -f private_feature_learning_for_single_index_models.bib %}.

Unexpectedly, it was shown that in overparameterised regime, privacy can be achieved with $$\varepsilon, \delta = o(1)$$ and the same performance for linear regression {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %} and random feature models {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %}, thus arguing for benefits of overparameterisation.

### 1.2. Feature learning.

That being said, one cannot expect linear regression to compete with large neural networks used in practice. It serves as a good proof-of-concept to suggest that privacy and performance might be attainable at the same time, but it remains an open question to show this for more complicated model.

Likewise, random feature models have been shown to be equivalent to _noisy linear models_ with Gaussian features and Gaussian noise, known as _Gaussian equivalence property_ {% cite Hu2023 Bosch2023 -f private_feature_learning_for_single_index_models.bib %}. This means that a random feature model can learn only components up to some degree $$q$$ of the true model, where $$q$$ depends on random feature model's size. Or, in another word, such a model has its performance comes purely from its size, and it does not learn any features of the data, thus lacks _feature learning_.

This prompts the question of this work, that is: 

_Can one simultaneously achieve both feature learning and differential privacy?_

### 1.3. Contribution and outline.

In this work, we consider
- a shallow neural network where the number of parameters $$p$$ potentially significantly larger than the sample size $$n$$ and the input dimension $$d$$
- with the training process comprising one step of gradient descent on the first layer, then gradient descent on the second layer, with some modification to ensure differential privacy
- over data generated by a single-index model.

We show that indeed, as with linear regression and random feature models, under some technical conditions and for all sufficiently over-parameterised models, it is possible to achieve both $$(\varepsilon, \delta)$$-differential privacy and $$o(1)$$ generalisation error, as $$p$$, $$n$$, and $$d$$ suitably tend to infinity, thus give an affirmative response to the above question.

Section 2 sets rigorously the problem and the main result, introduces necessary notations and definitions, as well as provides a review of related works, before stating the outline of the proof. Section 3 and Section 4 analyse first and second layer respectively. Section 5 concludes with some final remarks. 
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
## 2. Preliminaries.

### 2.1. Differential privacy.

Dwork and Roth {% cite Dwork2013 -f private_feature_learning_for_single_index_models.bib %} defined $$(\varepsilon, \delta)$$-differential privacy (DP) based on the idea that if regardless of the dataset, with high probability changing one sample does not change an algorithm's output, an attacker has hard time inferring any information about one sample from the final output alone. This is stated as follow.

**Definition 1**. A randomised algorithm $$\mathcal{A}$$ is said to be $$(\varepsilon, \delta)$$-differentially private if for any subset $$S$$ of output space $$\text{Range}(\mathcal{A})$$ and any pair of adjacent dataset $$D, D' \subseteq \text{Dom}(\mathcal{A})$$, we have that

$$
\begin{equation}
    \P[\mathcal{A}(D') \in S] \leq \exp(\varepsilon) \P[\mathcal{A}(D) \in S] + \delta
\end{equation}
$$

One caveat in machine learning is that the trained parameters are hard to characterise, and in practice, the result of a random process due to the use of stochastic gradient descent. As such, it is a priori non-trivial to enforce $$(\varepsilon, \delta)$$-DP; for numerical algorithms, a common approach, so-called Gaussian mechanism is to add Gaussian noise to the final output {% cite Dwork2013 -f private_feature_learning_for_single_index_models.bib %}, thus obscure its relationship with input.

Since one can compose two differentially private algorithms to obtain another {% cite Dwork2013 -f private_feature_learning_for_single_index_models.bib %}, one can apply this mechanism for every gradient step. Unfortunately the final result will have too much noise to guarantee any performance; fortunately, Abadi et al. {% cite Abadi2016 -f private_feature_learning_for_single_index_models.bib %} showed that the noise level needs not be too large, which allows to recover good performance. Bombari and Mondelli {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %} showed that indeed this is the case for random feature models. It essentially involves training only the second layer whilst leaving the first layer fixed, which is precisely what we will do to train the second layer. Details of training procedure are given in later sections.

### 2.2. Setting.

In what follows, we consider the input in $$\R^d$$ and output in $$\R$$.

- Model: as stated above, we consider a fully connected two-layer neural network with $$p$$ hidden neurons, which takes the form

$$
\begin{equation}
    \label{eq:model}
    f(x) = a^\intercal \sigma(W^\intercal x + b)
\end{equation}
$$

where $$\sigma = \tanh$$ is the non-linear activation function, $$W \in \R^{d \times p}$$ is coefficients of first layer, $$b \in \R^p$$ is the bias units, and $$a \in \R^p$$ is the coefficients of second layer.
-  Dataset: we have _two_ disjoint sets of data of same size $$n$$, namely $$X^W = (x^W_1, \ldots, x^W_n)$$ to train the first layer, and $$X^a = (x^a_1, \ldots, x^a_n)$$ to train the second layer. For $$1 \leq i \leq n$$, $$x^W_i$$ and $$x^a_i$$ have respective labels $$y^W_i$$ and $$y^a_i$$. This is merely to ensure independence of each layer's training and facilitate analysis later on, as was commonly done {% cite Dandi2024 Ba2023 Berthier2024 Abbe2023 Damian2022 -f private_feature_learning_for_single_index_models.bib %}.
- Loss function: we consider the square loss, which takes the form

$$
\begin{equation}
    \ell(x, y, f) = (f(x) - y)^2,
\end{equation}
$$

and we train with respect to empirical loss and regularisation on the second layer, meaning we minimise

$$
\begin{equation}
    \L(f) = \frac{1}{n} \sum_{i = 1}^n \left[f(x_i) - y_i\right]^2 + \lambda \|a \|^2.
\end{equation}
$$

Training procedure comprises two stages given below, denoted Algorithm 1. Note that we separate the training of each layer, as was commonly done {% cite Dandi2024 Ba2023 Berthier2024 Abbe2023 Damian2022 -f private_feature_learning_for_single_index_models.bib %}, for the same reason we have two datasets: it is merely to help the analysis later on, and we expect the general picture to hold without resampling, barring worse sample complexity and/or slower rate of convergence for test loss.

___

**Require.** Number of iterations $$T$$, step sizes $$\eta_W, \eta_a$$, clipping constants $$C_a$$, noise level $$\sigma_w, \sigma_a$$.

**$$\rhd$$ Initialisation.** $$(W_0)_{i, j} \sim \N\left(0, \frac{1}{d}\right), b \sim \N(0, \id), a_0 = \frac{1}{\sqrt{p}}(1, \ldots, 1)^\intercal$$.

**$$\rhd$$ First layer training.**

**for** $$i = 1, \ldots, n$$ **do**

&nbsp;&nbsp;&nbsp;&nbsp; Compute the gradient $$g_W(x^W_i, y^W_i, f) \gets \nabla_W \L = \nabla_W (a_0^\intercal \sigma(W_0^\intercal x^W_i + b) - y^W_i)^2$$.

**endfor**

Update the first layer $$W_1 \gets W_0 - \eta_W \sum_{i = 1}^n g^C_W(x^W_i, y^W_i, f)$$

**for** $$i = 1, \ldots, p$$ **do**

&nbsp;&nbsp;&nbsp;&nbsp; Normalise the neurons and add noise $$w_i^1 \gets \frac{w_i^1}{\|w_i^1\|} + \sigma_W \N(0, \id).$$

**endfor**

**$$\rhd$$ Second layer training.**

**for** $$t = 1, \ldots, T$$ **do**

&nbsp;&nbsp;&nbsp;&nbsp; **for** $$i = 1, \ldots, n$$ **do**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Compute the gradient 

$$g_a(x^a_i, y^a_i, f) \gets  \nabla_a \L = \nabla_a \left[(a_{i-1}^\intercal \sigma(W_1^\intercal x^a_i + b) - y^a_i)^2 + \lambda \|a_{i-1}\|^2\right]$$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Clip the gradient

$$g^C_a(x^a_i, y^a_i, f) \gets g_a(x^a_i, y^a_i, f) \cdot \max\left(1, \frac{\left\|g^C_a(x^a_i, y^a_i, f)\right\|}{C_a}\right).$$

&nbsp;&nbsp;&nbsp;&nbsp; **endfor**

Aggregate the gradient $g_a \gets \frac{1}{n} \sum_{i = 1}^n g^C_a(x^a_i, y^a_i, f)$.

Update the second layer layer $a_t \gets a_{t-1} - \eta_a g_a + \sqrt{\eta_a} \frac{C_a}{n} \sigma_a \N(0, 1)$.

**endfor**

**Return** $$W_1, a_T$$.

___

<div align="center">Algorithm 1.</div>

Finally, we make some technical assumptions.

**Assumption 1.**(Data distribution) The inputs $$x^W_i$$'s and $$x^a_i$$'s are i.i.d sampled from standard Gaussian distribution $$\N(0, 1)$$. The labels are generated from a single-index model, meaning for some fixed $$\mu \in \mathbb{S}^{d-1}$$ and some degree-$$q$$ polynomial $$\sigma_*$$ which define $$f_*(x) = \sigma_*(\langle x, \mu \rangle)$$, we have $$y^W_i = f_*(x^W_i)$$ and $$y^a_i = f_*(x^a_i)$$.

Moreover, we assume that $$\sigma_*$$ has information exponent {% cite Arous2021 -f private_feature_learning_for_single_index_models.bib %} of $$1$$, meaning that $$\alpha^*_0 = 0$$ and $$\alpha^*_1 \neq 0$$, where $$\alpha^*_k$$ denotes the $$k$$\textsuperscript{th} Hermite coefficient of $$f_*$$. This is standard {% cite Ba2022 -f private_feature_learning_for_single_index_models.bib %}, and in our regime and without the bias unit, it would have been necessary {% cite Bietti2022 Dandi2024a -f private_feature_learning_for_single_index_models.bib %}.

**Assumption 2.** (Regime) There exist $$\varepsilon_n, \varepsilon_p > 0$$ such that $$\frac{3}{16}\varepsilon_p < \varepsilon_n$$, $$p \geq d^{\varepsilon_p}$$ and $$n \geq d^{1 + 3\varepsilon_n}$$.

**Assumption 3.** (Privacy budget)

$$
\begin{equation}
    \delta \in (0, 1), \quad \varepsilon \in \left(0, 8 \log(1/\delta)\right), \quad \frac{\varepsilon}{\sqrt{\log \frac{1}{\delta}}} = \omega\left(\frac{1}{\log d}\right).
\end{equation}
$$

**Assumption 4.** (Hyper-parameter choice)

$$
\begin{equation}
    \begin{gathered}
        \lambda = \frac{p}{d^{2\varepsilon_n}}, \eta_W = d^{\frac{3\varepsilon_n}{2}}\sqrt{p}, \eta_a = \frac{d^{\varepsilon_n}\log^2 d}{p}, C_a = d^{\varepsilon_n}\sqrt{p} \log^{q+1} d, T = d^{\varepsilon_n} \\
        \sigma_W = \sqrt{\frac{d^{1 + \varepsilon_B}}{n}} \log^{q+3} d, \sigma_a = \sqrt{\eta_a T} \frac{\sqrt{8\log(1/\delta)}}{\varepsilon}.
    \end{gathered}
\end{equation}
$$

Now we can formally state the main result.

**Theorem 1.** Consider the two-layer neural network in $$\eqref{eq:model}$$ with input dimension $$d$$ and number of features $$p$$, training over two datasets as described in Algorithm $$\eqref{alg:training procedure}$$. Let $$f_{\text{NN}} = a_T^\intercal \sigma(W_1^\intercal x + b)$$ be the resulting neural network after training, then $$f_\text{NN}$$ is $$(\varepsilon, \delta)$$-differentially private, and moreover, we have

$$
\begin{equation}
    \begin{gathered}
        \mathcal{R}(f_{\text{NN}}) = \E_x[(f_{\text{NN}}(x) - f_*(x))^2] \\
        = O\left(d^{-\frac{\varepsilon_p}{8q}} + \frac{d^{1 + \varepsilon_n}}{n} \frac{\log(1.25/\delta)}{\varepsilon^2}\right) = o(1).
    \end{gathered}
\end{equation}
$$

with probability at least $$1 - \exp(-c \log^2 d)$$ where $$c$$ is an absolute constant.

## 2.3. Related works.
### 2.3.1. Feature learning.

Once the behaviour of random feature models is characterised {% cite Hu2023 Bosch2023 -f private_feature_learning_for_single_index_models.bib %} and it is clearly not capable of feature learning {% cite Yehudai2019 Refinetti2021 Yang2020 -f private_feature_learning_for_single_index_models.bib %}, two lines of research emerge. One is performing pertubative corrections to the large-with lazy regime {% cite Hanin2020 Dyer2020 Naveh2021 Naveh2020 -f private_feature_learning_for_single_index_models.bib %}.

This work belongs to another, where feature learning comes from one non-pertubative step of gradient descent. In particular, Ba et al. {% cite Ba2022 -f private_feature_learning_for_single_index_models.bib %} showed that for a slightly different model, $$f(x) = \frac{1}{\sqrt{p}} a^\intercal \sigma(W^\intercal x + b)$$, the first gradient step of first layer approximates linear term of target function. Moniri et al. {% cite Moniri2024 -f private_feature_learning_for_single_index_models.bib %} extended this analysis to show that with appropriately chosen step size $$\eta_W$$, one has $$\ell$$ in the spectrum of feature matrix corresponding to terms of degree at most $$\ell$$ of the target. Different from this work, in their setting, training of second layer was done via ridge regression. In another direction, Cui et al. {% cite Cui2024 -f private_feature_learning_for_single_index_models.bib %} characterised shallow neural networks after one step of gradient descent on the first layer, showed that they are equivalent to spiked random feature models, similar to Gaussian equivalence property of random feature models themselves.

### 2.3.2. Single- and multi-index models.

Due to its analytical tractability, single- and multi-index models have been subject of extensive research, for which Bruan and Hsu provided a survey {% cite Bruna2025 -f private_feature_learning_for_single_index_models.bib %}. In particular, Ba et al. {% cite Ba2023 -f private_feature_learning_for_single_index_models.bib %} considered one step of gradient descent for anisotropic single-index models. Based on this method, Oko et al. extended the analysis to sums-of-indices models {% cite Oko2024 -f private_feature_learning_for_single_index_models.bib %}, with multiple steps of gradient descent for the first layer. More generally, Damian et al. carried similar analysis for multi-index models {% cite Damian2022 -f private_feature_learning_for_single_index_models.bib %}. Dandi et al. {% cite Dandi2024 -f private_feature_learning_for_single_index_models.bib %} extended it further, showed the precise sample complexity, and hinted at a hierarchy of functions demanding increasingly large batch-size to learn. 

Beyond one step, similar analysis was carried out with gradient flow {% cite Mousavi-Hosseini2023 Bietti2022 Bietti2023 -f private_feature_learning_for_single_index_models.bib %}. As a culmination of much prior effort, for shallow neural networks without bias term, Dandi et al. characterised the classes of target learnable within finite-time of gradient-based methods as well as their sample complexity {% cite Dandi2024a -f private_feature_learning_for_single_index_models.bib %}.

### 2.3.3. Learning with differential privacy.

Whilst much has been discovered about feature learning and indexed models, much less is known about learning with differential privacy. For a long time, it was open whether gradient descent could be done differentially private without hurting the performance, until the breakthrough by Abadi et al. {% cite Abadi2016 -f private_feature_learning_for_single_index_models.bib %}. For example, some attempts led to too little privacy budget {% cite Shokri2015 -f private_feature_learning_for_single_index_models.bib %}, proven to be ineffective {% cite Hitaj2017 -f private_feature_learning_for_single_index_models.bib %}.

Likewise, for a long time, learning with differential privacy remained open. There were hints of a privacy-performance trade-off {% cite Chaudhuri2009 Chaudhuri2011 -f private_feature_learning_for_single_index_models.bib %}, which prompts a line of research focusing on low-dimensional subspaces via ad-hoc techniques {% cite Kairouz2020 Zhou2021 Yu2021 -f private_feature_learning_for_single_index_models.bib %}. Nevertheless in practice, even with only standard differentially private variants of gradient descents, one still achieves good performance {% cite Mehta2022 De2022 -f private_feature_learning_for_single_index_models.bib %}, and seemingly avoids the trade-off altogether.

Recently and independently, Bombari and Mondelli {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %}, and Brown et al. {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %} showed that it is possible to achieve differential privacy at no loss of performance by using overparameterised models, in random feature models and in linear regression respectively. The methodology is similar: one sets appropriate clipping constants and noise level, then shows with high probability, the clipping does not happen and the noise added is small enough not to affect the performance.

Interestingly, each team used different strategy: Bombari and Mondelli took $$\eta_a$$ small enough to approximate training trajectory by a differential equation, now stochastic due to Gaussian noise added, then bounded the solution by a leave-one-out argument. On the other hand, Brown et al. assumed that the target is linear, and used these hidden parameters to bound the trajectory directly via an induction argument. We shall borrow from both of these techniques to analyse second layer.

## 2.4. Proof strategy.

Our analysis of first layer follows closely that of Ba et al. {% cite Ba2023 -f private_feature_learning_for_single_index_models.bib %}, modified for isotropic data, comprising three steps.
- First we show that the first gradient step admits a simple rank-1 approximation similar to that by Ba et al. {% cite Ba2023 -f private_feature_learning_for_single_index_models.bib %} and Ba et al. {% cite Ba2022 -f private_feature_learning_for_single_index_models.bib %}, and further analysed by Sonthalia et al. {% cite Sonthalia2025 -f private_feature_learning_for_single_index_models.bib %}.
- Then, using an universal approximation theorem similar to that of Damian et al. {% cite Damian2022 Ba2023 Oko2024 -f private_feature_learning_for_single_index_models.bib %} but proven for $$\tanh$$, we construct a "certificate", that is, a choice of second layer with both small test loss and small magnitude.
- Finally, by a simplified argument similar to Ba et al. {% cite Ba2023 -f private_feature_learning_for_single_index_models.bib %}, we show that we would also achieve small test loss if we trained the second layer via ridge regression _non-privately_.
The analysis of second layer largely follows induction strategy of Brown et al. {% cite Brown2024 -f private_feature_learning_for_single_index_models.bib %}, but since the underlying model is not linear, we will replace part of the argument with a leave-one-out argument similar to that of Bombari and Mondelli {% cite Bombari2025 -f private_feature_learning_for_single_index_models.bib %}.

Section 3 is devoted to the first two steps of first layer's analysis; Section 4 finishes the rest of the argument.

## 2.5. Notations.

In what follows, all algorithms are in the natural basis. We denote by $$\| \cdot \|$$ the $$\ell_2$$ norm for vector, and the induced operator norm for matrices. For sub-Gaussian random variables, we denote by $$\|\cdot\|_{\psi_2}$$ the Gaussian norm {% cite Vershynin2018 -f private_feature_learning_for_single_index_models.bib %}.

Given to quantities $$a, b$$, we denote $$a \lesssim b$$ or $$b \gtrsim a$$ to mean there exists a hidden absolute constant $$C > 0$$ such that $$a \leq C b$$, and $$a \asymp b$$ to mean $$a \lesssim b$$ and $$b \lesssim a$$. All complexity notations are understood to be taken with the size of dataset $$n$$, input dimension $$d$$, and number of hidden neurons $$p$$ large enough.

The notation $$O_\P(\cdot), o_\P(\cdot), \lesssim_\P, \asymp_\P$$ are to be understood as $$O(\cdot), o(\cdot), \lesssim, \asymp$$ with probability at least $$1 - \exp(-c\log^2 d)$$ as $$d$$ (and thus $$n$$ and $$p$$) tends to infinity.

We write $$\he_n$$ and $$\h_n$$ for probabilist's and physicist's $$n$$<sup>th</sup> Hermite polynomial, given by

$$
\begin{equation}
	\he_n = (-1)^n e^{\frac{x^2}{2}} \frac{d^n}{dx^n} e^{-\frac{x^2}{2}}, \quad \h_n = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}.
\end{equation}
$$

Let $$\alpha^b_n$$ be the $$n$$<sup>th</sup> (probabilist) Hermite coefficient of $$\sigma(\cdot + b)$$, and $$\alpha_n = \alpha^0_n$$. Likewise, let $$\alpha^*_n$$ be the $$n$$<sup>th</sup> (probabilist) Hermite coefficient of $$f_*$$. Finally, we denote $$w_i \in \R^d$$ (resp. $$w^0_i, w^1_i$$) be the $$i$$<sup>th</sup> column of $$W$$ (resp. $$W_0, W_1$$).

## References

{% bibliography --cited -f private_feature_learning_for_single_index_models.bib %}