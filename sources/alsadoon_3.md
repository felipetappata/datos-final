## 3 Monte Carlo experiment of the pure AR(1) model

For the Monte Carlo experiment, we consider the following data-generating processes. First, we assume two different options for the selection equation:

$$
\begin{array}{r}
d_{i t}^{*}=a-z_{i t}-\eta_{i}-u_{i t} \\
d_{i t}^{*}=a-0.5 d_{i t-1}+z_{i t}-\eta_{i}-u_{i t} \\
d_{i t}=1\left[d_{i t}^{*}>0\right] \tag{42}
\end{array}
$$

where $a$ is set so $p\left(d_{i t}^{*}>0\right)=0.85$ and $z_{i t} \sim N\left(0, \sigma_{z}\right)$ with $\sigma_{z}=1$. Second, the outcome of interest is generated as follows:

$$
\begin{array}{r}
y_{i t}^{*}=\left(2+\alpha_{i}+\varepsilon_{i t}\right) /(1-\rho) \text { if } t=1 \\
y_{i t}^{*}=2+\rho y_{i t-1}^{*}+\alpha_{i}+\varepsilon_{i t} \text { if } t=2, \ldots, T \\
y_{i t}=y_{i t}^{*} \text { if } d_{i t}=1 \tag{45}
\end{array}
$$

We let $\rho$ vary between $0.25,0.50$ and 0.75 . We generate all variables for $T=17$ to $T=20$ and discard the first 13 observations to minimise any problem with initial conditions. ${ }^{6}$ Finally, we assume the following structure for the errors:

$$
\begin{array}{r}
\eta_{i} \sim N\left(0, \sigma_{\eta}\right) \text { with } \sigma_{\eta}=1 \\
u_{i t} \sim N\left(0, \sigma_{u}\right) \text { with } \sigma_{u}=1 \\
\alpha_{i}=\alpha_{i}^{0}+0.5 \eta_{i}, \alpha_{i}^{0} \sim N\left(0, \sigma_{\alpha^{0}}\right) \text { with } \sigma_{\alpha^{0}}=1 \\
\varepsilon_{i t}=\varepsilon_{i t}^{0}+0.5 u_{i t}, \varepsilon_{i t}^{0} \sim N\left(0, \sigma_{\varepsilon^{0}}\right) \text { with } \sigma_{\varepsilon^{0}}=1 \tag{49}
\end{array}
$$

These assumptions imply that $\operatorname{corr}\left(\varepsilon_{i t}, u_{i t}\right)=\operatorname{corr}\left(\alpha_{i}, \eta_{i}\right)=0.5 / \sqrt{1+0.5^{2}}=0.447$.

### 3.1 Description of the experiments

For each experiment, we set the initial sample size to $N=500$ or $N=5000$, and for each $i$, we draw up to 20 time series observations, from which the initial 13 are discarded. Once selection is applied, the unbalanced panels are formed. At least three consecutive observations of the same regime are needed to form an observation of the selected panel. This implies that a large fraction of the observations do not contribute to the identification of the parameters, even with a small

[^5]degree of sample selection. For example, a 15 per cent of initial selection implies that around $1 / 3$ of the observations are lost. For each combination of the parameters we perform 500 replications.

In each case, we evaluate the performance of two well-known GMM-IV estimators: AB and system. The structure of the model makes selection of the instruments a crucial step of this simulation study. We select the instruments as follows: we use lags from $t-2$ backwards for firstdifferenced equations, although we also evaluate the performance of the estimates with a restricted set of instruments. We use the lagged first difference of the outcome as an additional instrument for the equation in levels. Although we are aware of the instrument proliferation issue analysed by Roodman (2009), it does not constitute a problem here given the reduced number of periods (a maximum of 7) remaining for estimation. ${ }^{7}$

### 3.2 Simulation results for the pure autoregressive model

### 3.2.1 The basic results

Table 1 presents results for the AR(1) model for three values of the autoregressive parameter: $0.25,0.50$ and $0.75 .{ }^{8}$ We simulate two alternative selection models (static, A, and dynamic, B), as presented in equations (40) and (41). For each combination of selection model and autoregressive parameter, we report results for both the AB and the system estimators constructed under competing assumptions about the selection process: (a) non-endogenous selection; (b) endogenous selection without correction. The initial degree of sample selection is 15 per cent, while the fraction of the sample lost is much larger (around $1 / 3$ of the observations).

Let us start reviewing the results without endogenous selection, reported in columns (1) and (2). When the initial sample (before selecting the observations) is small $(N=500)$ the bias of the AB grows with the autoregressive parameter (for both selection models, A and B ) and becomes sizable when $\rho=0.75 .^{9}$ As we increase the sample size ( $N=5000$ ), the average bias of the AB estimator is reduced substantially and only remains noticeable for $\rho=0.75$. Alternatively, the system estimator, which is also consistent in this case, shows a very small bias for $N=500$ (never exceeding one per cent), even smaller when $N=5000$. Figure 1 confirms these results with a sample size varying from $N=200$ to $N=5000$ in absence of any sort of selection (estimators labelled AB all and system all).

When endogenous sample selection is considered (see columns (3) and (4) for AB and system estimators results) but we do not include any correction for endogenous sample selection in the model, we do not detect any significant change in the bias results for the AB estimator for both selection models, even when the initial sample is small . In fact, when the initial sample is small,

[^6]the difference between the cases with and without selection is practically undetectable (see Table A1). Alternatively, the results confirm the consistency of the AB and the small bias of the system estimator for $N=5000$. In contrast, the system estimator always shows a very small bias (between 1 per cent for $\rho=0.25$ and 2.25 per cent for $\rho=0.75$ ). Note that the bias becomes more evident as the sample size grows (see Figure 1). As a sort of compensation, the standard errors always tend to be substantially smaller.

Some additional conclusions can be drawn when varying sample size (Figure 1). When $N=200$, the AB estimator shows sizable bias, which decreases as $N$ increases. The system estimator has a very small bias, however. For a given $\rho$, it remains stable (between 1 and 2.5 per cent) as $N$ increases. We can find a threshold for $N$ for each combination of parameters. Below this threshold, the average bias of the system estimator is smaller, and it is larger above it. Therefore, we may conclude that for moderate and small samples (say, below the range 1000-1500), the system estimator is highly recommended because of the likely smaller small sample bias as well as smaller variance.

Finally, as shown in Figure 2, when the individual heterogeneity components, $\alpha_{i}$ and $\eta_{i}$, are not correlated, the bias of the system estimator practically disappears (in comparison with the previous case) due to the fact that the main source of bias is the correlation between the heterogeneous components of both the outcome and selection equations (see Table A. 1 for an illustration). ${ }^{10}$ This means that in those cases in which endogenous selection is not due to individual heterogeneity, all three estimators considered are, in essence, valid options for recovering the key parameters of the outcome equation.

[^7]Table 1: Average bias in the $\mathbf{A R ( 1 )}$ model ( $T=7,500$ replications)

|  |  |  | No endogenous selection |  | Endogenous selection |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Select. |  |  | (1) | (2) | (3) | (4) |
| Model | $\rho$ |  | AB | SYS | AB | SYS |
| $N=500$ |  |  |  |  |  |  |
| A | . 25 | bias | -. 01002 | . 00302 | -. 01054 | -. 00151 |
|  |  | s.e. | . 06004 | . 04378 | . 05447 | . 04209 |
| A | . 50 | bias | -. 03265 | . 00224 | -. 02964 | -. 00700 |
|  |  | s.e. | . 08833 | . 04990 | . 07807 | . 05027 |
| A | . 75 | bias | -. 19339 | . 00757 | -. 10176 | -. 00923 |
|  |  | s.e. | . 19402 | . 05879 | . 13350 | . 06437 |
| $N=5000$ |  |  |  |  |  |  |
| A | . 25 | bias | -. 00209 | . 00017 | -. 00059 | -. 00381 |
|  |  | s.e. | . 01661 | . 01272 | . 04358 | . 02387 |
| A | . 50 | bias | -. 00461 | . 00019 | -. 00264 | -. 00868 |
|  |  | s.e. | . 02479 | . 01491 | . 02285 | . 01493 |
| A | . 75 | bias | -. 02395 | . 00112 | -. 01118 | -. 01709 |
|  |  | s.e. | . 05792 | . 01791 | . 03863 | . 01894 |
| $N=500$ |  |  |  |  |  |  |
| B | . 25 | bias | -. 01076 | . 00216 | -. 01075 | -. 00024 |
|  |  | s.e. | . 05904 | . 04293 | . 05471 | . 04219 |
| B | . 50 | bias | -. 03324 | . 00135 | -. 03022 | -. 00556 |
|  |  | s.e. | . 08802 | . 04925 | . 07843 | . 05039 |
| B | . 75 | bias | -. 18966 | . 00709 | -. 10478 | -. 00823 |
|  |  | s.e. | . 19100 | . 05824 | . 13880 | . 06358 |
| $N=5000$ |  |  |  |  |  |  |
| B | . 25 | bias | -. 00208 | $-7.75 \mathrm{e}-07$ | -. 00103 | -. 00134 |
|  |  | s.e. | . 01680 | . 01267 | . 01700 | . 01279 |
| B | . 50 | bias | -. 00459 | . 00019 | -. 00263 | -. 00649 |
|  |  | s.e. | . 02540 | . 01481 | . 02380 | . 01482 |
| B | . 75 | bias | -. 02375 | . 00105 | -. 01129 | -. 01472 |
|  |  | s.e. | . 05964 | . 01780 | . 04192 | . 01871 |

Figure 1: Average bias of the $\mathbf{A B}$ and system estimators in the full sample ( $(N x T$ observations) and the endogenously selected sample
![](https://cdn.mathpix.com/cropped/2025_05_20_75ad17dbbf7cfcced0e8g-16.jpg?height=1191&width=1637&top_left_y=501&top_left_x=255)

Notes.
AB all: AB GMM-IV estimates using the full ( $N x T$ ) sample (no selection process). system all: System GMM-IV estimates using the full ( $N x T$ ) sample (no selection process).
AB select: Uncorrected for selection AB GMM-IV estimates on the selected sample under endogenous sample selection.
system select: Uncorrected system GMM-IV estimates on the selected sample under endogenous sample selection.

Figure 2: Average bias of system estimator in the full sample ( $N x T$ observations) and the endogenously selected sample when $\alpha_{i}$ and $\eta_{i}$ are not correlated
![](https://cdn.mathpix.com/cropped/2025_05_20_75ad17dbbf7cfcced0e8g-17.jpg?height=1191&width=1637&top_left_y=501&top_left_x=255)

Notes.
system all: System GMM estimates with the full sample (no-selection).
system select: Uncorrected system GMM estimates with the selected sample under endogenous selection due to correlation of the time-varying errors.

### 3.2.2 Sensitivity analysis

In this section, we comment on various departures from the basic assumptions. We consider the following representative cases: (a) varying the longitudinal dimension of the panel; (b) increasing the percentage of selection (from 0.15 to 0.25 ); (c) increasing the ratio of the variances to $\frac{\sigma_{\alpha}^{2}}{\rho_{\varepsilon}^{2}}=2$; (d) reducing the correlation between the errors (the correlation parameter is reduced from 0.5 to 0.25 ); (e) and, finally, non-stationary time varying errors and correlation of the time-varying error components. In particular, we allow the variance of the time-varying errors in (1) and (2) to vary over time ${ }^{11}$ and we also allow the correlation coefficient between the time-varying errors in (1) and (2) to vary over time. ${ }^{12}$ We present the simulation results for $N=500$ in Table $\mathbf{2}$ and for $N=5000$

[^8]in Table 3, for three values of the autoregressive parameters: $0.25,0.50$ and 0.75 .
Our first experiment reduces the maximum longitudinal dimension of the observed panel from $T=7$ to $T=4$. Apart from the expected increase in the estimated variance and regardless of the sample size considered, the effect on the average bias of the $\mathrm{AR}(1)$ coefficient implied by this change is very small for both estimators, all values of the autoregressive parameter considered and both sample sizes.

Increasing the degree of sample selection from 0.15 to 0.25 increases average bias of the autoregressive coefficient very mildly. In addition, it increases its variance due to the significant reduction in the number of observations (the average number of observations is reduced around 30 per cent).

The increase of the ratio of the variance of the individual heterogeneous component to the variance of the time-variant component of the outcome equation does not have an important effect on the average bias of the estimated parameters, either in the AB or the system GMM estimator. We can observe in Tables $\mathbf{2}$ and $\mathbf{3}$ that the effect is smaller when the sample size is larger.

We also consider a reduction in the correlation parameter of the errors. In particular, we assume the following structure for the errors: $\varepsilon_{i t}=\varepsilon_{i t}^{0}+0.25 u_{i t}$ and $\alpha_{i}=\alpha_{i}^{0}+0.25 \eta_{i}$, which implies a correlation coefficient of $0.2425\left(=0.25 / \sqrt{1+0.25^{2}}\right)$ in either case. As can be easily detected by comparing the results reported in Tables 1, $\mathbf{2}$ and 3, this change reduces the average bias of the estimators for both sample sizes and all autoregressive parameters. Finally, when we introduce in the model non-stationary time-varying error components and we allow the correlation between the errors in the outcome and the selection equation to vary over time, the average bias of the estimated parameters is significantly reduced, specially important when the initial sample size is small.

In sum, these sensitivity exercises confirms the main lessons we can draw from the analysis: (a) the AB (or the AH ) estimator is moderately biased when $N$ is small or moderate, and unbiased when $N$ is large. The system GMM estimator is always moderately biased. All these results imply that the system estimator is especially recommended when the sample size is small or even moderate (below one or one and a half thousand individuals) and less "important" when the sample size is larger.

Table 2: Average bias in the AR(1) model. Sensibility analysis for small $N$

| Model |  | $\rho=0.25$ |  | $\rho=0.5$ |  | $\rho=0.75$ |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  | AB | SYS | AB | SYS | AB | SYS |
| Experiment I: Very short $T(T=4)$ |  |  |  |  |  |  |  |
| A | bias | -. 00351 | -. 00096 | -. 01775 | -. 01349 | -. 08907 | -. 03272 |
|  | s.e. | . 13394 | . 07796 | . 20964 | . 09687 | . 39245 | . 13551 |
| B | bias | -. 00094 | . 00299 | -. 01147 | -. 00949 | -. 07759 | -. 02922 |
|  | s.e. | . 13261 | . 07737 | . 21371 | . 09551 | . 44549 | . 13314 |
|  | Experiment II: More sample selection (25\%) |  |  |  |  |  |  |  |
| A | bias | -. 01878 | -. 00135 | -. 04518 | -. 00622 | -. 13019 | -. 00601 |
|  | s.e. | . 06665 | . 05032 | . 09556 | . 06172 | . 15362 | . 07985 |
| B | bias | -. 01941 | . 00122 | -. 04817 | -. 00377 | -. 14386 | -. 00476 |
|  | s.e. | . 06811 | . 05026 | . 09801 | . 06047 | . 16296 | . 07769 |
|  | Experiment III: Increasing the ratio of variances: $\sigma_{\eta} / \sigma_{\epsilon}=2$ |  |  |  |  |  |  |  |
| A | bias | -. 01095 | . 00071 | -. 03149 | -. 00263 | -. 11253 | . 00535 |
|  | s.e. | . 05608 | . 04446 | . 08088 | . 05430 | . 14306 | . 07006 |
| B | bias | -. 01121 | . 00162 | -. 03211 | -. 00151 | -. 11640 | . 00565 |
|  | s.e. | . 05596 | . 04433 | . 0812 | . 05433 | . 14816 | . 06920 |
| Experiment IV: Reducing the correlation of the errors: $\rho=0.25$ |  |  |  |  |  |  |  |
| A | bias | -. 01018 | . 00109 | -. 03182 | -. 00069 | -. 14046 | . 00340 |
|  | s.e. | . 05720 | . 04293 | . 08392 | . 04999 | . 15923 | . 06081 |
| B | bias | -. 01077 | . 00077 | -. 03197 | -. 00050 | -. 14149 | . 00416 |
|  | s.e. | . 05673 | . 04243 | . 08376 | . 04961 | . 15990 | . 06031 |
| Experiment V: Non-stationary time-varying error components |  |  |  |  |  |  |  |
| A | bias | -. 01263 | -. 00333 | -. 02599 | -. 00784 | -. 07935 | -. 00925 |
|  | s.e. | . 05537 | . 03968 | . 07605 | . 04742 | . 12415 | . 05987 |
| B | bias | -. 01246 | -. 00179 | -. 02854 | -. 00593 | -. 07977 | -. 00781 |
|  | s.e. | . 05354 | . 03962 | . 07549 | . 04438 | . 12231 | . 05692 |

Notes.

1. $T=7$, except in experiment I.
2. $N=500$.
3. Number of replications: 500 .

Table 3: Average bias in the AR(1) model. Sensibility analysis for large $N$

| Model |  | $\rho=0.25$ |  | $\rho=0.5$ |  | $\rho=0.75$ |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  | AB | SYS | AB | SYS | AB | SYS |
| Experiment I: Very short $T(T=4)$ |  |  |  |  |  |  |  |
| A | bias | -. 00059 | -. 00381 | -. 00118 | -. 01070 | -. 00258 | -. 02123 |
|  | s.e. | . 04358 | . 02387 | . 06645 | . 02801 | . 11215 | . 03444 |
| B | bias | -. 00006 | -. 00191 | -. 00106 | -. 00858 | -. 00413 | -. 01891 |
|  | s.e. | . 04311 | . 02393 | . 06637 | . 02797 | . 11311 | . 03412 |
|  | Experiment II: More sample selection (25\%) |  |  |  |  |  |  |  |
| A | bias | -. 00226 | -. 00402 | -. 00462 | -. 01033 | -. 01433 | -. 01962 |
|  | s.e. | . 02124 | . 01549 | . 02968 | . 01774 | . 04639 | . 02241 |
| B | bias | -. 00218 | -. 00151 | -. 00463 | -. 00724 | -. 01566 | -. 01658 |
|  | s.e. | . 02158 | . 01572 | . 03027 | . 01803 | . 04780 | . 02257 |
|  | Experiment III: Increasing the ratio of variances: $\sigma_{\eta} / \sigma_{\epsilon}=2$ |  |  |  |  |  |  |  |
| A | bias | -. 00123 | -. 00248 | -. 00295 | -. 00721 | -. 01326 | -. 01454 |
|  | s.e. | . 01720 | . 01351 | . 02400 | . 01608 | . 04137 | . 02094 |
| B | bias | -. 00108 | -. 00099 | -. 00288 | -. 00539 | -. 01323 | -. 01253 |
|  | s.e. | . 01736 | . 01348 | . 02493 | . 01606 | . 04494 | . 02076 |
| Experiment IV: Reducing the correlation of the errors: $\rho=0.25$ |  |  |  |  |  |  |  |
| A | bias | -. 00165 | -. 00078 | -. 00373 | -. 00217 | -. 01594 | -. 00370 |
|  | s.e. | . 01653 | . 01275 | . 02423 | . 01506 | . 04627 | . 01838 |
| B | bias | -. 00164 | -. 00030 | -. 00368 | -. 00159 | -. 01615 | -. 00299 |
|  | s.e. | . 01679 | . 01264 | . 02489 | . 01499 | . 04916 | . 01812 |
| Experiment V: Non-stationary time-varying error components |  |  |  |  |  |  |  |
| A | bias | -. 00219 | -. 00246 | -. 00142 | -. 00586 | -. 01400 | -. 01263 |
|  | s.e. | . 01923 | . 01288 | . 03185 | . 01615 | . 05631 | . 02099 |
| B | bias | -. 00197 | -. 00089 | -. 00306 | -. 00484 | -. 01688 | -. 01152 |
|  | s.e. | . 01975 | . 012491 | . 02874 | . 01535 | . 05961 | . 02037 |

Notes.

1. $T=7$, except in experiment I.
2. $N=5000$.
3. Number of replications: 500 .