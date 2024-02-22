# saem-zibr: Estimation of the ZIBR model for microbiome data with the SAEM algorithm

## Definition

The zero-inflated beta regression (ZIBR) model was proposed by Chen and Li (2016), a two-stage mixed effects model that allows the inclusion of covariates both to
explain the presence or not of a certain bacterial taxon and, in case of presence, the influence of these covariates in the relative abundance of the taxon.

Its defintion is as follows: 

Let us define $Y_{it}$ as the relative abundance of a bacterial taxon in the individual $i$ at time $t$, $1\leq i\leq N$, $1\leq t\leq T_i$. The model assumes that
$$Y_{it}\leadsto
    \begin{cases}
    0&\mbox{with prob. } 1-p_{it},\\
    Beta(u_{it}\phi,(1-u_{it})\phi)&\mbox{with prob. } p_{it}.\\
    \end{cases}    
$$

with $0\leq Y_{it}<1$, $\phi>0$ and $0<u_{it},p_{it}<1$. These two last components are determined by

$$\begin{split}
    \log{\left(\frac{p_{it}}{1-p_{it}}\right)}&=a_i+X_{it}^T\alpha,\\
    \log{\left(\frac{u_{it}}{1-u_{it}}\right)}&=b_i+Z_{it}^T\beta,    
\end{split}$$

where $a_i$ and $b_i$ are individual specific intercepts, $\alpha$ and $\beta$ are vectors of regression coefficients and $X_{it}$ and $Z_{it}$ are covariates for each individual and time point. We further consider that $$a_i\leadsto N(a,\sigma^2_1),\hspace*{0.5cm}b_i\leadsto N(b,\sigma^2_2).$$


