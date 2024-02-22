# saem-zibr: Estimation of the ZIBR model for microbiome data with the SAEM algorithm

## Definition

The zero-inflated beta regression (ZIBR) model was proposed by Chen and Li (2016), a two-stage mixed effects model that allows the inclusion of covariates both to
explain the presence or not of a certain bacterial taxon and, in case of presence, the influence of these covariates in the relative abundance of the taxon.

Its defintion is as follows: let us define $Y_{it}$ as the relative abundance of a bacterial taxon in the individual $i$ at time $t$, $1\leq i\leq N$, $1\leq t\leq T_i$. The model assumes that

```math
Y_{it} \leadsto
    \begin{cases}
    0&\mbox{with prob. } 1-p_{it},\\
    Beta(u_{it}\phi,(1-u_{it})\phi)&\mbox{with prob. } p_{it}.
    \end{cases}    
```

with $0 \leq Y_{it} < 1$, $\phi>0$ and $0 < u_{it}, p_{it}<1$. These two last components are determined by

$$\begin{split}
    \log{\left(\frac{p_{it}}{1-p_{it}}\right)}&=a_i+X_{it}^T\alpha,\\
    \log{\left(\frac{u_{it}}{1-u_{it}}\right)}&=b_i+Z_{it}^T\beta,    
\end{split}$$

where $a_i$ and $b_i$ are individual specific intercepts, $\alpha$ and $\beta$ are vectors of regression coefficients and $X_{it}$ and $Z_{it}$ are covariates for each individual and time point. We further consider that $$a_i\leadsto N(a,\sigma^2_1),\ b_i\leadsto N(b,\sigma^2_2).$$

The main difference that we address is the estimation method. Chen and Li use Gauss-Hermite quadrature to approximate the log-likelihood and optimize to find the estimators, and is implemented in the package [**ZIBR**](https://github.com/PennChopMicrobiomeProgram/ZIBR). We propose the usage of the *Stochastic Approximation Expectation Maximization (SAEM) algorithm* (Delyon et al., 1999) that is a useful tool in other complex mixed effects models.

## Contents

This repository has two files: **saem-zibr.R**, which contains the syntaxis for auxiliar functions and the main function *saem_zibr* that carries out the SAEM estimation of the ZIBR model given the data in a specific format; and **saem-results.R**, which shows the results of the application of the SAEM estimation in two real datasets from Lee et al. (2015) and Romero et al. (2014) that will be published in a future paper.

## Usage of saem_zibr

```r
saem_zibr<-function(Y,X=NULL,Z=NULL,index,v0,a0,b0,seed,iter,ncad=5,a.fix=NULL,b.fix=NULL)
```

- *Y:* a vector that contains the relative abundance of a bacterial taxon.
- *X:* a matrix or data frame that contains the covariates defined for the logistic part of the model $p_{it}$. Must have the same number of rows as Y. Default is NULL, which assumes only a random intercept in the corresponding part of the model.
- *Z:* a matrix or data frame that contains the covariates defined for the beta part of the model $u_{it}$. Must have the same number of rows as Y. Default is NULL, which assumes only a random intercept in the corresponding part of the model.
- *index:* a vector with labels that identify the different individuals in the data. Must have the same number of rows as Y, X and Z. Can be factor, character or integer.
-  *v0:* initial value for $\phi$. Must be a number.
-  *a0:* initial value for the coefficient vector associated to $p_{it}$. If X is NULL, is a number. If X has $m$ columns, a0 must have $m+1$ components, and the first is for the intercept.
-  *b0:* initial value for the coefficient vector associated to $u_{it}$. If Z is NULL, is a number. If Z has $m$ columns, b0 must have $m+1$ components, and the first is for the intercept.
-  *seed:* for the replicability.
-  *iter:* number of iterations for the SAEM algorithm.
-  *ncad:* number of Markov chains to use in the SAEM algorithm. Default is 5.
-  *a.fix:*  vector of ones and zeros, same length as a0. Indicates if the effect associated to a column in X is random (0) or fixed (1). Default is NULL, which takes only a random intercept and all the columns of X with fixed effects (if there is any).
- *b.fix:*  vector of ones and zeros, same length as b0. Indicates if the effect associated to a column in Z is random (0) or fixed (1). Default is NULL, which takes only a random intercept and all the columns of Z with fixed effects (if there is any).

The application of *saem_zibr* returns an object of the class *SAEM_ZIBR_result*, with print an plot methods.

## References

Chen, E.Z., Li, H.: A two-part mixed-effects model for analyzing longitudinal microbiome compositional data. 
Bioinformatics 32(17), 2611–2617 (2016)

Delyon, B., Lavielle, M., Moulines, E.: Convergence of a stochastic approximation version of
the EM algorithm. Annals of Statistics, 94–128 (1999)

Lewis, J.D., Chen, E.Z., Baldassano, R.N., Otley, A.R., Griffiths, A.M., Lee, D., Bittinger, K.,
Bailey, A., Friedman, E.S., Hoffmann, C., et al.: Inflammation, antibiotics, and diet as environmental stressors of 
the gut microbiome in pediatric Crohn’s disease. Cell Host & Microbe 18(4), 489–500 (2015)

Romero, R., Hassan, S.S., Gajer, P., Tarca, A.L., Fadrosh, D.W., Nikita, L., Galuppi, M., Lamont, R.F., 
Chaemsaithong, P., Miranda, J., et al.: The composition and stability of the vaginal microbiota of normal pregnant 
women is different from that of non-pregnant women. Microbiome 2(1), 1–19 (2014)

