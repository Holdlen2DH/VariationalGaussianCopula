# Variational Gaussian Copula Inference

We use Gaussian copulas (combined with fixed/free-form margins) as **automated inference engines** for variational approximation in generic hierarchical Bayesian models (The only two model-specific terms are the log likelihood & prior term and its derivatives). We evaluate the **peculiarities** reproduced in the univariate margins and the **posterior dependence** captured broadly across latent variables.

## Matlab code for the paper

Shaobo Han, Xuejun Liao, David B. Dunson, and Lawrence Carin, <a href="http://people.ee.duke.edu/~lcarin/VGC_AISTATS2016.pdf"> "Variational Gaussian Copula Inference"</a>, *The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)*, Cadiz, Spain, May, 2016

## Examples

#### Demo 1: Marginal Adaptation (Skew normal, Student's t, Beta, Gamma) 

```Matlab
>> demo_SkewNormal
>> demo_StudentT
>> demo_Gamma
>> demo_Beta
```
The accuracy of marginal approximation for real, positive real, and truncated [0,1] variables, 

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/margins.png" align="center" height="200" width="800"></a>


---
#### Demo 2: Bivariate Log-Normal

```Matlab
>> demo_BivariateLN
```
Approximate bivariate log-normal distributions using a bivariate Gaussian copula with (1) fixed-form log-normal distributed margins (2) free-form Bernstein polynomial based margins,

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/lognormal.png" align="center" height="300" width="800"></a>

---
#### Demo 3: Horseshoe Shrinkage

Baseline comparisons:  
1.Gibbs sampler 
2.Mean-field VB  
3.VGC-LN-full: Gaussian Copula with Log-normal margins  
4.VGC-LN-diag: Independent Copula with Log-normal margins 
5.VGC-BP-full: Gaussian Copula with Bernstein polynomial margins

```Matlab
>> demo_Horseshoe
```
---
#### Demo 4: Poisson Log-Linear Regression

MCMC sampler implemented in RJAGS:

```r
>> demo_JAGS_PoissonLogLinear
```
VGC-BP method: 

```Matlab
>> demo_VGC_PoissonLogLinear
```

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/VGC-JAGS.png" align="center" height="550" width="800"></a>

Univaraite Margins and Pairwise Posteriors (JAGS v.s. VGC-BP)

---

## Citations

If you find this code helpful, please cite the work using the following information:

    @inproceedings{VGC_2016,
      title={Variational Gaussian Copula Inference},
      author={Shaobo Han and Xuejun Liao and David B. Dunson and Lawrence Carin},
      booktitle={The 19th International Conference on Artificial Intelligence and Statistics (AISTATS)},
      year={2016},
    }
