# Variational Gaussian Copula Inference

We use Gaussian copulas (combined with fixed/free-form margins) as **automated inference engines** for variational approximation in generic hierarchical Bayesian models (The only two model-specific terms are the log likelihood & prior term and its derivatives). We evaluate the peculiarities reproduced in the univariate margins and the posterior dependence captured broadly across latent variables.

## Matlab code for the paper

Shaobo Han, Xuejun Liao, David B. Dunson, and Lawrence Carin, <a href="http://people.ee.duke.edu/~lcarin/VGC_AISTATS2016.pdf"> "Variational Gaussian Copula Inference"</a>, *The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)*, Cadiz, Spain, May, 2016

## Examples

#### Demo 1: Flexible Margins (Skew normal, Student's t, Beta, Gamma) 
<figure>
<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/margins.png" align="center" height="200" width="800"></a>
<figcaption>Fig1. - A view of the pulpit rock in Norway.</figcaption>
</figure>



#### Demo 2: Bivariate Log-Normal

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/lognormal.png" align="center" height="300" width="800"></a>


#### Demo 3: Horseshoe Shrinkage


#### Demo 4: Poisson Log-Linear Regression

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/VGC-JAGS.png" align="center" height="550" width="800"></a>


## Citations

If you find this code helpful, please cite the work using the following information:

    @inproceedings{VGC_2016,
      title={Variational Gaussian Copula Inference},
      author={Shaobo Han and Xuejun Liao and David B. Dunson and Lawrence Carin},
      booktitle={The 19th International Conference on Artificial Intelligence and Statistics (AISTATS)},
      year={2016},
    }
