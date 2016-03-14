# Variational Gaussian Copula Inference

We use Gaussian copulas with fixed/free-form margins as automated inference engines for variational approximation in generic hierarchical Bayesian models. The only two model-specific terms are the log likelihood & prior term and its derivatives. We evaluate the peculiarities reproduced in the univariate margins and the posterior dependence captured broadly across latent variables.

# Matlab code for the paper

Shaobo Han, Xuejun Liao, David B. Dunson, and Lawrence Carin, "Variational Gaussian Copula Inference", The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016), Cadiz, Spain, May, 2016

# Details

* Demo 1: Flexible Margins (Skew normal, Student's t, Beta, Gamma) 

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/flexiblemargins.png" align="left" height="450" width="450" ></a>


* Demo 2: Bivariate Log-Normal

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/lognormal.png" align="left" height="450" width="900" ></a>


* Demo 3: Horseshoe Shrinkage

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/horseshoe.png" align="left" height="450" width="500" ></a>


* Demo 4: Poisson Log-Linear Regression

<a href="url"><img src="https://github.com/shaobohan/VariationalGaussianCopula/blob/master/figure/VGC-JAGS.png" align="left" height="600" width="900" ></a>

# Citations

If you find this code helpful, please cite the work using the following information:

    @inproceedings{VGC_2016,
      title={Variational Gaussian Copula Inference},
      author={Shaobo Han and Xuejun Liao and David B. Dunson and Lawrence Carin},
      booktitle={The 19th International Conference on Artificial Intelligence and Statistics (AISTATS)},
      year={2016},
    }
