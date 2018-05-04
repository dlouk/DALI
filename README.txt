Dimension-Adaptive Leja Interpolation (DALI)
(alternatively: DArmstadt's Leja Interpolation) 
--------------------------------------------------------------------------------

Development/maintenance: Dimitrios Loukrezis (loukrezis@temf.tu-darmstadt.de)
--------------------------------------------------------------------------------

DALI is a Python software for multivariate approximation with a 
dimension-adaptive stochastic collocation algorithm, based on univariate Leja 
interpolation rules. The software has been developed during my PhD research at 
the Institute for Theory of Electromagnetic Fields (TEMF) of TU Darmstadt, under 
the supervision of Prof. Dr.-Ing. Herbert De Gersem (TU Darmstadt) and 
Jun.-Prof. Dr.-Ing. Ulrich Roemer (TU Braunschweig).
--------------------------------------------------------------------------------

Theory papers:
- "High-Dimensional Adaptive Sparse Polynomial Interpolation and Applications 
to Parametric PDEs", Chkifa, Cohen, and Schwab.
- "Adaptive Leja sparse grid constructions for stochastic collocation and 
high-dimensional approximation", Narayan and Jakeman.
--------------------------------------------------------------------------------

Implementation papers:
- "Numerical Comparison of Leja and Clenshaw-Curtis Dimension-Adaptive 
Collocation for Stochastic Parametric Electromagnetic Field Problems", 
Loukrezis, Roemer and De Gersem
- "Uncertainty Quantification for an Optical Grating Coupler with an Adjoint 
Error-Based Leja Adaptive Collocation Method", Georg, Loukrezis, Roemer and 
Schoeps.
--------------------------------------------------------------------------------

The 1D Leja points are computed using the Chaospy Python toolbox.
- https://github.com/jonathf/chaospy 
- "Chaospy: An open source tool for designing methods of uncertainty 
quantification", Feinberg and Langtangen
--------------------------------------------------------------------------------
