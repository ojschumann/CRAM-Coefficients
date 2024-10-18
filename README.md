# CRAM-Coefficients

Generation of Coefficients for Chebyshev Rational Approximation Method (CRAM) for burnup applications

The code is based on the two papers:

    [1] Carpenter, A.J., Ruttan, A., Varga, R.S. (1984). Extended numerical computations on the “1/9” conjecture in rational approximation theory.
        In: Graves-Morris, P.R., Saff, E.B., Varga, R.S. (eds) Rational Approximation and Interpolation.
        Lecture Notes in Mathematics, vol 1105. Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0072427

    [2] Pusa, Maria (2016). Higher-Order Chebyshev Rational Approximation Method and Application to Burnup Equations.
        Nuclear Science and Engineering, 182(3), 297–318. https://doi.org/10.13182/NSE15-26


The fitRational.py script fits a rational function of order N to exp(-x) and saves it as a python pickle object. The calcCRAMCoefficients.py generates a table
like those in M.Pusa [2]. Rationale to develope this script was a flaw in the data of [2] for order 24. It was found to be a type in the coefficient theta_5. In [2] it is 
given as -7.15331... but the actual value has to be -7.15531...

All math is done with the mpmath multiprecision library. For order 64 a precission of > 1024 bits is needed. No check was done, if results depend on the used precission.
Runtime for order 64 was in about 10 minutes.

The results might be used to generate source code for CRAM coefficients in nuclear burnup codes.
