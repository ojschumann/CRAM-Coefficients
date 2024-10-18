from remez import Remez
import glob
import mpmath
from numpy.polynomial import Polynomial
import pickle
mpmath.mp.prec = 1024


def findRoots(p, z=1j):
  '''simple root finding with newton algorithm'''
  roots = []
  while p.degree() > 0:
    dp = p.deriv()
    while True:
      dz = -p(z)/dp(z)
      z += dz
      if abs(dz) < mpmath.mpf('1e-250'):
        break

    roots.append(z)
    if abs(z.imag) < mpmath.mpf('1e-250'):
      # real root
      p = p // Polynomial((-z, 1))
    else:
      p = p // (Polynomial((-z, 1)) * Polynomial((-z.conjugate(), 1)))
  return roots


def fmt(z):
  '''nice formatin of complex numbers for the tables'''
  _fmt = lambda z: mpmath.nstr(z, 18, show_zero_exponent=True, min_fixed=1, max_fixed=0, strip_zeros=False)
  if isinstance(z, mpmath.mpf):
    return _fmt(z)
  else:
    return (" " if z.real>0 else "") + _fmt(z.real) + (" + " if z.imag >=0 else " - ") + _fmt(abs(z.imag))





def calcAlpha(theta, t1, t2):
  '''calculation of the residual. [2] eq. 20'''
  alpha = -0.5j * (theta - t1)*(theta - t2) / theta.imag
  return alpha


# Load all calculated fits and output theta and alpha.
# the thetas are just the roots of the denominator polynomial q
# alphas are generated by combining zwo roots of the nominator polinomial p with a
# conjugated root from q (e.g. one theta value)
# This combination is entirely arbirtrary, if one choice is farourable in terms of
# accurary of the resulting CRAM is unknown to me...
# I tried to generate alphas, which are small by magnitude
#
# First I sorted the thetas by absolut value of the imaginary part, because that is
# in the denominator for the the resuidual. Then I split the roots of p into real values
# and complex values. The reals are sorted by absolute value and (artificially)
# converted to complex by using the larges and smalles real value, the second largest
# and second smales value and so on. By that, no very large numbers are combined.
# the complex roots and the "artificial" real roots are sorted by absoulte value.
# here it is remebered, if the number is a real complex value or just the combination
# of two real roots. Then the alphas are calculated from this sorted list, splitting
# thet complex number to two reals if nessessary. By that, the small number is paired
# with a theta with small imaginary part and the large number is paired with
# a theta with large imaginary part.


for fn in sorted(glob.glob("result/R*.dat")):
  R = pickle.load(open(fn, 'rb')) # Load fit result

  # Print order and limiting value for x-> -oo
  print(f'Order = {R.N}')
  print(f'alpha0: {fmt(R.p.coef[-1]/R.q.coef[-1])}')

  # find the roots of poly q
  theta = findRoots(R.q)

  # sort them according to the absolut imaginary value
  theta.sort(key = lambda z: abs(z.imag))

  # Print thetas
  for t in theta:
    print(f'Theta: {fmt(t)}')


  # find the roots of poly p
  L = findRoots(R.p)

  # Split to real and complex roots (complex roots in fact two complex conjudated roots, from which only one is present in the array)
  Lc = [z for z in L if abs(z.imag) > mpmath.mpf('1e-150')]
  Lr = [z.real for z in L if abs(z.imag) < mpmath.mpf('1e-150')]

  # Check that the splitting worked
  assert (len(Lc) + len(Lr) == len(L)) # number of roots has to be conserved
  assert (2*len(Lc) + len(Lr) == 2*len(theta)) # total number of roots to match the poly degree
  assert (len(Lr) % 2 == 0) # has to be an even number of real roots

  # sort real roots according to magnitude
  Lr.sort(key = lambda x: abs(x))

  # Combine larges and smalles element to one complex value and add all complex roots. Remeber, what were the complex roots
  N = len(Lr)-1
  U = [ (Lr[n] + 1j*Lr[N-n], False) for n in range(len(Lr)//2) ] + [(z, True) for z in Lc ]

  # Sort according to magnitude
  U.sort(key = lambda key: abs(key[0]))

  # Calculate and print alphas
  for n, (z, isRealComplex) in enumerate(U):
    if isRealComplex:
      alpha = calcAlpha(theta[n], z, z.conjugate())
    else:
      alpha = calcAlpha(theta[n], z.real, z.imag)
    print(f'alpha: {fmt(alpha)}')
  print()