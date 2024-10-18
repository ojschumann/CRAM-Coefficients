from numpy import *
import scipy
import scipy.linalg
import scipy.optimize
from numpy.polynomial import Polynomial
import mpmath
from mpmath import mpf


''' Calculation of coefficients for Chebyshev Rational Approximation Method (CRAM)

    Based on the two papers:

    [1] Carpenter, A.J., Ruttan, A., Varga, R.S. (1984). Extended numerical computations on the “1/9” conjecture in rational approximation theory.
        In: Graves-Morris, P.R., Saff, E.B., Varga, R.S. (eds) Rational Approximation and Interpolation.
        Lecture Notes in Mathematics, vol 1105. Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0072427

    [2] Maria, P. (2016). Higher-Order Chebyshev Rational Approximation Method and Application to Burnup Equations.
        Nuclear Science and Engineering, 182(3), 297–318. https://doi.org/10.13182/NSE15-26

'''


class Remez:
  def __init__(self, N, sc=0.9):
    self.N = N
    M = 2*N+1

    t = lambda j: -cos(-pi*j/M) * (1 - 3.36*(j/M)**2*(1-j/M)**2)  # Nearly roots of chebyshev poly [1] eq. 3.5
    self.t = [mpf(t(j)) for j in range(M+1)]
    self.t0 = array(self.t) # save initial state for later

    self.rho = mpf(10)**-N # initial guess for final error. Sussecc of fit STRONGLY depends on this parameter!


    # Guess for initial polynom. Based on the observation [2] eq. 13 and 14
    # sc parameter just choosen arbitrarely
    self.p = [1-self.rho]
    self.q = [mpf(1)]

    for _ in range(N):
      self.p.append(mpf( sc*self.p[-1]))
      self.q.append(mpf(-sc*self.q[-1]))

    # Use the new numpy Polynomial class
    self.p = Polynomial(self.p)
    self.q = Polynomial(self.q)

    # Update the p_0 and p_N
    self.updateP()

    # Squeeze Parameter, to have the alternants more evenly spaced in [-1:1] see [1] eq. 3.2
    self.C = N / (1.41492705366909 - 0.649837326265066 / N**0.974869345680685) # Empirical formula


  def phi(self, t):
    '''Mapping of [-1:1] to negative real axis'''
    return self.C * (t+1) / (t-1)


  def Dphi(self, t):
    '''first derivative of phi(t)'''
    return -2 * self.C /(t-1)**2


  def D2phi(self, t):
    '''second derivative phi(t)'''
    return 4 * self.C / (t-1)**3


  def score(self, t):
    '''The score function'''
    if t == 1: # x = -oo
      return self.p.coef[-1] / self.q.coef[-1]

    x = self.phi(t)
    return self.p(x) / self.q(x) - mpmath.exp(x)


  def Dscore(self, t):
    '''first derivative of score(t)'''
    x = self.phi(t)

    P = self.p
    dP = P.deriv()

    Q = self.q
    dQ = Q.deriv()
    return (dP(x) / Q(x) - dQ(x) * P(x) / Q(x)**2 - mpmath.exp(x)) * self.Dphi(t)


  def D2score(self, t):
    '''first derivative of score(t)'''
    x = self.phi(t)

    P = self.p
    dP = P.deriv()
    d2P = dP.deriv()

    Q = self.q
    dQ = Q.deriv()
    d2Q = dQ.deriv()

    val  = ( dP(x) / Q(x)                           - dQ(x) *  P(x) / Q(x)**2                                                        - mpmath.exp(x)) * self.D2phi(t)
    val += (d2P(x) / Q(x) - dP(x) * dQ(x) / Q(x)**2 - dQ(x) * dP(x) / Q(x)**2 - d2Q(x) *P(x) / Q(x)**2 + 2*dQ(x)**2 * P(x) / Q(x)**3 - mpmath.exp(x)) * self.Dphi(t)**2

    return val



  def updateP(self):
    '''Set the depending coefficients in the polynomial p'''
    # see [2] eq. 10 and 11
    self.p.coef[0] = 1-self.rho
    self.p.coef[-1] = self.rho*self.q.coef[-1]


  def calcM(self):
    '''calc the linear system for the newton iteration'''

    N = self.N
    M = N*2

    J = mpmath.matrix(M)
    b = mpmath.matrix(M, 1)

    # exp(xi) - self.p(xi)/self.q(xi) - (-1)**i*self.rho == 0    @ xi=0
    # => 1 - p_0 / q_0 - rho = 0 => p_0 = 1 - rho

    # exp(xi) - self.p(xi)/self.q(xi) - (-1)**i*self.rho == 0    @ xi -> -oo
    # => 0 - p_N / q_N + rho = 0 => p_M = rho * q_N

    # p(x) = (1-rho) + p_1*x + p_2*x**2 + ... + p_N-1 * x**N-1 + rho*q_N * x**N
    # q(x) =     1   + q_1*x + q_2*x**2 + ... + q_N-1 * x**N-1 +     q_N * x**N

    for i in range(M):
      ti = self.t[i+1]
      xi = self.phi(ti)

      S = (-1)**i
      b[i] = (mpmath.exp(xi) + S*self.rho)*self.q(xi) - self.p(xi) # the [2] equation after eq. 11

      for j in range(1, N):
        J[i, j-1] = -xi**j # dS/dp_i

      for j in range(1, N+1):
        J[i, N+j-2] =  xi**j * ( mpmath.exp(xi) + S*self.rho ) # dS/dq_i
      J[i, 2*N-2] -=  xi**j * self.rho # p_N is q_N*rho

      J[i, M-1] = S*self.q(xi) + 1 - self.q.coef[-1]*xi**N # dS/drho

    return J, b

  def update(self, s0=1.0):
    '''Updates the polynomial coefficients'''
    N = self.N

    # get the linear system
    M, b = self.calcM()

    # solve it
    s = mpmath.lu_solve(M, b)

    # and update the state
    self.p.coef[1:-1] -= s0*s[:N-1]
    self.q.coef[1:] -= s0*s[N-1:2*N-1]
    self.rho -= s0*s[2*N-1]
    self.updateP()

    # return the absolut sum of changes for convegence control
    return sum([abs(d) for d in s])


  def updateAlternands(self, s0=1.0):
    '''Update the alternands positions.
    The score should have alternating minima & maxima at the positions,
    so we search for zeros of the first derivative of the score.
    Classical newton iteration to search for the zero.'''

    # Remember total shift
    shift = 0
    # loop over all inner alternands
    for n in range(1, 2*self.N+1):
      t = self.t[n]
      # Calculate the shift
      ds = self.Dscore(t)
      d2s = self.D2score(t)
      if (d2s>0) == (n%2==0): # min/max position correct
        d = -s0 * ds / d2s
      else:
        # mofe point to middle of neighbours
        d = 0.5*(self.t[n-1] + self.t[n+1]) - self.t[n]

      # limit the shift to maximal 10% of the distance to the neigbouring position. Alternnandes need to stay sorted!!!
      d = max(d, 0.4*(self.t[n-1]-self.t[n]))
      d = min(d, 0.4*(self.t[n+1]-self.t[n]))

      # Sum absolut shift
      shift += abs(d)

      # Apply shift
      self.t[n] += d
    return shift

  def updateSqueeze(self):
    '''Update the squeeze parameter, to have mean(alternands) == initial guess'''
    dt = 10 * mean(self.t - self.t0)
    self.C += dt
    return abs(dt)
