from remez import Remez
import mpmath
import pickle

# Use high precission arithmetics
mpmath.mp.prec = 1536


def refine(R):

  score = 0
  alternand = 0
  position = 0
  squeeze = 0

  lastShift = None


  while True:
    sq_update = False
    al_update = False

    if score > 20:
      if alternand > 10:
        shift=sum((R.t-R.t0)**3)
        position = -float(mpmath.log10(abs(shift)))
        if lastShift is None:
          dC = 1e-9*R.C
        else:
          dC = -dC/(shift - lastShift)*shift

        dC = min(dC, 0.05*R.C)
        dC = max(dC, -0.05*R.C)
        lastShift = shift
        R.C += dC

        squeeze = -float(mpmath.log10(abs(dC)))
        sq_update = True

      alternand = -float(mpmath.log10(R.updateAlternands(1)))
      al_update = True

    score = -float(mpmath.log10( R.update(1) ))

    print(f'N:{N} score:{score:7.2f} alternand:{alternand:7.2f}{"*" if al_update else " "} shift:{position:7.2f} squeeze:{squeeze:7.2f}{f" dC={float(dC):.3e}" if sq_update else " "}')


    if score > 250 and alternand > 250 and squeeze > 10:
       break

  return R


if __name__ == '__main__':
  for N in range(4,65,4):
    R = Remez(N)
    R = refine(R)
    pickle.dump(R, open(f'result/R{R.N:03}.dat', 'wb'))
