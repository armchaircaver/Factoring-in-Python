from math import gcd
from time import time,perf_counter
from operator import mul
import random
import itertools
import os
import subprocess
import csv
from functools import reduce
import binascii

verbose=False

try:
  import pyecm
  ECM_available=True
except:
  ECM_available=False
  

try:
  from factorsdll import factoriseCPP, allfactorCPP, isprimeCPP, TotientCPP
  CPP_available=True
except:
  CPP_available=False

def setCPP( x ):
  global CPP_available
  CPP_available = x
  

YAFU_available = os.path.isfile("yafu-x64.exe")

# allows switching on and off of yafu option
def setYAFU( x ):
  global YAFU_available
  YAFU_available = x

# allows switching on and off of CPP option
def setCPP( x ):
  global CPP_available
  CPP_available = x


try:  
  from gmpy2 import mpz
  import gmpy2
  gmpy2_available=True
  if verbose: print("gmpy2 available")
except:
  gmpy2_available=False
  if verbose: print("gmpy2 not available")
  
# assume factordb is available, until it is explicitly turned off
# with setfactordb or fails when it is attempted to be used

factordb_available=False
"""  
try:
  from factordb.factordb import FactorDB
  f = FactorDB(6)
  f.connect()
  facs = f.get_factor_list()
  factordb_available=True
  if verbose: print("factorDB available")
except:
  factordb_available=False
  if verbose: print("factorDB not available")
"""

def setfactordb( x ):
  global factordb_available
  factordb_available = x



other_methods_avaliable = factordb_available or YAFU_available or ECM_available

MAXLOOKUP=1000000
MAXSMALLPRIME=1200


# prime_or_factor[n]==0 indicates n is a prime,
# otherwise prime_or_factor[n] holds a prime factor of n
start=time()
prime_or_factor=[0]*(MAXLOOKUP+1)
max_p = int(MAXLOOKUP**0.5)+1
prime_or_factor[4::2]=[2]*len(prime_or_factor[4::2])
prime_or_factor[9::6]=[3]*len(prime_or_factor[9::6])
p=5
while (p <= max_p):
  if prime_or_factor[p]==0:
    prime_or_factor[p*p::2*p]=[p]*len(prime_or_factor[p*p::2*p])
  if p%6==5:
    p+=2
  else:
    p+=4

smallprimes = [p for p in range(2,MAXSMALLPRIME) if prime_or_factor[p]==0]
#--------------------------------------------------------------------------


def yafu_factorise(n, timeout):
  """
  run yafu to factorise a number, and read the output from -of <file>
  into a list.
  if the number is prime, -of <file> option doesn't appear to be created
  """

  filename = "yafu_out.log"
  try:
      os.remove(filename)
  except OSError:
      pass

  cmd = "yafu-x64 ""factor(" + str(n) + ")"" -of "+ filename 
  #print cmd

  # https://stackoverflow.com/questions/7006238/how-do-i-hide-the-console-when-i-use-os-system-or-subprocess-call
  with open(os.devnull, 'wb') as shutup:
    p = subprocess.call(cmd, stdout=shutup, stderr=shutup,
                        creationflags=subprocess.CREATE_NO_WINDOW,
                        shell=False, timeout=timeout)

  if os.path.exists(filename):
    factors = []
    with open(filename) as csvfile:
      yafureader = csv.reader(csvfile, delimiter='/', quotechar='|')
      for row in yafureader:
        for pp in row[1:]:
          caret = pp.find('^')
          if caret > 0 :
            factors += [int(pp[:caret])]*int(pp[caret+1:])
          else:
            factors += [int(pp)]
    v=1
    for f in factors:
      v*=f
    if v==n:  
      return factors
    else:
      if verbose: print("not all factors found, calling yafu again")
      return factors + yafu_factorise(n//v, timeout)

  else:
    return [n]
#--------------------------------------------------------------------------


def factordb_factorise(n):
  global factordb_available
  try:
    from factordb.factordb import FactorDB
    f = FactorDB(n)
    f.connect()
    st = f.get_status()
    if st != 'FF':
      if verbose: print("factorDB failed for n={}, status={}".format(n,st))
      if verbose: print("partial factorisation", f.get_factor_list())
      raise Exception("FactorDB fail")
      factordb_available=True
    if verbose: print("factorDB available")
    return f.get_factor_list()
  except:
    factordb_available=False
    if verbose: print("factorDB not available")
#--------------------------------------------------------------------------


# returns all the factors (aka divisors) of a number
# e.g. all_factors(12) returns [1,2,3,4,6,12]

def all_factors( x ):
  if type(x)==list:
    allf = [1]
    for p in set(x):
      e=x.count(p)
      le = len(allf)
      pn = 1
      for i in range(e):
        pn *= p
        allf.extend([a*pn for a in allf[0:le] ])        
    return allf

  elif type(x)==int or type(x)==int:
    if CPP_available and x < 2**64:
      return allfactorCPP(x)
    else:
      return all_factors(factorise(x))

  
#-----------------------------------------------
def totient(n):
  if CPP_available and n < 2**64:
    return TotientCPP(n)

  t=n
  for p in set(factorise(n)):
    t -= t//p
  return t    

#-----------------------------------------------
def lcm(a,b):
  return (a*b)//gcd(a,b)

#-----------------------------------------------
# carmichael function, see https://en.wikipedia.org/wiki/Carmichael_function
def carmichael(n):
  c = 1
  f = factorise(n)
  for p in set(f):
    if p==2 and f.count(p)>=3:
      c = lcm(c, (p-1)*p**(f.count(p)-2))
    else:
      c = lcm(c, (p-1)*p**(f.count(p)-1))
      
  return c    


#-----------------------------------------------
def factorise_lookup(n):
  i = prime_or_factor[n]
  f=[]  
  while i :
    f.append(i)
    n//=i
    i = prime_or_factor[n]
  if n>1:
    f.append(n)
  return f

#-----------------------------------------------
def brent_gmp(n,maxBrentTime=1.0):
  if n%2==0:
    return 2

  start=perf_counter()
  # http://xn--2-umb.com/09/12/brent-pollard-rho-factorisation suggests that
  #   m=1000 is the optimum value
  m=1000

  # values to match c++ tests
  rs = gmpy2.random_state( int(binascii.hexlify(os.urandom(32)), 16) )
  #rs = gmpy2.random_state(int(datetime.now().microsecond) )

 
  y=gmpy2.mpz_random(rs, n-2) + 1
  y0=y

  # c=-5 from http://www.diva-portal.se/smash/get/diva2:1330398/FULLTEXT01.pdf
  c=gmpy2.mpz_random(rs, n-3) + 1 # random number in [1,n-3] (need to avoid n-2)

  N=mpz(n)
   
  g,q = mpz(1),mpz(1)
  r=1
  while g==1:             
    x = y
    
    for _ in range(r):
      y = (y*y+c)%N

    k = 0
    while (k<r and g==1):
      ys = y
      for _ in range(min(m,r-k)):
          y = (y*y+c)%N
          q = (q*(abs(x-y)))%N

      g = gmpy2.gcd(q,N)
      k +=  m
      
    r = r*2
    if other_methods_avaliable and r>2048:
      pcount = perf_counter()
      if pcount-start > maxBrentTime:
        if verbose: print ("pollard rho gmp timed out ", pcount-start, f",c=(c), y0={y0}")
        raise Exception("Brent gmp out of time, N={0}".format(N))  
 
  if g==N:
    while True:
      ys = (ys*ys+c)%N

      g = gmpy2.gcd(abs(x-ys),N)
      if g>1:
        break
  #print ( "brent_gmp found factor {} in {} microsec".format(int(g),int(1000000*(perf_counter()-start))  )  ,end='\n')
  if g==n:
    if verbose: print (f"pollard rho gmp cycle={n} for c={c}, y0={y0}")
  return int(g)  
#-----------------------------------------------

def brent(N,maxBrentTime=1.0):
  if gmpy2_available and N>2**64:
    return brent_gmp(N,maxBrentTime)
  
  if verbose: print("calling brent for N=",N)
  
  start=perf_counter()
  if N%2==0:
    return 2
  # http://xn--2-umb.com/09/12/brent-pollard-rho-factorisation suggests that
  #   m=1000 is the optimum value
  y = random.randint(1, N-1)
  y0=y
  c = random.randint(1, N-3)
  m=1000
  g,r,q = 1,1,1
  
  while g==1:             
    x = y
    
    for _ in range(r):
      y = (y*y+c)%N
      
    k = 0
    while (k<r and g==1):
      ys = y
      for _ in range(min(m,r-k)):
          y = (y*y+c)%N
          q = (q*(abs(x-y)))%N
      g = gcd(q,N)
      k +=  m
      
    r = r*2
    if other_methods_avaliable and r>2048:
      if perf_counter()-start > maxBrentTime:
        if verbose: print("Brent out of time, N={0}".format(N))  
        raise Exception("Brent out of time, N={0}".format(N),end='\n') 

  if g==N:
    while True:
      ys = (ys*ys+c)%N
      g = gcd(abs(x-ys),N)
      if g>1:
        break
  #if verbose: print ( "brent found factor {} in {} microsec".format(int(g),int(1000000*(perf_counter()-start))  ) ,end='\n')
  if g==N:
    if verbose: print (f"pollard rho cycle={N} for c={c}, y0={y0}")
  return g  
#-----------------------------------------------

def millerRabin(n, factor=[1]):
  global factorFromMillerRabin 

  if n == 2: return True
  if n<2: return False
  if n%2==0:
    factor[0]=2  
    return False

  two64 = pow(2,64)
  if n < 2047: P=(2,)
  elif n < 1373653: P=(2,3)
  elif n < 9080191: P=(31,73)
  elif n < 4759123141: P=(2,7,61)
  elif n < 2152302898747: P=(2,3,5,7,11)
  elif n < 3474749660383: P=(2,3,5,7,11,13)
  elif n < 341550071728321: P=(2,3,5,7,11,13,17)
  elif n < 3825123056546413051: P=(2,3,5,7,11,13,17,19,23)
  elif n < two64: P=(2, 325, 9375, 28178, 450775, 9780504, 1795265022) #http://miller-rabin.appspot.com/
  else: P=[2,3,5,7,11,13,17,19,23]+[random.randint(24,n-2) for i in range(40)]

  if n < 3825123056546413051 and n in P: return True

  # find s,d so that d is odd and (2^s)d = n-1
  d = (n-1)>>1
  s = 1
  while d&1 == 0:
    d >>= 1
    s += 1


  for w in P:
    a = pow(w, d, n)
    if a == 1: continue
    r = 0
    while r < s:
      if a == n-1: break
      prev_a = a
      a = (a*a)%n
      if a==1:
        # Previous value of a was a square root of unity not equal to -1.
        # We can use this to find a factor of n
        factor[0]=gcd(prev_a - 1, n)
        if verbose:
          if verbose: print("miller rabing found factor",factor[0],"for n=",n,"using witness",w)
        return False
      r += 1
    if r == s:
      # If we arrive here, then a^(n-1) != 1 mod n,
      # so n cannot be a prime, but we can't find a factor of n
      factor[0]=1
      return False

  return True

def is_prime(n):
  return millerRabin(n)

brent_failures=0
brent_calls=0
brent_avoid=0
#-----------------------------------------------

def factorise_large(n, maxBrentTime=1.0):
  global brent_failures,  brent_calls, brent_avoid, verbose

  if verbose:
    if verbose: print("calling factorise_large for n=",n)
  
  if n<=MAXLOOKUP:
    return factorise_lookup(n)

  # for recursive calls to factorise_large
  # when a factor has been found that makes one or both of the components < 2^64
  if CPP_available and n<2**64:
    return factoriseCPP(n)

  millerrabinfactor = [1]
  if millerRabin(n,millerrabinfactor):
    return [n]  

  if 1 < millerrabinfactor[0] < n : 
    p = millerrabinfactor[0] 
    brent_avoid+=1
    if verbose:
      if verbose: print("miller rabin found factor",p ,"for n=",n)
  else:
    #if verbose: print ("calling brent({})".format(n) ,end='\n')
    p=brent(n, maxBrentTime)
    brent_calls+=1

  q=n//p
  if p>q:
    p,q = q,p
  if p==1:
    if millerRabin(q):  
      return [q]
    else:
      brent_failures+=1
      if verbose: print ("re-trying brent for q={}".format(q))
      return factorise_large(q, maxBrentTime) # Try again (sometimes brent doesn't find factors)
  else:
    return factorise_large(p, maxBrentTime)+factorise_large(q, maxBrentTime)
#-----------------------------------------------

def factorise(n, maxBrentTime=1.0):
  if n<=MAXLOOKUP:
    return factorise_lookup(n)

  if CPP_available and n<2**64:
    return factoriseCPP(n)
  
  factors=[]
  for p in smallprimes:
    while n%p==0:
      factors.append(p)
      n//=p
    
  if n ==1:
    return factors
  
  if n <= MAXLOOKUP:
    return factors+factorise_lookup(n) 

  if CPP_available and n<2**64:
    return factors+factoriseCPP(n)

  try:
    return factors + factorise_large(n, maxBrentTime)
  except:
    if YAFU_available:
      if factordb_available:
        # try yafu for a few seconds, then factordb
        if verbose: print("trying yafu for n={}".format(n))
        try:
          return factors+yafu_factorise(n, timeout=5.0)
        except:
          if verbose: print("trying factordb for n={}".format(n))
          try:
            return factors+factordb_factorise(n)
          except:
            if verbose: print("trying yafu again for n={}".format(n))
            return factors+yafu_factorise(n, timeout=10000000.0)
           
      else:  
        # try yafu for a long time
        if verbose: print("trying yafu for n={}".format(n))
        return factors+yafu_factorise(n, timeout=10000000.0)

    if ECM_available and gmpy2_available and n < 10**50:
      if verbose: print("trying ecm for n={}".format(n))
      return factors+[int(x) for x in pyecm.factors( n,ra=True,ov=10,veb=False,pr=1)]

    if factordb_available:
      if verbose: print("trying factordb for n={}".format(n))
      return factors+factordb_factorise(n)

    raise Exception("unable to factorise {}".format(n) )
#-----------------------------------------------
def factorise_try(n, maxBrentTime=1.0):
  #attempt factorisation in a limitied time. Don't use ECM
  if n<=MAXLOOKUP:
    return factorise_lookup(n)

  if CPP_available and n<2**64:
    return factoriseCPP(n)
  
  factors=[]
  for p in smallprimes:
    while n%p==0:
      factors.append(p)
      n//=p
    
  if n ==1:
    return factors
  
  if n <= MAXLOOKUP:
    return factors+factorise_lookup(n) 

  if CPP_available and n<2**64:
    return factors+factoriseCPP(n)

  try:
    return factors+factorise_large(n, maxBrentTime)
  except:
    return factors
#-----------------------------------------------

def unfactorise( d ):
  return reduce(mul, d, 1)

def radical( d ):
  if type(d)==list:  
    return reduce(mul, set(d), 1)
  else:
    return reduce(mul, set(factorise(d)), 1)
#-----------------------------------------------
def multipart(n):
  # Multiplicative partition, a.k.a. factorisatio numerorum
  divisors = sorted(all_factors(n))
  divisors.remove(1)

  def multipart_inner(m, prev):
    if m==1:
      yield []
    elif m < prev:
      return
    else:  
      for d in divisors:
        if d >= prev:
          if m%d==0:
            for i in  multipart_inner(m//d, d):
              yield [d]+i
          if d==m:
            break

  for x in multipart_inner(n, 1):
    yield x
#-----------------------------------------------
 

