from factors import factorise, all_factors, prime_or_factor, brent_gmp, totient, carmichael, unfactorise, factoriseCPP, millerRabin as isprime, setCPP, brent, factorise_large, ECM_available, YAFU_available
from time import perf_counter
import random



def pathological():
  # a few pathological cases
  
  tests = (
    1023371,    2**155-1 ,    2**163-1,    2**198-1 ,    2**157-1 ,
    10000000000000000000000000000006
           )

  brenttests = ( 1351369, 1023371 )
  for n in brenttests:
    p=factorise_large(n)
    print ("brent test n=",n, ", p=",p)

  totaltime=0
  for n in tests:
    
    start=perf_counter()
    f=factorise(n,3.0)
    end=perf_counter()

    print ("\n",n)
    print (f)
    print (int( (end-start)*1000000 ), "microsec")
    totaltime += int( (end-start)*1000000 )
                      
  print ("totaltime", totaltime)
#-----------------------------------------------------------------------------

def factorise_trial_division(n):
  if n==1:
    return []

  factors=[]
  for p in (2,3):
    while n%p==0:
      factors.append(p)
      n//=p
  q=5
  while q*q <= n:
    for p in (q,q+2):
      while n%p==0:
        factors.append(p)
        n//=p
    q+=6
  if n>1:
    factors.append(n)
  return factors

#----------------------------------------------------------------------------------
def test_cases():

  global brent_failures,  brent_calls, brent_avoid

  print ("brent_gmp tests...")
  n = brent_gmp(21399691250168634880751)
  print (n)
  
  print("\nTotient tests...")  
  A000010 = (0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8, 12, 10,
            22, 8, 20, 12, 18, 12, 28, 8, 30, 16, 20, 16, 24, 12, 36, 18, 24, 16, 40,
            12, 42, 20, 24, 22, 46, 16, 42, 20, 32, 24, 52, 18, 40, 24, 36, 28, 58,
            16, 60, 30, 36, 32, 48, 20, 66, 32, 44)
  for n in range(1,10000):
    t=totient(n)
    if (n<len(A000010) and t !=  A000010[n]):
      print("***MISMATCH***, n=",n,", totient=",totient(n),", A000010=",A000010[n],", factors =",factorise(n))
      raise Exception("Mismatch")
  print("Totient tests completed")  


  print("\nCarmichael tests...")  
  A002322 = (0, 1, 1, 2, 2, 4, 2, 6, 2, 6, 4, 10, 2, 12, 6, 4, 4, 16, 6, 18, 4, 6,
             10, 22, 2, 20, 12, 18, 6, 28, 4, 30, 8, 10, 16, 12, 6, 36, 18, 12, 4,
             40, 6, 42, 10, 12, 22, 46, 4, 42, 20, 16, 12, 52, 18, 20, 6, 18, 28,
             58, 4, 60, 30, 6, 16, 12, 10, 66, 16, 22, 12, 70, 6, 72, 36, 20, 18,
             30, 12, 78, 4, 54)
  for n in range(1,10000):
    c=carmichael(n)
    if (n<len(A002322) and c !=  A002322[n]):
      print("***MISMATCH***, n=",n,", carmichael=",carmichael(n),", A002322=",A002322[n],", factors =",factorise(n))
      raise Exception("Mismatch")
  print("Carmichael tests completed")  

  
  print("\n64 bit numbers:")
  for N in ( 18446743721522234449, 10000000000000000049, 100000000000000000092, 100000000000000000092//12):
    print(N, factorise(N))

  
  
  print("\nmersenne test ")  
  A000043 = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127] # mersenne primes
  
  for N in range(2,200):
    x=2**N-1
    if N<61:
      start=perf_counter()  
      f=sorted(factorise_trial_division(x))
      end=perf_counter()
      if (x != unfactorise(f)):
        print("***Trial division MISMATCH***, f=",f,", n=",x)
        raise Exception("Trial division Mismatch ")
    else:
      start=0.0
      end=0.0
      
    if N in []: # [137,193]: # difficult semiprimes
      print("2^{0}-1 \t skipped".format(N))
    else:  
      start2=perf_counter()  
      f2=sorted(factorise(x))
      end2=perf_counter()

      if N<61:
        if f!=f2:
          print("***MISMATCH***, f=",f,", f2=",f2)
          raise Exception("Mismatch between factorise and trial division")
      assert( (N in A000043) == (len(f2)==1) )
      print("2^{0}-1".format(N),"\t","\t",f2,"\t",end='')
      #print(int(round((end-start)*1000000,0)), end=' ')
      print("\t",int(round((end2-start2)*1000000,0)),"micro s")

  
  print("mersenne test completed")  

  print("\npowers of primes")
  for p in [2,3,5]:
      x=1
      for n in range(100):
        x*=p
        assert(x == unfactorise(factorise(x)))
  print("powers of primes tests completed")


  print("\nfirst n digit prime ")
  A003617 = set((  2,11,101,1009,10007,100003,1000003,10000019,100000007,1000000007,
                 10000000019,100000000003,1000000000039,10000000000037,
                 100000000000031,1000000000000037,10000000000000061,100000000000000003,
                 1000000000000000003))
  for e in range(1,19):
      n=10**e+1
      while True:
        f=factorise(n)
        if len(f)==1:
          print(n,"is prime")
          assert(n in A003617)
          break
        n+=2
  print("first n digit prime test completed")

  print("\nCarmichael numbers with 9 factors")
  A112432 = (9746347772161, 11537919313921, 11985185775745, 14292786468961,
             23239986511105, 24465723528961, 26491881502801, 27607174936705,
             30614445878401, 30912473358481, 34830684315505, 51620128928641)
  for n in A112432:
    f=factorise(n)
    print(n, f)
    assert(len(f)==9)
  print("Carmichael numbers with 9 factors test completed")

  print("\n3...31 is prime?")
  p=31
  i=1
  stillPrime=True
  while stillPrime:
    f=factorise(p)
    print(p, f)
    stillPrime = (len(f)==1)
    i+=1
    p+=3*10**i
  print("3...31 tests completed")

  print("\nFermat numbers, 2^(2^n)+1")
  p=3
  while p<2**65:
    f=factorise(p)
    print(p, f)
    p=(p-1)**2+1
  print("\nFermat numbers tests completed")
  
  A066386 = set((0, 3906, 4620, 5166, 5376, 5460, 8190))
  print("\nn^6+1091 composite for n=1..3905")
  n=0
  while n<=3906:
    p=n**6+1091
    f=factorise(p)
    if len(f)==1:
      print(n, p, f)  
    assert( (len(f)==1) == (n in A066386) )
    n+=1
  print("n^6+1091 tests completed")

  print("\nCarmichael numbers A202562 ")
  A202562 =[561,84350561,851703301,2436691321,34138047673, 60246018673,63280622521,
            83946864769,110296864801, 114919915021,155999871721,225593397919,
            342267565249,534919693681,660950414671, 733547013841,1079942171239,1301203515361,
            1333189866793]
  for n in A202562:
    start=perf_counter()  
    f=factorise(n)
    end=perf_counter()
    print(n, f, int( (end-start)*1000000 )," micro s")
    assert(len(f)==3)
  print("Carmichael numbers A202562 test completed")

  assert( radical(2*2*3*5*5*7)==2*3*5*7 )
  
  print(brent_failures,"brent failures out of", brent_calls, "(",float(brent_failures)*100.0/float(brent_calls),"%)")
  print(brent_avoid, "brent avoid")
#----------------------------------------------------------------------

def maxsmallprime_test():
  ##############################################
  #  need to disable CPP for this test to work
  ##############################################
  global CPP_available, YAFU_available, ECM_available, factordb_available
  global smallprimes, MAXSMALLPRIME
  
  CPP_available=False
  YAFU_available=False
  ECM_available=False
  factordb_available=False

  for MAXSMALLPRIME in (100,300,500,1000,1500,3000,5000,10000,50000):
    print("microsec/number for MAXSMALLPRIME=",MAXSMALLPRIME,":", end=' ')
    smallprimes = [p for p in range(2,MAXSMALLPRIME) if prime_or_factor[p]==0]

    brent_failures = brent_calls = brent_avoid = 0
    # small numbers are dealt with via factirose_lookup
    for e in range(6,19,2):
      N=10**e  
      start=perf_counter()
      tim=0.0
      n=N
      while tim<=2.0 :
        facs=factorise(n,10)
        #print n, facs
        #if not(n == unfactorise(facs)):
        #  print n,"is not the product of its factors",facs
        #  exit(1)
        n+=1
        tim=perf_counter()-start
      microsec_per_number =round((tim*1000000.0)/(n-N),1)
      print(microsec_per_number, end=' ')
    print("")
  print("re-enabling CPP")
  try:
    from factorsdll import factoriseCPP, allfactorCPP, isprimeCPP, TotientCPP
    CPP_available=True
  except:
    pass
#------------------------------------------------------------------------------
    
def performance_test():

  global brent_failures,  brent_calls, brent_avoid
  global smallprimes, MAXSMALLPRIME

  print("\nfactor 10^k+i, highlighting harder cases")
  for k in range(27,33):
    print("k={}".format(k))
    n=10**k
    while n< 10**k+100:
      start=perf_counter()
      f=factorise(n,0.1)
      end=perf_counter()
      if end-start > 0.5:
        print(end-start, "sec to factorise",n,f)
      n+=1
    
  

  if CPP_available:
    print("\nmiller rabin timings for radom samples in range [2^n ..2^n+1] (times in microsec)")
    print("n\tPython\tC++")
    samples = 1000
    y=0
    for i in range(7,64):
      start=perf_counter()
      for j in range(samples):
        x = random.randint( pow(2,i),pow(2,i+1)-1 )
        y += millerRabin(x)
      python_time = (perf_counter()-start)*1000000.0/samples

      y=0
      start=perf_counter()
      for j in range(samples):
        x = random.randint( pow(2,i),pow(2,i+1)-1 )
        y += isprimeCPP(x)
      cpp_time = (perf_counter()-start)*1000000.0/samples

      print(i,"\t", round(python_time,1), "\t", round(cpp_time,1) )


  print("miller rabing timings complete")
    #print brent_failures,"brent failures out of", brent_calls, "(",float(brent_failures)*100.0/float(brent_calls),"%)"
    #print brent_avoid, "brent avoid"
#-----------------------------------------------------------------------------------------------
  

def CPP_comparison_test():
  
  print("\nperformance timings for random samples in range [2^n ..2^n+1] (times in microsec)")
  print("n\tPython\tC++\tratio")
  samples = 1000

  for i in range(7,65,2):
    tests=[]
    for j in range(samples):
      tests.append( random.randint( pow(2,i),pow(2,i+1)-1 ) )

    y_python = []
    setCPP(False)

    start=perf_counter()
    for x in tests:
      y_python.append(factorise(x))
    python_time = (perf_counter()-start)*1000000.0/samples

    y_cpp = []
    setCPP(True)
    start=perf_counter()
    for x in tests:
      y_cpp.append(factoriseCPP(x))
    cpp_time = (perf_counter()-start)*1000000.0/samples

    print(i,"\t", round(python_time,1), "\t",
          round(cpp_time,1),"\t",round(python_time/cpp_time,1) )

    for i in range(len(y_python)):
      if not ( sorted(y_python[i]) == sorted(y_cpp[i])):
        print (y_python[i], y_cpp[i])
        raise Exception("mismatch in python vs cpp factorisation")
      

  assert(len(y_python)==len(y_cpp))
  for i in range(len(y_cpp)):
    assert ( sorted(y_python[i]) == sorted(y_cpp[i]) )
    
  print("tests completed")
#----------------------------------------------------------------------------------------------

def mersenne_test():
  
  print("\nperformance timings for mersenne numbers 2^i-1")
  samples = 1000

  for i in range(63,300,2):
    n=2**i-1
    start=perf_counter()
    f=factorise(n)
    python_time = (perf_counter()-start)*1000000.0/samples

    print (f"2^{i}-1, {round(python_time,2)}ms :", f)
    
  print("tests completed")
  

if __name__=="__main__":

    print("ECM available", ECM_available)
    print("Yafu available", YAFU_available)

    verbose=False;

    #mersenne_test()
    #CPP_comparison_test()
    #maxsmallprime_test()     
    test_cases()
    performance_test()

  
