# Factoring in Python

A module to factor numbers in Python. The module uses several techniques to factor an integer, including 

   - Lookup tables for numbers below 10^6
   
   - Trial division of primes below 1200
   
   - Miller-Rabin primality test (deterministic up to the published limit)
   
   - Brent's variant of the Pollard Rho algorithm, by default running for up to 1 second
   
   - An interface to a DLL of a fast C++ implementation, for 64 bit numbers. The C++ implementation can be found at https://github.com/armchaircaver/Factors
   
   - Python interface to YAFU, running YAFU as a subprocess, to factorise very large numbers (YAFU can be obtained from https://sourceforge.net/projects/yafu/)
   
   - Using Martin G Kelly's implementation of the Elliptic Curve Method (ECM) if the Brent algorithm fails to find a factor within its time limit. In practice, ECM is only useful if GMP arithmetic is available (via module gmpy2). pyecm can be found at https://github.com/martingkelly/pyecm

The main functions of the module are:
- factorise(n, maxBrentTime=1.0) return a list of the prime factors of n. maxBrentTime specifies the number of seconds that the Pollard Rho algorithm should run for, before resorting to ECM or YAFU

- millerRabin(n) returns True if the Miller Rabin algorithm determines that n is probably prime, False if it identifies n as a composite. The algorithm is deterministic up to 2^64

- brent(N,maxBrentTime=1.0) implements the Pollard Rho algorithm with Brent's vaiation, and returns a prime factor of N, or raises and exception if it times out without finding a factor

- brent_gmp(n,maxBrentTime=1.0) GMP implemtation of the Pollard Rho algorithm. Requires module gmpy2, available at https://pypi.org/project/gmpy2/

- carmichael(n) returns the Carmichael number of n

- totient(n) returns the Euler totient number of n
- all_factors(n) returns all the factors, aka divisors of n, e.g. all_factors(12) returns [1,2,3,4,6,12]


A test suite is also provided in "factors test cases.py"



