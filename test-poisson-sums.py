import bayes
import time
import math

def test_Poisson_sums():
    print("Checking Poisson sums. Manually check that these are ~1:")
    print("Looping over k:")
    s=math.fsum([bayes.PoissonApprox2(x,10) for x in range(300)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(x,100) for x in range(3000)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(x,1000) for x in range(10000)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(x,2000) for x in range(20000)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(x,3000) for x in range(30000)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(x,200000) for x in range(800000)])
    print(s)

    print("Looping over lambda:")
    s=math.fsum([bayes.PoissonBetterApprox(1000,x) for x in range(3000)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(1000,x) for x in range(3000)])
    print(s)
    s=[bayes.PoissonApprox2(1000,x) for x in range(3000)]
    s.sort()
    s=math.fsum(s)
    print(s)
    s=[bayes.PoissonBetterApprox(100,x) for x in range(250)]
    s.sort()
    s=math.fsum(s)
    print(s)
    s=math.fsum([bayes.PoissonBetterApprox(100,x) for x in range(250)])
    print(s)
    s=math.fsum([bayes.PoissonBetterApprox(200,x) for x in range(350)])
    print(s)
    s=math.fsum([bayes.PoissonBetterApprox(10,x) for x in range(50)])
    print(s)
    s=[bayes.PoissonBetterApprox(10,x) for x in range(50)]
    s.sort()
    s=math.fsum(s)
    print(s)
    s=math.fsum([bayes.PoissonApprox2(10,x) for x in range(50)])
    print(s)
    s=math.fsum([bayes.PoissonBetterApprox(3,x) for x in range(50)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(3,x) for x in range(50)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(5000,x) for x in range(25000)])
    print(s)
    s=[bayes.PoissonApprox2(5000,x) for x in range(25000)]
    s.sort()
    s=math.fsum(s)
    print(s)
    s=math.fsum([bayes.PoissonApprox2(50000,x) for x in range(250000)])
    print(s)
    s=math.fsum([bayes.PoissonApprox2(100000,x) for x in range(300000)])
    print(s)
    s=math.fsum([bayes.Poisson(10,x) for x in range(200)])
    print(s)
    s=math.fsum([bayes.Poisson(4,x) for x in range(100)])
    print(s)
    s=math.fsum([bayes.Poisson(2,x) for x in range(100)])
    print(s)
    s=math.fsum([bayes.Poisson(2,x) for x in range(1000)])
    print(s)
    s=math.fsum([bayes.Poisson(1,x) for x in range(100)])
    print(s)
    s=math.fsum([bayes.Poisson(0,x) for x in range(100)])
    print(s)
    print("Done with Poisson sums.")

test_Poisson_sums()
