import unittest
import bayes
import math
import time

class TestMain(unittest.TestCase):

	def test_Poisson(self):
		print("Checking Poisson probabilities")
		print(bayes.Poisson(100, 1.4))
		print(bayes.Poisson(100, 1.4))
		print(bayes.Poisson(5, 1.4))
		self.assertEqual(bayes.Poisson(5, 1.4),0.011052147127910878)
		print("Done with Poisson probabilities")

	def test_Poisson_sums(self):
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
		print("Done with Poisson sums.")

	def test_Poisson_sum_times(self):
		print("Testing Poisson sum times")
		k_max=400000
		lmbda=100000
		print("k_max=%s, lmbda=%s" % (k_max,lmbda))
		gt1=time.clock()
		s=math.fsum([bayes.PoissonApprox2(x,lmbda) for x in range(k_max)])
		gt2=time.clock()
		print("Time: %s seconds" % (gt2-gt1,))
		k_max=400000
		lmbda=200000
		print("k_max=%s, lmbda=%s" % (k_max,lmbda))
		gt1=time.clock()
		s=math.fsum([bayes.PoissonApprox2(x,lmbda) for x in range(k_max)])
		gt2=time.clock()
		print("Time: %s seconds" % (gt2-gt1,))
		k_max=600000
		lmbda=200000
		print("k_max=%s, lmbda=%s" % (k_max,lmbda))
		gt1=time.clock()
		s=math.fsum([bayes.PoissonApprox2(x,lmbda) for x in range(k_max)])
		gt2=time.clock()
		print("Time: %s seconds" % (gt2-gt1,))
		k_max=800000
		lmbda=200000
		print("k_max=%s, lmbda=%s" % (k_max,lmbda))
		gt1=time.clock()
		s=math.fsum([bayes.PoissonApprox2(x,lmbda) for x in range(k_max)])
		gt2=time.clock()
		print("Time: %s seconds" % (gt2-gt1,))

if __name__ == '__main__':
    unittest.main()
