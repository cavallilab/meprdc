from numpy import *
from sys import argv,exit
from pdb import set_trace
from copy import deepcopy as dc
from scipy.optimize import fmin
def exp_pdf(x,a,p,sign):
     return exp(sign*dot(a,x**p))

def lexp_pdf(x,a,p,sign):
     return sign*dot(a,x**p)


class ScaleEstimator:
     """  stochastic optimization of scale vectors """
     def __init__(self, data, samples, errors=None, verbose=False, scale_invariant_optimization=False):
          self.iter=0
          self.sio=scale_invariant_optimization
          self.powavg=1.
          self.sign=-1.
          self.samples=samples
          self.data=data
          if errors:
            self.errors=errors
          else:
            self.errors=[ones(shape(d)) for d in data]
          self.verbose=verbose
          self.means=[mean(s,axis=1) for s in self.samples] 
          self.means_old=dc(self.means)
          self.alla=[s*0+1e-8 for s in self.data]
          self.old_alla=dc(self.alla)
          self.calc_mean()
          self.update_chi()
          self.accept()
          self.chis=[]
          self.small=1e-2
          self.allas=[]
          if self.verbose:
             print "chisq", self.chisq, self.chisq_old

     def calc_mean(self):
         lprobs=[lexp_pdf(sample,m,self.powavg,self.sign) for sample,m in zip(self.samples, self.alla)]
         probs=[exp(l-l.max()) for l in lprobs]
         self.probs=probs
         try:
           self.means=[average(s**(self.powavg), weights=p[0], axis=1)**(1./self.powavg) for a,s,p in zip(self.alla,self.samples,probs)]
         except:
           if self.verbose:
              print "Warning at iteration "+str(self.iter)+", weights possibly at extreme values, consider stopping estimation"
              print self.probs
     def sample_as(self):
          if self.sio:
            s=25.585
            self.alla=[ self.alla[i]+abs(random.normal(0,self.small, 1))*(self.means[i]*s-self.data[i]) for i in range(len(self.data))]
          else:
            self.alla=[ self.alla[i]+abs(random.normal(0,self.small, 1))*(self.means[i]-self.data[i]) for i in range(len(self.data))]
     def update_chi(self):
          if self.sio:
            self.chisq=sum([ corrcoef(self.means[i],self.data[i])[0,1] for i in range(len(self.data))])
          else:
            self.chisq=sum([ dot((self.means[i]-self.data[i]),dot(diag(1./self.errors[i][0]), (self.means[i]-self.data[i]).T ))/0.1 for i in range(len(self.data))])

     def mh_step(self):
          self.accept()
          self.sample_as()
          self.calc_mean()
          self.update_chi()
          if self.verbose:
            print "=================="
            print "proposed chi", self.chisq
          if self.sio:
            if self.chisq>self.chisq_old:
              self.accept()
            else:
              self.reject()                                                                               
          else:  
            aprob=min(1, exp(-self.chisq+self.chisq_old))
            if random.rand()<aprob:
                 self.accept()
                 if self.verbose:
                    print "accept"
                    print "aprob", aprob
            else:
                 self.reject()
                 if self.verbose:
                    print "reject"
                    print "aprob", aprob
          if self.iter%10==0:
              self.chis.append(self.chisq)
              self.allas.append(dc(self.alla))
          self.iter=self.iter+1

     def accept(self):
          self.old_alla=dc(self.alla)
          self.chisq_old=dc(self.chisq)
          self.means_old=dc(self.means)

     def reject(self):
          self.alla=dc(self.old_alla)
          self.chisq=dc(self.chisq_old)
          self.means=dc(self.means_old)

     def get_scales(self):
         return self.alla

def print_usage():
   print "Scale estimator."
   print argv[0]+" experimental-data samples-table <uncertainty>"

if __name__=="__main__":
  if len(argv)<2: 
    print_usage()
    exit()
  print "loading data and samples"
  data=loadtxt(argv[1])
  samples=loadtxt(argv[2])
  try:
      scale_estimator_instance = ScaleEstimator([data], [samples], [loadtxt(argv[3])], verbose=False)
  except:
      scale_estimator_instance = ScaleEstimator([data], [samples], verbose=False)

  for i in range(250):
      scale_estimator_instance.mh_step()

  [savetxt(argv[1]+".scales"+str(i), scale) for i,scale in enumerate(scale_estimator_instance.get_scales())]

