import unittest
import arrayfire as af
import magnumaf
import math

class sp4(unittest.TestCase):
  mesh=magnumaf.Mesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
  m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)
  
  material=magnumaf.Material()
  A =1.3e-11
  
  m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
  m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
  m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
  state=magnumaf.State(mesh, Ms = 8e5, m = m)
  
  demag=magnumaf.DemagField(mesh)
  exch=magnumaf.ExchangeField(A)
  Llg=magnumaf.LLGIntegrator(alpha = 1, terms = [demag,exch])

  def test_relaxation(self):
    intx=0
    inty=0
    intz=0
    while self.state.t < 1e-9:
      t1=self.state.t
      self.Llg.step(self.state)
      t2=self.state.t
      stepsize=t2-t1
      intx+=self.state.m_mean(0)*stepsize
      inty+=self.state.m_mean(1)*stepsize
      intz+=self.state.m_mean(2)*stepsize

    self.assertLess(math.fabs(intx - 9.81206172824e-10), 1e-15)
    self.assertLess(math.fabs(inty - 9.14350283169e-11), 1e-15)
    self.assertLess(math.fabs(intz + 5.74381785359e-13), 1e-15)

  def test_switch(self):
    self.Llg.alpha = 0.02
    
    zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
    zeeswitch[0,0,0,0]=-24.6e-3/magnumaf.Constants.mu0
    zeeswitch[0,0,0,1]= +4.3e-3/magnumaf.Constants.mu0
    zeeswitch[0,0,0,2]=0.0
    zeeswitch = af.tile(zeeswitch,100,25,1)
    zee=magnumaf.ExternalField(zeeswitch)
    self.Llg.add_terms(zee)
    intx=0
    inty=0
    intz=0
    while self.state.t < 2e-9:
      t1=self.state.t
      self.Llg.step(self.state)
      t2=self.state.t
      stepsize=t2-t1
      intx+= self.state.m_mean(0) * stepsize
      inty+= self.state.m_mean(1) * stepsize
      intz+= self.state.m_mean(2) * stepsize

    self.assertLess(math.fabs(intx + 6.41261165705e-10), 1e-15)
    self.assertLess(math.fabs(inty - 1.47353233738e-10), 1e-15)
    self.assertLess(math.fabs(intz + 1.78144535231e-11), 1e-15)

if __name__ == '__main__':
  unittest.main()
