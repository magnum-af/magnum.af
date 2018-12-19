import unittest
import arrayfire as af
import magnum_af
import math

class sp4(unittest.TestCase):
  meshvar=magnum_af.pyMesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
  m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)
  
  param=magnum_af.pyParam()
  param.set_ms    (8e5)
  param.set_A     (1.3e-11)
  param.set_alpha (1)
  
  m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
  m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
  m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
  pystate=magnum_af.pyState(meshvar,param,m)
  
  demag=magnum_af.pyDemagSolver(meshvar,param)
  exch=magnum_af.pyExchSolver(meshvar,param)
  Llg=magnum_af.pyLLG(demag,exch)

  def test_relaxation(self):
    intx=0
    inty=0
    intz=0
    while self.pystate.t() < 1e-9:
      t1=self.pystate.t()
      self.Llg.llgstep(self.pystate)
      t2=self.pystate.t()
      stepsize=t2-t1
      intx+=self.pystate.meanxyz(0)*stepsize
      inty+=self.pystate.meanxyz(1)*stepsize
      intz+=self.pystate.meanxyz(2)*stepsize

    self.assertLess(math.fabs(intx - 9.81206172824e-10), 1e-15)
    self.assertLess(math.fabs(inty - 9.14350283169e-11), 1e-15)
    self.assertLess(math.fabs(intz + 5.74381785359e-13), 1e-15)

  def test_switch(self):
    self.pystate.set_alpha(0.02)
    
    zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
    zeeswitch[0,0,0,0]=-24.6e-3/self.param.get_mu0()
    zeeswitch[0,0,0,1]= +4.3e-3/self.param.get_mu0()
    zeeswitch[0,0,0,2]=0.0
    zeeswitch = af.tile(zeeswitch,100,25,1)
    zee=magnum_af.pyZee(zeeswitch)
    self.Llg.add_terms(zee)
    intx=0
    inty=0
    intz=0
    while self.pystate.t() < 2e-9:
      t1=self.pystate.t()
      self.Llg.llgstep(self.pystate)
      t2=self.pystate.t()
      stepsize=t2-t1
      intx+= self.pystate.meanxyz(0) * stepsize
      inty+= self.pystate.meanxyz(1) * stepsize
      intz+= self.pystate.meanxyz(2) * stepsize

    self.assertLess(math.fabs(intx + 6.41261165705e-10), 1e-15)
    self.assertLess(math.fabs(inty - 1.47353233738e-10), 1e-15)
    self.assertLess(math.fabs(intz + 1.78144535231e-11), 1e-15)

if __name__ == '__main__':
  unittest.main()
