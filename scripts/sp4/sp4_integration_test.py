import unittest
import arrayfire as af
import magnum_af
import math

class sp4(unittest.TestCase):
  meshvar=magnum_af.Mesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
  m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)
  
  material=magnum_af.Material()
  material.ms    (8e5)
  material.A     (1.3e-11)
  material.alpha (1)
  
  m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
  m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
  m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
  pystate=magnum_af.State(meshvar,material,m)
  
  demag=magnum_af.DemagField(meshvar,material)
  exch=magnum_af.ExchangeField(meshvar,material)
  Llg=magnum_af.LLGIntegrator([pystate,demag,exch])

  def test_relaxation(self):
    intx=0
    inty=0
    intz=0
    while self.pystate.t() < 1e-9:
      self.Llg.llgstep(self.pystate)
      intx+=self.pystate.meanxyz(0)*self.Llg.print_stepsize()
      inty+=self.pystate.meanxyz(1)*self.Llg.print_stepsize()
      intz+=self.pystate.meanxyz(2)*self.Llg.print_stepsize()

    self.assertLess(math.fabs(intx - 9.81206172824e-10), 1e-15)
    self.assertLess(math.fabs(inty - 9.14350283169e-11), 1e-15)
    self.assertLess(math.fabs(intz + 5.74381785359e-13), 1e-15)

  def test_switch(self):
    self.Llg.set_state0_alpha(0.02)
    
    zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
    zeeswitch[0,0,0,0]=-24.6e-3/self.material.print_mu0()
    zeeswitch[0,0,0,1]=+4.3e-3/self.material.print_mu0()
    zeeswitch[0,0,0,2]=0.0
    zeeswitch = af.tile(zeeswitch,100,25,1)
    zee=magnum_af.ExternalField(zeeswitch)
    self.Llg.add_terms(zee)
    intx=0
    inty=0
    intz=0
     
    with open('../../Data/sp4.dat', 'w') as f:
      while self.pystate.t() < 2e-9:
        self.Llg.llgstep(self.pystate)
        intx+= self.pystate.meanxyz(0) * self.Llg.print_stepsize()
        inty+= self.pystate.meanxyz(1) * self.Llg.print_stepsize()
        intz+= self.pystate.meanxyz(2) * self.Llg.print_stepsize()
        f.write("%10.12f %10.12f %10.12f %10.12f\n" % (self.pystate.t(), self.pystate.meanxyz(0), self.pystate.meanxyz(1), self.pystate.meanxyz(2)))

    self.assertLess(math.fabs(intx + 6.41261165705e-10), 1e-15)
    self.assertLess(math.fabs(inty - 1.47353233738e-10), 1e-15)
    self.assertLess(math.fabs(intz + 1.78144535231e-11), 1e-15)

if __name__ == '__main__':
  unittest.main()
