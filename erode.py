# ##### BEGIN GPL LICENSE BLOCK #####
#
#  erode.py  -- a script to simulate erosion of height fields
#  (c) 2014 Michel J. Anders (varkenvarken)
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

from time import time
import unittest
import sys
import os
# import resource     # so much for platform independence, this only works on unix :-(
from random import random as rand, shuffle
import numpy as np

numexpr_available = False
try:
    import numexpr as ne
    numexpr_available = True
except ImportError:
    pass

def getmemsize():
  return 0.0
  #return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*resource.getpagesize()/(1024.0*1024.0)
  
def getptime():
  #r = resource.getrusage(resource.RUSAGE_SELF)
  #return r.ru_utime + r.ru_stime
  return time()
      
class Grid:

    def __init__(self, size=10, dtype=np.single):
        self.center = np.zeros([size,size], dtype)
        self.water = None
        self.sediment = None
        self.minx=None
        self.miny=None
        self.maxx=None
        self.maxy=None
        self.zscale=1
        self.maxrss=0.0
        self.sequence=[0,1,2,3]
 
    def init_water_and_sediment(self):
        if self.water is None:
            self.water = np.zeros(self.center.shape, dtype=np.single)
        if self.sediment is None:
            self.sediment = np.zeros(self.center.shape, dtype=np.single)

    def __str__(self):
        ''.join(self.__str_iter__(fmt=".3f"))
        
    def __str_iter__(self, fmt):
        for row in self.center[::]:
            values=[]
            for v in row:
                values.append(fmt%v)
            yield  ' '.join(values) + '\n'

    @staticmethod
    def fromFile(filename):
        if filename == '-' : filename = sys.stdin
        g=Grid()
        g.center=np.loadtxt(filename,np.single)
        return g

    def toFile(self, filename, fmt="%.3f"):
        if filename == '-' : 
            filename = sys.stdout.fileno()
        with open(filename,"w") as f:
            for line in self.__str_iter__(fmt):
                f.write(line)
    
    def raw(self,format="%.3f"):
        fstr=format+" "+ format+" "+ format+" "
        a=self.center / self.zscale
        minx=0.0 if self.minx is None else self.minx
        miny=0.0 if self.miny is None else self.miny
        maxx=1.0 if self.maxx is None else self.maxx
        maxy=1.0 if self.maxy is None else self.maxy
        dx=(maxx-minx)/a.shape[0]
        dy=(maxy-miny)/a.shape[1]
        for row in range(a.shape[0]-1):
            row0=miny+row*dy
            row1=row0+dy
            for col in range(a.shape[1]-1):
                col0=minx+col*dx
                col1=col0+dx
                yield (fstr%(row0 ,col0 ,a[row  ][col  ])+
                       fstr%(row0 ,col1 ,a[row  ][col+1])+
                       fstr%(row1 ,col0 ,a[row+1][col  ])+"\n")
                yield (fstr%(row0 ,col1 ,a[row  ][col+1])+
                       fstr%(row1 ,col0 ,a[row+1][col  ])+
                       fstr%(row1 ,col1 ,a[row+1][col+1])+"\n")

    def toRaw(self, filename, infomap=None):
        with open(filename if type(filename) == str else sys.stdout.fileno() , "w") as f:
            f.writelines(self.raw())
        if infomap:
            with open(os.path.splitext(filename)[0]+".inf" if type(filename) == str else sys.stdout.fileno() , "w") as f:
                f.writelines("\n".join("%-15s: %s"%t for t in sorted(infomap.items())))
            
    @staticmethod
    def fromRaw(filename):
        """initialize a grid from a Blender .raw file.
        currenly suports just rectangular grids of all triangles
        """
        g=Grid.fromFile(filename)
        #  we assume tris and an axis aligned grid
        g.center=np.reshape(g.center,(-1,3))
        #  keep unique vertices only by creating a set and sort first on x then on y coordinate
        #  using rather slow python sort but couldn;t wrap my head around np.lexsort
        verts = sorted(list({ tuple(t) for t in g.center[::] }))
        x=set(c[0] for c in verts)
        y=set(c[1] for c in verts)
        nx=len(x)
        ny=len(y)
        g.minx=min(x)
        g.maxx=max(x)
        g.miny=min(y)
        g.maxy=max(y)
        xscale=(g.maxx-g.minx)/(nx-1)
        yscale=(g.maxy-g.miny)/(ny-1)
        if abs(xscale/yscale) - 1.0 > 1e-3 : raise ValueError("Mesh spacing not square %d x %d  %.4f x %4.f"%(nx,ny,xscale,yscale))
        g.zscale=1.0/xscale 
        #  keep just the z-values and null any ofsset
        g.center=np.array([c[2] for c in verts],dtype=np.single).reshape(nx,ny)
        g.center=(g.center-np.amin(g.center))*g.zscale
        return g
        
    def peak(self, value=1):
        nx,ny = self.center.shape
        self.center[int(nx/2),int(ny/2)] += value

    def shelf(self, value=1):
        nx,ny = self.center.shape
        self.center[:nx/2] += value

    def mesa(self, value=1):
        nx,ny = self.center.shape
        self.center[nx/4:3*nx/4,ny/4:3*ny/4] += value

    def random(self, value=1):
        self.center += np.random.random_sample(self.center.shape)*value

    def neighborgrid(self):
        self.up=np.roll(self.center,-1,0)
        self.down=np.roll(self.center,1,0)
        self.left=np.roll(self.center,-1,1)
        self.right=np.roll(self.center,1,1)

    def zeroedge(self, quantity=None):
        c = self.center if quantity is None else quantity
        c[0,:]=0
        c[-1,:]=0
        c[:,0]=0
        c[:,-1]=0

    def diffuse(self, Kd, numexpr):
        self.zeroedge()
        c     = self.center[1:-1,1:-1]
        up    = self.center[ :-2,1:-1]
        down  = self.center[2:  ,1:-1]
        left  = self.center[1:-1, :-2]
        right = self.center[1:-1,2:  ]
        if(numexpr and numexpr_available):
            self.center[1:-1,1:-1] = ne.evaluate('c + Kd * (up + down + left + right - 4.0 * c)')
        else:
            self.center[1:-1,1:-1] = c + Kd * (up + down + left + right - 4.0 * c)
        self.maxrss = max(getmemsize(), self.maxrss)
        return self.center

    def avalanche(self, delta, numexpr):
        self.zeroedge()
        #print(self.center)
        
        c     = self.center[1:-1,1:-1]
        up    = self.center[ :-2,1:-1]
        down  = self.center[2:  ,1:-1]
        left  = self.center[1:-1, :-2]
        right = self.center[1:-1,2:  ]
        where = np.where
       
        if(numexpr and numexpr_available):
            self.center[1:-1,1:-1] = ne.evaluate('c + where((up   -c) > delta ,(up   -c -delta)/2, 0) \
                 + where((down -c) > delta ,(down -c -delta)/2, 0)  \
                 + where((left -c) > delta ,(left -c -delta)/2, 0)  \
                 + where((right-c) > delta ,(right-c -delta)/2, 0)  \
                 + where((up   -c) < -delta,(up   -c +delta)/2, 0)  \
                 + where((down -c) < -delta,(down -c +delta)/2, 0)  \
                 + where((left -c) < -delta,(left -c +delta)/2, 0)  \
                 + where((right-c) < -delta,(right-c +delta)/2, 0)')
        else:
            c = c + (
            # incoming
                   where((up   -c) > delta ,(up   -c -delta)/2, 0) 
                 + where((down -c) > delta ,(down -c -delta)/2, 0)
                 + where((left -c) > delta ,(left -c -delta)/2, 0)
                 + where((right-c) > delta ,(right-c -delta)/2, 0)
            # outgoing
                 + where((up   -c) < -delta,(up   -c +delta)/2, 0) 
                 + where((down -c) < -delta,(down -c +delta)/2, 0)
                 + where((left -c) < -delta,(left -c +delta)/2, 0)
                 + where((right-c) < -delta,(right-c +delta)/2, 0)
                 )
            self.center[1:-1,1:-1] = c
        
        #print(self.center)
        self.maxrss = max(getmemsize(), self.maxrss)
        return self.center

    def rain(self, amount, variance=0):
        self.water += (1.0 - np.random.random(self.water.shape) * variance) * amount

    def flow3(self, Kc, Ks, Kd, Ka, numexpr):

      zeros = np.zeros
      where = np.where
      min = np.minimum
      max = np.maximum
      abs = np.absolute
      arctan = np.arctan
      sin = np.sin
      
      center = (slice(   1,   -1,None),slice(   1,  -1,None))
      up     = (slice(None,   -2,None),slice(   1,  -1,None))
      down   = (slice(   2, None,None),slice(   1,  -1,None))
      left   = (slice(   1,   -1,None),slice(None,  -2,None))
      right  = (slice(   1,   -1,None),slice(   2,None,None))
        
      water = self.water
      rock = self.center
      sediment = self.sediment
      height = rock + water
      sc = where(water>0, sediment/water, 0)  ##!! this gives a runtime warning for division by zero
      sdw = zeros(water[center].shape)
      svdw = zeros(water[center].shape)
      sds = zeros(water[center].shape)
      angle = zeros(water[center].shape)
      #print(height[center])
      #print(water[center])
      for d in (up,down,left,right):
          if(numexpr and numexpr_available):
            hdd = height[d]
            hcc = height[center]
            dw = ne.evaluate('hdd-hcc')
            inflow = ne.evaluate('dw > 0')
            wdd = water[d]
            wcc = water[center]
            dw = ne.evaluate('where(inflow, where(wdd<dw, wdd, dw), where(-wcc>dw, -wcc, dw))/4.0') # nested where() represent min() and max()
            sdw  = ne.evaluate('sdw + dw')
            scd  = sc[d]
            scc  = sc[center]
            rockd= rock[d]
            rockc= rock[center]
            sds  = ne.evaluate('sds + dw * where(inflow, scd, scc)')
            svdw = ne.evaluate('svdw + abs(dw)')
            angle= ne.evaluate('angle + arctan(abs(rockd-rockc))')
          else:
            dw = (height[d]-height[center])
            inflow = dw > 0
            dw = where(inflow, min(water[d], dw), max(-water[center], dw))/4.0
            sdw  = sdw + dw
            sds  = sds + dw * where(inflow, sc[d], sc[center])
            svdw = svdw + abs(dw)
            angle= angle + np.arctan(abs(rock[d]-rock[center]))
          
      if(numexpr and numexpr_available):
          wcc = water[center]
          scc = sediment[center]
          rcc = rock[center]
          water[center] = ne.evaluate('wcc + sdw')
          sediment[center] = ne.evaluate('scc + sds')
          sc = ne.evaluate('where(wcc>0, scc/wcc, 2000*Kc)')
          fKc = ne.evaluate('Kc*sin(Ka*angle)*svdw')
          ds = ne.evaluate('where(sc > fKc, -Kd * scc, Ks * svdw)')
          rock[center] = ne.evaluate('rcc - ds')
          rock[center] = ne.evaluate('where(rcc<0,0,rcc)') # there isn't really a bottom to the rock but negative values look ugly
          sediment[center] = ne.evaluate('scc + ds')
      else:
          wcc = water[center]
          scc = sediment[center]
          rcc = rock[center]
          water[center] = wcc + sdw
          sediment[center] = scc + sds
          sc = where(wcc>0, scc/wcc, 2*Kc)
          fKc = Kc*np.sin(Ka*angle)*svdw
          ds = where(sc > fKc, -Kd * scc, Ks * svdw)
          rock[center] = rcc - ds
          rock[center] = where(rcc<0,0,rcc) # there isn't really a bottom to the rock but negative values look ugly
          sediment[center] = scc + ds
      
    def flow2(self, Kc, Ks, Kd, numexpr):
      
        min = np.minimum
        max = np.maximum
        where = np.where
        logical_not = np.logical_not
        logical_and = np.logical_and
        
        center = (slice(   1,   -1,None),slice(   1,  -1,None))
        up     = (slice(None,   -2,None),slice(   1,  -1,None))
        down   = (slice(   2, None,None),slice(   1,  -1,None))
        left   = (slice(   1,   -1,None),slice(None,  -2,None))
        right  = (slice(   1,   -1,None),slice(   2,None,None))
        
        
        Wv = self.water[center]
        Av = self.center[center]
        Sv = self.sediment[center]

        dirs = up,down,left,right
        shuffle(self.sequence)
        sdw = np.zeros(Wv.shape)
        for u in [ dirs[i] for i in self.sequence]:
            H  = self.center + self.water
            Hv = H[center]
            Hu = H[u]
            Wu = self.water[u]
            
            dw = max(min(Wv, Hv - Hu), -Wu)
            self.water[center] = Wv - dw
            sdw = sdw + dw
        
        #print(dw)
        #print('Water after\n',self.water)
        
        dwnegative = sdw<0
        # inflow
        self.center[center]   = where(dwnegative,  Av + Kd * Sv, Av)
        self.sediment[center] = where(dwnegative, (1 - Kd) * Sv, Sv)

        # outflow
        if 1:
            Cs = Kc * sdw
            dwpositive = logical_not(dwnegative)
            saturated = (Sv > Cs)
            sp = logical_and(dwpositive, saturated)
            usp= logical_and(dwpositive, logical_not(saturated))

            self.sediment[up] = where(sp, Su + Cs            , Su)
            self.center[center] = where(sp, Av + Kd * (Sv - Cs), Av)
            self.sediment[center] = where(sp, (1 - Kd)* (Sv - Cs), Sv)

            self.sediment[up] = where(usp, Su + Sv + Ks * (Cs - Sv), Su)
            self.center[center] = where(usp, Av + Ks * (Cs - Sv)     , Av)
            self.sediment[center] = where(usp, 0                       , Sv)

    def flow(self, Kc, Ks, Kd, numexpr):
        # some local variables to reduce typing and because
        # numexpr doesn't work with qualifiers (I.e. A.b )
        center= self.center
        water = self.water
        sediment=self.sediment
        height= center + water

        c     = height[1:-1,1:-1]
        up    = height[ :-2,1:-1]
        down  = height[2:  ,1:-1]
        left  = height[1:-1, :-2]
        right = height[1:-1,2:  ]
        where = np.where
        min   = np.minimum
        max   = np.maximum
        wc    = water[1:-1,1:-1]
        wup   = water[ :-2,1:-1]
        wdown = water[2:  ,1:-1]
        wleft = water[1:-1, :-2]
        wright= water[1:-1,2:  ]
        sc    = sediment[1:-1,1:-1]
        sup   = sediment[ :-2,1:-1]
        sdown = sediment[2:  ,1:-1]
        sleft = sediment[1:-1, :-2]
        sright= sediment[1:-1,2:  ]

        if(numexpr and numexpr_available):
            raise NotImplementedError("fluvial erosion not available in numexpr optimized form")
        else:
            sumda = np.zeros(c.shape, dtype=np.single)
            sumdw = np.zeros(c.shape, dtype=np.single)
            sumds = np.zeros(center.shape, dtype=np.single)
            sumdsc    = sumds[1:-1,1:-1]
            sumdsup   = sumds[ :-2,1:-1]
            sumdsdown = sumds[2:  ,1:-1]
            sumdsleft = sumds[1:-1, :-2]
            sumdsright= sumds[1:-1,2:  ]
            # waterflow depends on height difference but is limited
            # by the available amount of water.
            # division by 4 guarantees mass balance
            # positive value are inflow, negative ones outflow
            dw =  (where((up - c) > 0, min(up - c, wup), max(up - c, -wc) ) )
            sumdw = sumdw + dw/4.0
            Cs = -Kc * dw # for outflow (negative dw) this is a positive number!
            #                                 ++ inflow ++|++ outflow ++  ++ abrasion ++       ++ extra deposition ++
            sumdsc  = sumdsc  + where(dw > 0,    -Kd * sc, where(sc < Cs,                - sc ,   -Kd * (sc - Cs) - Cs) )/4.0
            sumda   = sumda   + where(dw > 0,     Kd * sc, where(sc < Cs,      -Ks * (Cs - sc),    Kd * (sc - Cs)     ) )/4.0
            # sediment transport to  neighboring cell
            sumdsup = sumdsup + where(dw > 0,           0, where(sc < Cs,  sc + Ks * (Cs - sc),                     Cs) )/4.0  # last part should be clamped to Cs as well I think
            
            dw =  (where((down - c) > 0, min(down - c, wdown), max(down - c, -wc) ) )
            sumdw = sumdw + dw/4.0
            Cs = -Kc * dw
            sumdsc    = sumdsc    + where(dw > 0,    -Kd * sc, where(sc < Cs,                - sc ,   -Kd * (sc - Cs) - Cs) )/4.0
            sumda     = sumda     + where(dw > 0,     Kd * sc, where(sc < Cs,      -Ks * (Cs - sc),    Kd * (sc - Cs)     ) )/4.0
            sumdsdown = sumdsdown + where(dw > 0,           0, where(sc < Cs,  sc + Ks * (Cs - sc),                     Cs) )/4.0
            
            dw =  (where((left - c) > 0, min(left - c, wleft), max(left - c, -wc) ) )
            sumdw = sumdw + dw/4.0
            Cs = -Kc * dw
            sumdsc    = sumdsc    + where(dw > 0,    -Kd * sc, where(sc < Cs,                - sc ,   -Kd * (sc - Cs) - Cs) )/4.0
            sumda     = sumda     + where(dw > 0,     Kd * sc, where(sc < Cs,      -Ks * (Cs - sc),    Kd * (sc - Cs)     ) )/4.0
            sumdsleft = sumdsleft + where(dw > 0,           0, where(sc < Cs,  sc + Ks * (Cs - sc),                     Cs) )/4.0
            
            dw =  (where((right - c) > 0, min(right - c, wright), max(right - c, -wc) ) )
            sumdw = sumdw + dw/4.0
            Cs = -Kc * dw
            sumdsc    = sumdsc    + where(dw > 0,    -Kd * sc, where(sc < Cs,                - sc ,   -Kd * (sc - Cs) - Cs) )/4.0
            sumda     = sumda     + where(dw > 0,     Kd * sc, where(sc < Cs,      -Ks * (Cs - sc),    Kd * (sc - Cs)     ) )/4.0
            sumdsright= sumdsright+ where(dw > 0,           0, where(sc < Cs,  sc + Ks * (Cs - sc),                     Cs) )/4.0
            
            #print(Kc, Ks, Kd)
            #print(water)
            #print(height)
            #print(sediment)
            #print("dw\n",sumdw)
            #print("ds\n",sumdsc)
            #print("du\n",sumdsup)
            #print("da\n",sumda)
            
            water   [1:-1,1:-1] += sumdw
            sediment[1:-1,1:-1] += sumdsc
            sediment[ :-2,1:-1] += sumdsup
            sediment[2:  ,1:-1] += sumdsdown
            sediment[1:-1, :-2] += sumdsleft
            sediment[1:-1,2:  ] += sumdsright
            center  [1:-1,1:-1] += sumda
            center  [1:-1,1:-1]  = where(center  [1:-1,1:-1]<0, 0, center  [1:-1,1:-1])
            
        self.maxrss = max(getmemsize(), self.maxrss)

    def fluvial_erosion(self, rain, Kc, Ks, Kd, Ka, numexpr):
        self.init_water_and_sediment()
        self.rain(rain)
        self.zeroedge(self.water)
        self.zeroedge(self.sediment)
        self.flow3(Kc, Ks, Kd, Ka, numexpr)
        
    def analyze(self):
      self.neighborgrid()
      # just looking at up and left to avoid needless doubel calculations
      slopes=np.concatenate((np.abs(self.left - self.center),np.abs(self.up - self.center)))
      return '\n'.join(["%-15s: %.3f"%t for t in [
        ('height average', np.average(self.center)),
        ('height median', np.median(self.center)),
        ('height max', np.max(self.center)),
        ('height min', np.min(self.center)),
        ('height std', np.std(self.center)),
        ('slope average', np.average(slopes)),
        ('slope median', np.median(slopes)),
        ('slope max', np.max(slopes)),
        ('slope min', np.min(slopes)),
        ('slope std', np.std(slopes))
        ]]
      )

class TestGrid(unittest.TestCase):

  def test_diffuse(self):
    g=Grid(5)
    g.peak(1)
    self.assertEqual(g.center[2,2],1.0)
    g.diffuse(0.1, numexpr=False)
    for n in [(2,1),(2,3),(1,2),(3,2)]:
      self.assertAlmostEqual(g.center[n],0.1)
    self.assertAlmostEqual(g.center[2,2],0.6)

  def test_diffuse_numexpr(self):
    g=Grid(5)
    g.peak(1)
    g.diffuse(0.1, numexpr=False)
    h=Grid(5)
    h.peak(1)
    h.diffuse(0.1, numexpr=True)
    self.assertEqual(list(g.center.flat),list(h.center.flat))

  def test_avalanche_numexpr(self):
    g=Grid(5)
    g.peak(1)
    g.avalanche(0.1, numexpr=False)
    h=Grid(5)
    h.peak(1)
    h.avalanche(0.1, numexpr=True)
    print(g)
    print(h)
    np.testing.assert_almost_equal(g.center,h.center)

if __name__ == "__main__":

  import argparse

  parser = argparse.ArgumentParser(description='Erode a terrain while assuming zero boundary conditions.')
  parser.add_argument('-I', dest='iterations', type=int, default=1, help='the number of iterations')
  parser.add_argument('-Kd', dest='Kd', type=float, default=0.01, help='Diffusion constant')
  parser.add_argument('-Kh', dest='Kh', type=float, default=6, help='Maximum stable cliff height')
  parser.add_argument('-Kp', dest='Kp', type=float, default=0.1, help='Avalanche probability for unstable cliffs')
  parser.add_argument('-Kr', dest='Kr', type=float, default=0.1, help='Average amount of rain per iteration')
  parser.add_argument('-Kdep', dest='Kdep', type=float, default=0.1, help='Sediment deposition constant')
  parser.add_argument('-Ks', dest='Ks', type=float, default=0.1, help='Soil softness constant')
  parser.add_argument('-Kc', dest='Kc', type=float, default=1.0, help='Sediment capacity')
  parser.add_argument('-Ka', dest='Ka', type=float, default=2.0, help='Slope dependency of erosion')
  parser.add_argument('-ri', action='store_true', dest='rawin', default=False, help='use Blender raw format for input')
  parser.add_argument('-ro', action='store_true', dest='rawout', default=False, help='use Blender raw format for output')
  parser.add_argument('-i',  action='store_true', dest='useinputfile', default=False, help='use an inputfile (instead of just a synthesized grid)')
  parser.add_argument('-t',  action='store_true', dest='timingonly', default=False, help='do not write anything to an output file')
  parser.add_argument('-infile', type=str, default="-", help='input filename')
  parser.add_argument('-outfile', type=str, default="-", help='output filename')
  parser.add_argument('-Gn', dest='gridsize', type=int, default=20, help='Gridsize (always square)')
  parser.add_argument('-Gp', dest='gridpeak', type=float, default=0, help='Add peak with given height')
  parser.add_argument('-Gs', dest='gridshelf', type=float, default=0, help='Add shelve with given height')
  parser.add_argument('-Gm', dest='gridmesa', type=float, default=0, help='Add mesa with given height')
  parser.add_argument('-Gr', dest='gridrandom', type=float, default=0, help='Add random values between 0 and given value')
  parser.add_argument('-m', dest='threads', type=int, default=1, help='number of threads to use')
  parser.add_argument('-u', action='store_true', dest='unittest', default=False, help='perfom unittests')
  parser.add_argument('-a', action='store_true', dest='analyze', default=False, help='show some statistics of input and output meshes')
  parser.add_argument('-d', action='store_true', dest='dump', default=False, help='show sediment and water meshes at end of run')
  parser.add_argument('-n', action='store_true', dest='usenumexpr', default=False, help='use numexpr optimizations')

  args = parser.parse_args()
  print("\nInput arguments:")
  print("\n".join("%-15s: %s"%t for t in sorted(vars(args).items())), file=sys.stderr)

  if args.unittest:
    unittest.main(argv=[sys.argv[0]])
    sys.exit(0)
  
  if args.useinputfile:
    if args.rawin:
      grid = Grid.fromRaw(args.infile)
    else:
      grid = Grid.fromFile(args.infile)
  else:
    grid = Grid(args.gridsize)

  if args.gridpeak   > 0 : grid.peak(args.gridpeak)
  if args.gridmesa   > 0 : grid.mesa(args.gridmesa)
  if args.gridshelf  > 0 : grid.shelf(args.gridshelf)
  if args.gridrandom > 0 : grid.random(args.gridrandom)

  if args.analyze:
    print('\nstatistics of the input grid:\n\n', grid.analyze(), file=sys.stderr, sep='' )
  t = getptime()
  for g in range(args.iterations):
    if args.Kd > 0:
      grid.diffuse(args.Kd, args.usenumexpr)
    if args.Kh > 0 and args.Kp > rand():
      grid.avalanche(args.Kh, args.usenumexpr)
    if args.Kr > 0:
      grid.fluvial_erosion(args.Kr, args.Kc, args.Ks, args.Kdep, args.Ka, args.usenumexpr)
  t = getptime() - t
  print("\nElapsed time: %.1f seconds, max memory %.1f Mb.\n"%(t,grid.maxrss), file=sys.stderr)
  if args.analyze:
    print('\nstatistics of the output grid:\n\n', grid.analyze(), file=sys.stderr, sep='')

  if not args.timingonly:
    if args.rawout:
      grid.toRaw(args.outfile, vars(args))
    else:
      grid.toFile(args.outfile)

  if args.dump:
    print("sediment\n", np.array_str(grid.sediment,precision=3), file=sys.stderr)
    print("water\n", np.array_str(grid.water,precision=3), file=sys.stderr)
    print("sediment concentration\n", np.array_str(grid.sediment/grid.water,precision=3), file=sys.stderr)
