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
from random import random as rand
import numpy as np

numexpr_available = False
try:
    import numexpr as ne
    numexpr_available = True
except ImportError:
    pass
    
class Grid:

    def __init__(self, size=10, dtype=np.single):
        self.center = np.zeros([size,size], dtype)
        self.minx=None
        self.miny=None
        self.maxx=None
        self.maxy=None
 
    def __str__(self):
        values=[]
        for row in self.center[::]:
            for v in row:
                values.append('%.4f'%v+' ')
            values.append('\n')
        return ''.join(values)

    @staticmethod
    def fromFile(filename):
        if filename == '-' : filename = sys.stdin
        g=Grid()
        g.center=np.loadtxt(filename,np.single)
        return g

    def toFile(self, filename, fmt="%.3f"):
        if filename == '-' : 
            filename = sys.stdout
        print(filename)
        np.savetxt(filename,self.center)
    
    @staticmethod        
    def raw(self,format="%.3f"):
        fstr=format+" "+ format+" "+ format+" "
        a=self.center
        minx=0.0 if self.minx is None else self.minx
        miny=0.0 if self.miny is None else self.miny
        maxx=0.0 if self.maxx is None else self.maxx
        maxy=0.0 if self.maxy is None else self.maxy
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

    def toRaw(self, filename):
        with open(filename if type(filename) == str else sys.stdout.fileno() , "w") as f:
            f.writelines(self.raw(self))
    
    @staticmethod
    def fromRaw(filename):
        g=Grid.fromFile(filename)
        #  we assume tris and an axis aligned grid
        g.center=np.reshape(g.center,(-1,3))
        #  keep unique vertices only and sort first on x then on y coordinate
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
        #  keep just the z-values
        g.center=np.array([c[2] for c in verts],dtype=np.single).reshape(nx,ny)
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

    def zeroedge(self):
        c = self.center
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
        return self.center
    
    def avalanche_old(self, delta, prob, numexpr):
        self.neighborgrid()
        if numexpr and numexpr_available:
            c=self.center
            up=self.up
            down=self.down
            right=self.right
            left=self.left
            
            dp=np.random.random_sample(c.shape, dtype=np.single)
            
            du=ne.evaluate('up - c')
            dd=ne.evaluate('down - c')
            dl=ne.evaluate('left - c')
            dr=ne.evaluate('right - c')
            
            c=ne.evaluate('c+where(((abs(du) - delta)>0) & (dp<prob), du/2, du)+where(((abs(dd) - delta)>0) & (dp<prob), dd/2, dd)+where(((abs(dr) - delta)>0) & (dp<prob), dr/2, dr)+where(((abs(dl) - delta)>0) & (dp<prob), dl/2, dl)')

        else:
            c=self.center
            dp=np.random.random_sample(c.shape, dtype=np.single) < prob
            
            du=self.up - c
            dd=self.down - c
            dl=self.left - c
            dr=self.right - c
            
            ddd=(np.absolute(du) - delta)>0
            np.divide(du,2,du,where=ddd)
            np.add(c,du,c,where=dp)
            
            ddd=(np.absolute(dd) - delta)>0
            np.divide(dd,2,dd,where=ddd)
            np.add(c,dd,c,where=dp)
            
            ddd=(np.absolute(dl) - delta)>0
            np.divide(dl,2,dl,where=ddd)
            np.add(c,dl,c,where=dp)
            
            ddd=(np.absolute(dr) - delta)>0
            np.divide(dr,2,dr,where=ddd)
            np.add(c,dr,c,where=dp)
        
        
    def diffuse_old(self, Kd, numexpr=False):
        self.neighborgrid()
        center = self.center
        up = self.up
        down = self.down
        right = self.right
        left = self.left
        if(numexpr and numexpr_available):
            center = ne.evaluate('center+Kd*(up+down+left+right-4.0*center)')
        else:
            center += Kd*(up+down+left+right-4.0*center)
            
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
  t = time()
  for g in range(args.iterations):
    if args.Kd > 0:
      grid.diffuse(args.Kd, args.usenumexpr)
    if args.Kh > 0 and args.Kp > rand():
      grid.avalanche(args.Kh, args.usenumexpr)
  t = time()-t
  print("\nElapsed time: %.1f seconds.\n"%t, file=sys.stderr)
  if args.analyze:
    print('\nstatistics of the output grid:\n\n', grid.analyze(), file=sys.stderr, sep='')

  if not args.timingonly:
    if args.rawout:
      grid.toRaw(args.outfile)
    else:
      grid.toFile(args.outfile)
