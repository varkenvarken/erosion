bl_info = {
    "name": "Erosion",
    "author": "Michel Anders (varkenvarken)",
    "version": (0, 0, 2),
    "blender": (2, 69, 0),
    "location": "View3D > Object > Erode",
    "description": "Apply various kinds of erosion to a mesh",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Mesh"}

from random import random as rand
from math import tan, radians
import bpy
from bpy.props import FloatProperty, IntProperty, BoolProperty, EnumProperty, StringProperty    
from .erode import Grid
from .stats import Stats
from .utils import numexpr_available

def availableVertexGroupsOrNone(self, context):
    groups = [ ('None', 'None', 'None', 1) ]
    return groups + [(name, name, name, n+1) for n,name in enumerate(context.active_object.vertex_groups.keys())]
    
class Erode(bpy.types.Operator):
    bl_idname = "mesh.erode"
    bl_label = "Erode"
    bl_options = {'REGISTER', 'UNDO', 'PRESET'}

    Iterations = IntProperty(name="Iterations", description="Number of iterations", default=10, min=0, soft_max=100)
    
    Kd = FloatProperty(name="Kd", description="Thermal diffusion rate", default=0.005, min=0)

    Kt = FloatProperty(name="Kt", description="Maximum stable talus angle", default=radians(60), min=0, max=radians(90), subtype='ANGLE')

    Kr   = FloatProperty(name="Rain amount"      , description="Rain amount                                ", default=1   , min=0, soft_max=1)
    Kv   = FloatProperty(name="Rain variance"    , description="Rain variance (0 is constant, 1 is uniform)", default=0   , min=0, max=1)
    userainmap = BoolProperty(name="Use rain map", description="Use active vertex group as a rain map"      , default=False)
    
    Ks   = FloatProperty(name="Soil solubility"  , description="Soil solubity"                              , default=0.07, min=0, soft_max=1)
    Kdep = FloatProperty(name="Deposition rate"  , description="Sediment deposition rate"                   , default=0.07, min=0, soft_max=1)
    Kc   = FloatProperty(name="Carrying capacity", description="Base sediment carrying capacity"            , default=0.9 , min=0, soft_max=1)
    Ka   = FloatProperty(name="Slope dependence" , description="Slope dependence of carrying capacity"      , default=1.0 , min=0, soft_max=2)

    numexpr = BoolProperty(name="Numexpr", description="Use numexpr module (if available)", default=True)

    Pd = FloatProperty(name="Pd", description="Diffusion probability"    , default=0.1, min=0, max=1)
    Pa = FloatProperty(name="Pa", description="Avalanche probability"    , default=0.1, min=0, max=1)
    Pw = FloatProperty(name="Pw", description="Water erosion probability", default=1  , min=0, max=1)
    
    smooth = BoolProperty(name="Smooth", description="Set smooth shading", default=True)

    showiterstats = BoolProperty(name="Iteration Stats", description="Show iteraration statistics", default=False)
    showmeshstats = BoolProperty(name="Mesh Stats"     , description="Show mesh statistics"       , default=False)

    stats = Stats()
    counts= {}

    # add poll function to restrict action to mesh object in object mode
    
    def execute(self, context):
        ob = context.active_object
        me = ob.data
        self.stats.reset()
        vg=ob.vertex_groups.active
        g = Grid.fromBlenderMesh(me, vg)
            
        me = bpy.data.meshes.new(me.name)

        self.counts['diffuse']=0
        self.counts['avalanche']=0
        self.counts['water']=0
        for i in range(self.Iterations):
            if self.Kd > 0.0 and rand() < self.Pd:
                g.diffuse(self.Kd, self.numexpr)
                self.counts['diffuse']+=1
            if self.Kt < radians(90) and rand() < self.Pa:
                # since dx and dy are scaled to 1, tan(Kt) is the height for a given angle
                g.avalanche(tan(self.Kt), self.numexpr)
                self.counts['avalanche']+=1
            if self.Kr > 0 and rand() < self.Pw:
                g.fluvial_erosion(self.Kr, self.Kv, self.userainmap, self.Kc, self.Ks, self.Kdep, self.Ka, 0,0,0,0, self.numexpr)
                self.counts['water']+=1

        g.toBlenderMesh(me)
        ob.data = me
        if vg:
            for row in range(g.rainmap.shape[0]):
                for col in range(g.rainmap.shape[1]):
                    i = row * g.rainmap.shape[1] + col
                    vg.add([i],g.rainmap[row,col],'ADD')
            
        if self.smooth:
            bpy.ops.object.shade_smooth()
        self.stats.time()
        self.stats.memory()
        if self.showmeshstats:
            self.stats.meshstats = g.analyze()
            
        return {'FINISHED'}
    
    def draw(self,context):
        layout = self.layout

        layout.operator('screen.repeat_last', text="Repeat", icon='FILE_REFRESH' )

        layout.prop(self, 'Iterations')
        
        box = layout.box()
        box.label("Thermal erosion")
        box.prop(self, 'Kd')
        box.prop(self, 'Kt')
        
        box = layout.box()
        box.label("Hydraulic erosion")
        box.prop(self, 'Ks')
        box.prop(self, 'Kc')
        box.prop(self, 'Kdep')
        box.prop(self, 'Kr')
        box.prop(self, 'Kv')
        box2 = box.box()
        box2.prop(self, 'userainmap')
        box2.enabled = context.active_object.vertex_groups.active is not None
        box.prop(self, 'Ka')
        
        box = layout.box()
        box.label("Probabilities")
        box.prop(self, 'Pd')
        box.prop(self, 'Pa')
        box.prop(self, 'Pw')
        
        layout.prop(self,'smooth')
        
        if numexpr_available:
          layout.prop(self, 'numexpr')
        else:
          box = layout.box()
          box.alert=True
          box.label("Numexpr not available. Will slow down large meshes")

        box = layout.box()
        box.prop(self,'showiterstats')
        if self.showiterstats:
            row = box.row()
            col1 = row.column()
            col2 = row.column()
            col1.label("Time"); col2.label("%.1f s"%self.stats.elapsedtime)
            if self.stats.memstats_available:
                col1.label("Memory"); col2.label("%.1f Mb"%(self.stats.maxmem/(1024.0*1024.0)))
            col1.label("Diffusions"); col2.label("%d"% self.counts['diffuse'])
            col1.label("Avalanches"); col2.label("%d"% self.counts['avalanche'])
            col1.label("Water movements"); col2.label("%d"% self.counts['water'])
        box = layout.box()
        box.prop(self,'showmeshstats')
        if self.showmeshstats:
            row = box.row()
            col1 = row.column()
            col2 = row.column()
            for line in self.stats.meshstats.split('\n'):
                label, value = line.split(':')
                col1.label(label)
                col2.label(value)

def menu_func_erode(self, context):
    self.layout.operator(Erode.bl_idname, text="Erode",
                         icon='RNDCURVE')

def register():
    bpy.utils.register_class(Erode)
    bpy.types.VIEW3D_MT_object.append(menu_func_erode)


def unregister():
    bpy.types.VIEW3D_MT_object.remove(menu_func_erode)
    bpy.utils.unregister_class(Erode)