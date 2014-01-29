bl_info = {
    "name": "Erosion",
    "author": "Michel Anders (varkenvarken)",
    "version": (0, 0, 1),
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

class Erode(bpy.types.Operator):
    bl_idname = "mesh.erode"
    bl_label = "Erode"
    bl_options = {'REGISTER', 'UNDO', 'PRESET'}

    Iterations = IntProperty(name="Iterations", description="Number of iterations", default=10, min=0, soft_max=100)
    Kd = FloatProperty(name="Kd", description="Thermal diffusion rate", default=0.005, min=0)

    Kt = FloatProperty(name="Kt", description="Maximum stable talus angle", default=radians(60), min=0, max=radians(90), subtype='ANGLE')

    Kr   = FloatProperty(name="Rain amount"      , description="Amount of rain per iteration"         , default=0.1 , min=0, soft_max=2)
    Ks   = FloatProperty(name="Soil solubility"  , description="Soil solubity"                        , default=0.07, min=0, soft_max=1)
    Kdep = FloatProperty(name="Deposition rate"  , description="Sediment deposition rate"             , default=0.07, min=0, soft_max=1)
    Kc   = FloatProperty(name="Carrying capacity", description="Base sediment carrying capacity"      , default=0.9 , min=0, soft_max=1)
    Ka   = FloatProperty(name="Slope dependence" , description="Slope dependence of carrying capacity", default=1.0 , min=0, soft_max=2)

    numexpr = BoolProperty(name="Numexpr", description="Use numexpr module (if available)", default=True)

    Pd = FloatProperty(name="Pd", description="Diffusion probability"    , default=0.1, min=0, max=1)
    Pa = FloatProperty(name="Pa", description="Avalanche probability"    , default=0.1, min=0, max=1)
    Pw = FloatProperty(name="Pw", description="Water erosion probability", default=1  , min=0, max=1)
    
    smooth = BoolProperty(name="Smooth", description="Set smooth shading", default=True)

    showstats = BoolProperty(name="Stats", description="Show statistics", default=False)

    stats = Stats()
    counts= {}

    # add poll function to restrict action to mesh object in object mode
    
    def execute(self, context):
        ob = context.active_object
        me = ob.data
        self.stats.reset()
        g = Grid.fromBlenderMesh(me)
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
                g.fluvial_erosion(self.Kr, self.Kc, self.Ks, self.Kdep, self.Ka, 0,0,0,0, self.numexpr)
                self.counts['water']+=1

        g.toBlenderMesh(me)
        ob.data = me
        if self.smooth:
            bpy.ops.object.shade_smooth()
        self.stats.time()
        self.stats.memory()
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
        box.prop(self,'showstats')
        if self.showstats:
            row = box.row()
            col1 = row.column()
            col2 = row.column()
            col1.label("Time"); col2.label("%.1f s"%self.stats.elapsedtime)
            if self.stats.memstats_available:
                col1.label("Memory"); col2.label("%.1f Mb"%(self.stats.maxmem/(1024.0*1024.0)))
            col1.label("Diffusions"); col2.label("%d"% self.counts['diffuse'])
            col1.label("Avalanches"); col2.label("%d"% self.counts['avalanche'])
            col1.label("Water movements"); col2.label("%d"% self.counts['water'])
        

def menu_func_erode(self, context):
    self.layout.operator(Erode.bl_idname, text="Erode",
                         icon='PLUGIN')

def register():
    bpy.utils.register_class(Erode)
    bpy.types.VIEW3D_MT_object.append(menu_func_erode)


def unregister():
    bpy.types.VIEW3D_MT_object.remove(menu_func_erode)
    bpy.utils.unregister_class(Erode)