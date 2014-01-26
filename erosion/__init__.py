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

import bpy
from bpy.props import FloatProperty, IntProperty, BoolProperty, EnumProperty, StringProperty    
from .erode import Grid
from .stats import Stats

class Erode(bpy.types.Operator):
    bl_idname = "mesh.erode"
    bl_label = "Erode"
    bl_options = {'REGISTER', 'UNDO'}

    Iterations = IntProperty(name="Iterations", description="Number of iterations", default=10, min=0, soft_max=100)
    Kd = FloatProperty(name="Kd", description="Thermal diffusion rate", default=0.01, min=0)
    numexpr = BoolProperty(name="Numexpr", description="Use numexpr module (if available)", default=True)
    
    stats = Stats()
    
    # add poll function to restrict action to mesh object in object mode
    
    def execute(self, context):
        ob = context.active_object
        me = ob.data
        self.stats.reset()
        g = Grid.fromBlenderMesh(me)
        me = bpy.data.meshes.new(me.name)
        for i in range(self.Iterations):
            g.diffuse(self.Kd, self.numexpr)
        g.toBlenderMesh(me)
        ob.data = me
        self.stats.time()
        self.stats.memory()
        return {'FINISHED'}
    
    def draw(self,context):
        layout = self.layout

        layout.prop(self, 'Iterations')
        layout.prop(self, 'Kd')
        layout.prop(self, 'numexpr')

        box = layout.box()
        box.label("Stats:")
        box.label("Time  : %.1f s"%self.stats.elapsedtime)
        if self.stats.memstats_available:
            box.label("Memory: %.1f Mb"%(self.stats.maxmem/(1024.0*1024.0)))

def menu_func_erode(self, context):
    self.layout.operator(Erode.bl_idname, text="Erode",
                         icon='PLUGIN')

def register():
    bpy.utils.register_class(Erode)
    bpy.types.VIEW3D_MT_object.append(menu_func_erode)


def unregister():
    bpy.types.VIEW3D_MT_object.remove(menu_func_erode)
    bpy.utils.unregister_class(Erode)