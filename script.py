import bpy
import bmesh
from mathutils import Vector

import numpy as np
from scipy.integrate import odeint
from numpy.random import choice

import math
from itertools import product

from sympy.physics.hydrogen import Psi_nlm
from sympy import Symbol, simplify
from sympy.vector import gradient, CoordSys3D
from sympy.utilities.lambdify import lambdify
from sympy.functions.elementary.complexes import Abs

# Parameters of hydrogen orbital
n = 5
l = 4
m = 1   

# Number of points that will be sampled
num_points = 100

# Some animation parameters
num_frames = 1800 # Number of frames
end_time = 100 # The physical time that passes in num_frames
framestep = 10 # Every framestep-th frame gets a keyframe, the rest is interpolation

# Discretization of space
r = np.linspace(0, 100, 100)
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)



# Conversion between spherical and cartesian coordinates

def to_spherical(x,y,z):
    XsqPlusYSq = x**2 +y**2
    r = np.sqrt(XsqPlusYSq+z**2)
    theta = np.arctan2(np.sqrt(XsqPlusYSq), z)
    phi = np.arctan2(y,x)
    return r, theta, phi


def to_cartesian(co):
    r, theta, phi = co[0], co[1], co[2]
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

######################################################
# Sampling the wavefunction and creating the meshes
######################################################

locations = np.array(list(product(r, theta, phi)))

c = CoordSys3D('c', transformation='spherical', variable_names=('r','theta', 'phi'))

prob = Abs(Psi_nlm(n=n, l=l, m=m, r=c.r, phi=c.phi, theta=c.theta))**2
prob = lambdify((c.r, c.theta, c.phi), prob)
prob = np.vectorize(prob)

print('Calculating probabilities')
probabilities = prob(locations[:,0], locations[:,1], locations[:,2])
probabilities = probabilities/np.sum(probabilities)

print('Sampling probabilities')

inds = choice(np.arange(len(locations)), size=num_points, p=probabilities, replace=False)

rs = locations[inds, 0]
thetas = locations[inds, 1]
phis = locations[inds, 2]

xs = rs * np.sin(thetas) * np.cos(phis)
ys = rs * np.sin(thetas) * np.sin(phis)
zs = rs * np.cos(thetas)

# Preparing the 'orb' that will represent the electron locations

bpy.ops.mesh.primitive_ico_sphere_add(radius=0.56, subdivisions=1)
orig_cube = bpy.context.active_object
bpy.ops.object.modifier_add(type='SUBSURF')
orig_cube.modifiers['Subdivision'].levels = 0
orig_cube.modifiers['Subdivision'].render_levels = 2

# Make the mesh that holds all the positions

bpy.ops.mesh.primitive_plane_add()
o = bpy.context.active_object
me = o.data
bm = bmesh.new()
for i in range(num_points):
    bm.verts.new().co=Vector((xs[i],ys[i], zs[i]))
bm.to_mesh(me)
bm.free()

# Make the 'position' mesh the parent of the 'orb' and use 
# verts instancing, so each vertex of the 'position' mesh
# gets rendered as the 'orb'

orig_cube.parent = o
o.instance_type = 'VERTS'

######################################################
# Creating the animation data and animating
######################################################

# Preparing the mesh

mesh = o.data
action = bpy.data.actions.new("MeshAnimation")

mesh.animation_data_create()
mesh.animation_data.action = action

data_path = "vertices[%d].co"

# Using sympy's hydrogen wavefunction and 
# defining a spherical coordinate system in sympy

c = CoordSys3D('c', transformation='spherical', variable_names=('r','theta', 'phi'))

psi = Psi_nlm(n=n, l=l, m=m, r=c.r, phi=c.phi, theta=c.theta)
# This is (almost) the right hand side of the guiding equation
gr = gradient(psi)/psi
# Split the vector into its components
gr1 = simplify(gr.dot(c.i)) # Simplify helps with making the 
gr2 = simplify(gr.dot(c.j)) # numerical solution more stable
gr3 = simplify(gr.dot(c.k))
# Sympy is for symbolic manipulation, but we want numerical functions
# -> lambdify turns sympy functions into function we can plug numbers into
f1 = lambdify((c.r, c.theta, c.phi), gr1)
f2 = lambdify((c.r, c.theta, c.phi), gr2)
f3 = lambdify((c.r, c.theta, c.phi), gr3)


frames = np.arange(0, num_frames, framestep)
t = np.linspace(0, end_time, len(frames))

# This is the right hand side of the guiding equation
def dSdt(S, t):
    return [f1(S[0], S[1], S[2]).imag,
            f2(S[0], S[1], S[2]).imag,
            f3(S[0], S[1], S[2]).imag]

# For each vertex in the mesh, solve the guiding equation numerically
# and then animate the vertex with the solution
for  v in mesh.vertices:
    S0 = [*to_spherical(*v.co.copy())]
    S = odeint(dSdt, S0, t)
    values = to_cartesian([S[:,0],S[:,1],S[:,2]])
    # Progress
    print(v.index/num_points * 100 ,"%")
    
    fcu = [action.fcurves.new(data_path % v.index, index=i) for i in range(3)]
    for i,fc in enumerate(fcu):
        fc.keyframe_points.add(count=len(t))
        fc.keyframe_points.foreach_set('co', [x for co in zip(frames, values[i]) for x in co])

# Set interpolation type of keyframes to Bezier curves
area = bpy.context.area.type
bpy.context.area.type = 'DOPESHEET_EDITOR'
bpy.ops.action.interpolation_type(type='BEZIER')
bpy.context.area.type = area


######################################################
# Setting up material
######################################################


# Create and assign material
mat_name = "MyMat"
# Test if material exists
# If it does not exist, create it:
mat = (bpy.data.materials.get(mat_name) or 
       bpy.data.materials.new(mat_name))
       
# Create the material in the Node editor
# In summary: The materials is invisible
# for y<0 and z>0
mat.use_nodes = True
nodes = mat.node_tree.nodes
nodes.clear()
node_bsdf = nodes.new(type='ShaderNodeBsdfPrincipled')
node_objectinfo = nodes.new(type='ShaderNodeObjectInfo')
node_sep = nodes.new(type='ShaderNodeSeparateXYZ')
node_mathY = nodes.new(type='ShaderNodeMath')
node_mathY.operation = 'GREATER_THAN'
node_mathY.inputs[1].default_value = 0.0
node_mathZ = nodes.new(type='ShaderNodeMath')
node_mathZ.operation = 'GREATER_THAN'
node_mathZ.inputs[0].default_value = 0.0
node_add = nodes.new(type='ShaderNodeMath')
node_add.operation = 'ADD'
node_add.use_clamp = True
node_output = nodes.new(type='ShaderNodeOutputMaterial')

# link nodes
links = mat.node_tree.links
link = links.new(node_objectinfo.outputs[0], node_sep.inputs[0])
link = links.new(node_sep.outputs[1], node_mathY.inputs[0]) # y position
link = links.new(node_sep.outputs[2], node_mathZ.inputs[1]) # z positiob
link = links.new(node_mathY.outputs[0], node_add.inputs[0])
link = links.new(node_mathZ.outputs[0], node_add.inputs[1])
link = links.new(node_add.outputs[0], node_bsdf.inputs[19]) #link to alpha
link = links.new(node_bsdf.outputs[0], node_output.inputs[0])

mat.blend_method = 'CLIP'
mat.shadow_method = 'CLIP'
# assign material
if orig_cube.data.materials:
    # assign to 1st material slot
    orig_cube.data.materials[0] = mat
else:
    # no slots
    orig_cube.data.materials.append(mat)
