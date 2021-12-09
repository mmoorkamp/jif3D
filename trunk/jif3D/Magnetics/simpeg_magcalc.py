"""
Forward Simulation of Total Magnetic Intensity Data
===================================================

Here we use the module *SimPEG.potential_fields.magnetics* to predict magnetic
data for a magnetic susceptibility model. We simulate the data on a tensor mesh.
For this tutorial, we focus on the following:

    - How to define the survey
    - How to predict magnetic data for a susceptibility model
    - How to include surface topography
    - The units of the physical property model and resulting data


"""

#########################################################################
# Import Modules
# --------------
#

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

from discretize import TensorMesh
from discretize.utils import mkvc
from SimPEG.utils import plot2Ddata, model_builder, surface2ind_topo
from SimPEG import maps
from SimPEG.potential_fields import magnetics
from SimPEG import utils, data

save_file = False

# sphinx_gallery_thumbnail_number = 2


#############################################
# Topography
# ----------
#
# Surface topography is defined as an (N, 3) numpy array. We create it here but
# topography could also be loaded from a file.
#

[x_topo, y_topo] = np.meshgrid(np.linspace(-200, 200, 41), np.linspace(-200, 200, 41))
z_topo = np.zeros(x_topo.shape)
x_topo, y_topo, z_topo = mkvc(x_topo), mkvc(y_topo), mkvc(z_topo)
xyz_topo = np.c_[x_topo, y_topo, z_topo]

#############################################
# Defining the Survey
# -------------------
#
# Here, we define survey that will be used for the simulation. Magnetic
# surveys are simple to create. The user only needs an (N, 3) array to define
# the xyz locations of the observation locations, the list of field components
# which are to be modeled and the properties of the Earth's field.
#

# Define the observation locations as an (N, 3) numpy array or load them.
y = np.arange(-20.0, 20, 0.2)
x = np.ones(y.shape) * 10
z = np.ones(x.shape)  # Flight height 1 m above surface.
receiver_locations = np.c_[x, y, z]

# Define the component(s) of the field we want to simulate as a list of strings.
# Here we simulation total magnetic intensity data.
components = ["by","bx","bz", "tmi"]

# Use the observation locations and components to define the receivers. To
# simulate data, the receivers must be defined as a list.
receiver_list = magnetics.receivers.Point(receiver_locations, components=components)

receiver_list = [receiver_list]

# Define the inducing field H0 = (intensity [nT], inclination [deg], declination [deg])
inclination = 15
declination = 10
strength = 50000
inducing_field = (strength, inclination, declination)

source_field = magnetics.sources.SourceField(
    receiver_list=receiver_list, parameters=inducing_field
)

# Define the survey
survey = magnetics.survey.Survey(source_field)


#############################################
# Defining a Tensor Mesh
# ----------------------
#
# Here, we create the tensor mesh that will be used for the forward simulation.
#

dh = 5.0
hx = [(dh, 5, -1.3), (dh, 40), (dh, 5, 1.3)]
hy = [(dh, 5, -1.3), (dh, 40), (dh, 5, 1.3)]
hz = [(dh, 5, -1.3), (dh, 15)]
mesh = TensorMesh([hx, hy, hz], "CCN")


#############################################
# Defining a Susceptibility Model
# -------------------------------
#
# Here, we create the model that will be used to predict magnetic data
# and the mapping from the model to the mesh. The model
# consists of a susceptible sphere in a less susceptible host.
#

# Define susceptibility values for each unit in SI
background_susceptibility = 0.0
sphere_susceptibility = 0.2

# Find cells that are active in the forward modeling (cells below surface)
ind_active = surface2ind_topo(mesh, xyz_topo)

# Define mapping from model to active cells
nC = int(ind_active.sum())
model_map = maps.IdentityMap(nP=nC)  # model is a vlue for each active cell

# Define model. Models in SimPEG are vector arrays
model = background_susceptibility * np.ones(ind_active.sum())
ind_sphere = model_builder.getIndicesBlock([20,10,0], [-20,-10,-5], mesh.gridCC)
model[ind_sphere] = sphere_susceptibility

# Plot Model
fig = plt.figure(figsize=(9, 4))

plotting_map = maps.InjectActiveCells(mesh, ind_active, np.nan)
ax1 = fig.add_axes([0.1, 0.12, 0.73, 0.78])
mesh.plotSlice(
    plotting_map * model,
    normal="Y",
    ax=ax1,
    ind=int(mesh.nCy / 2),
    grid=True,
    clim=(np.min(model), np.max(model)),
)
ax1.set_title("Model slice at y = 0 m")
ax1.set_xlabel("x (m)")
ax1.set_ylabel("z (m)")

ax2 = fig.add_axes([0.85, 0.12, 0.05, 0.78])
norm = mpl.colors.Normalize(vmin=np.min(model), vmax=np.max(model))
cbar = mpl.colorbar.ColorbarBase(ax2, norm=norm, orientation="vertical")
cbar.set_label("Magnetic Susceptibility (SI)", rotation=270, labelpad=15, size=12)

plt.show()



fig = plt.figure(figsize=(9, 4))

plotting_map = maps.InjectActiveCells(mesh, ind_active, np.nan)
ax1 = fig.add_axes([0.1, 0.12, 0.73, 0.78])
mesh.plotSlice(
    plotting_map * model,
    normal="X",
    ax=ax1,
    ind=int(mesh.nCy / 2),
    grid=True,
    clim=(np.min(model), np.max(model)),
)
ax1.set_title("Model slice at x = 0 m")
ax1.set_xlabel("y (m)")
ax1.set_ylabel("z (m)")

ax2 = fig.add_axes([0.85, 0.12, 0.05, 0.78])
norm = mpl.colors.Normalize(vmin=np.min(model), vmax=np.max(model))
cbar = mpl.colorbar.ColorbarBase(ax2, norm=norm, orientation="vertical")
cbar.set_label("Magnetic Susceptibility (SI)", rotation=270, labelpad=15, size=12)

plt.show()



###################################################################
# Simulation: TMI Data for a Susceptibility Model
# -----------------------------------------------
#
# Here we demonstrate how to predict magnetic data for a magnetic
# susceptibility model using the integral formulation.
#

# Define the forward simulation. By setting the 'store_sensitivities' keyword
# argument to "forward_only", we simulate the data without storing the sensitivities
simulation = magnetics.simulation.Simulation3DIntegral(
    survey=survey,
    mesh=mesh,
    modelType="susceptibility",
    chiMap=model_map,
    actInd=ind_active,
    store_sensitivities="forward_only",
)

# Compute predicted data for a susceptibility model
dpred = simulation.dpred(model)
np.savetxt("magx.simpeg",np.transpose([dpred[::4]]))
np.savetxt("magy.simpeg",np.transpose([dpred[1::4]]))
np.savetxt("magz.simpeg",np.transpose([-dpred[2::4]]))
np.savetxt("magt.simpeg",np.transpose([dpred[3::4]]))

jif3D = np.genfromtxt("/home/max/workspace/jif3D/magt.out")
jif3Dx = np.genfromtxt("/home/max/workspace/jif3D/magx.out")
jif3Dy = np.genfromtxt("/home/max/workspace/jif3D/magy.out")
jif3Dz = np.genfromtxt("/home/max/workspace/jif3D/magz.out")


      
def MagBox(xpos, width, depth, height, magnetization):
    m = width / 2.0
    factor = magnetization / (2 * np.pi)
    Bz = factor * (np.arctan2(xpos+m,depth) - np.arctan2(xpos-m,depth) - np.arctan2(xpos+m,depth+height) + np.arctan2(xpos-m,depth+height))
    return Bz

bbox = MagBox(x, 20, z[0], 5, strength * sphere_susceptibility)

plt.figure()
plt.subplot(221)
plt.plot(dpred[::4])
plt.plot(jif3Dx[:,1])
plt.title("Bx")

plt.subplot(222)
plt.plot(dpred[1::4])
plt.plot(jif3Dy[:,1])
plt.title("By")

plt.subplot(223)
plt.plot(-dpred[2::4])
plt.plot(jif3Dz[:,1])
plt.title("Bz")

#plt.plot(bbox)
plt.show()