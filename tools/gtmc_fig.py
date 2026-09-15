# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | wangyudong@buaa.edu.cn
"""Render real GTMC central slice and translucent 3D field demonstration."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cmcrameri import cm
import pyvista as pv
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--geometry',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);b=out;d=np.load(args.data)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,2,figsize=(10.5,10),layout='constrained')
speed=np.linalg.norm(d['Us'],axis=-1);solid=d['ms']==0
for ax,a,title,palette,limits,ticks,label in [(axs[0],speed,'Velocity magnitude',cm.batlow,(0,150),[0,30,60,90,120,150],'|U| (m/s)'),(axs[1],d['Ts'],'Temperature',cm.lajolla,(290,2300),[300,700,1100,1500,1900,2300],'T (K)')]:
 im=ax.pcolormesh(d['X']*1000,d['Y']*1000,np.ma.array(a,mask=solid),cmap=palette,vmin=limits[0],vmax=limits[1],rasterized=True)
 ax.set_facecolor('#e2e7ea');ax.set_aspect('equal');ax.set_title(title,loc='left',fontweight='bold',pad=12);ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)');ax.spines[['top','right']].set_visible(False)
 cb=fig.colorbar(im,ax=ax,orientation='horizontal',pad=.06,shrink=.98,aspect=26,ticks=ticks);cb.set_label(label)
fig.suptitle('HUNDUN-FLOW  /  GTMC reacting flow',x=.055,ha='left',fontweight='bold',fontsize=19)
fig.text(.055,-.025,'Step %d   •   t = %.6f ms   •   z = %.4f mm   |   Gray: IBM solid'%(d['step'],d['time']*1000,d['slice_z']*1000),fontsize=10,color='#52616b')
for ext in ('png','svg','pdf'):fig.savefig(out/('gtmc-mid.'+ext),dpi=220,bbox_inches='tight',pad_inches=.16)
fig.canvas.draw();plt.close(fig)
# True 3D sampled field: physical coordinates; point samples every second cell.
grid=pv.RectilinearGrid(d['x']*1000,d['y']*1000,d['z']*1000)
grid['T']=d['T'].ravel(order='C');grid['fluid']=d['mask'].ravel(order='C').astype(float)
grid['U']=d['U'].reshape(-1,3);grid['speed']=np.linalg.norm(d['U'],axis=-1).ravel(order='C')
vertices=np.array([list(map(float,line.split()[1:])) for line in args.geometry.read_text().splitlines() if line.strip().startswith('vertex ')])*1000
assert len(vertices)%3==0 and np.isfinite(vertices).all()
facets=np.column_stack((np.full(len(vertices)//3,3),np.arange(len(vertices)).reshape(-1,3)))
geom=pv.PolyData(vertices,facets).clean();assert geom.n_cells==len(vertices)//3
# Complete physical shell; transparency reveals the internal field.
pl=pv.Plotter(off_screen=True,window_size=(1700,1500));pl.set_background('white')
pl.add_mesh(geom,color='#8a9aab',opacity=.09,smooth_shading=True,show_scalar_bar=False)
hot=grid.contour([1400,1800],scalars='T')
pl.add_mesh(hot,scalars='T',cmap=cm.lajolla,clim=(290,2300),opacity=.45,smooth_shading=True,show_scalar_bar=False)
# Spatially distributed seeds in actual fluid cells throughout the swirler.
flow=grid.threshold(.5,scalars='fluid',preference='point',all_scalars=True)
seed_points=[];seed_levels=[]
for height in (-44,-32,-18,-7,2):
 j=int(np.argmin(abs(d['y']*1000-height)))
 speed_layer=np.linalg.norm(d['U'][:,j,:,:],axis=-1)
 kk,ii=np.nonzero((d['mask'][:,j,:]==1)&(speed_layer>.5))
 candidates=np.column_stack([d['x'][ii]*1000,np.full(len(ii),d['y'][j]*1000),d['z'][kk]*1000])
 candidates=candidates[flow.find_containing_cell(candidates)>=0]
 selected=[]
 if len(candidates):
  selected=[0];distance=np.sum((candidates-candidates[0])**2,axis=1)
  while len(selected)<90 and distance.max()>2.0**2:
   q=int(np.argmax(distance));selected.append(q)
   distance=np.minimum(distance,np.sum((candidates-candidates[q])**2,axis=1))
 seed_points.extend(candidates[selected].tolist())
 seed_levels.append(dict(y_mm=float(d['y'][j]*1000),seeds=len(selected)))
assert len(seed_points)>200
seeds=pv.PolyData(np.array(seed_points))
lines=flow.streamlines_from_source(seeds,vectors='U',integration_direction='both',max_steps=1600,initial_step_length=.12,max_step_length=.25,terminal_speed=.2,max_length=180)
if lines.n_points:
 lines['speed']=np.linalg.norm(lines['U'],axis=1)
 pl.add_mesh(lines.tube(radius=.10),scalars='speed',cmap=cm.batlow,clim=(0,150),show_scalar_bar=False,opacity=.9)
pl.camera_position=[(235,145,305),(0,49,0),(0,1,0)];pl.camera.parallel_projection=True;pl.camera.parallel_scale=124
pl.enable_depth_peeling(number_of_peels=16,occlusion_ratio=0)
pl.add_text('HUNDUN-FLOW',position=(85,1410),font_size=30,color='#1b2c40',font='arial')
pl.add_text('GTMC / Reacting flow',position=(88,1360),font_size=19,color='#4b6276',font='arial')
pl.add_text('t = %.6f ms  |  Step %d'%(d['time']*1000,d['step']),position=(88,1315),font_size=13,color='#647381',font='arial')
pl.add_text('Temperature: 1400 / 1800 K isosurfaces\nVelocity: distributed swirler-channel streamlines\nGeometry: full translucent STL shell',position=(85,35),font_size=13,color='#52616b',font='arial')
pl.screenshot(str(b/'3d.png'));pl.close()
# Compose honest shared color scales alongside the rendered 3D scene.
from PIL import Image
fig=plt.figure(figsize=(11.4,10.8));ax=fig.add_axes([0,0,1,1]);ax.imshow(Image.open(b/'3d.png'));ax.axis('off')
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
for pos,limits,palette,label,ticks in [([.80,.38,.014,.31],(290,2300),cm.lajolla,'T (K)',[300,700,1100,1500,1900,2300]),([.88,.38,.014,.31],(0,150),cm.batlow,'|U| (m/s)',[0,30,60,90,120,150])]:
 cax=fig.add_axes(pos);cb=fig.colorbar(ScalarMappable(norm=Normalize(*limits),cmap=palette),cax=cax,ticks=ticks);cb.ax.set_title(label,fontsize=10,pad=12);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=9,length=2)
for ext in ('png','svg','pdf'):fig.savefig(out/('gtmc-3d.'+ext),dpi=190,bbox_inches='tight',pad_inches=.04)
plt.close(fig)
(b/'render.json').write_text(json.dumps(dict(geometry=str(args.geometry),isosurfaces_K=[1400,1800],streamline_points=int(lines.n_points),streamline_cells=int(lines.n_cells),seed_count=len(seed_points),seed_levels=seed_levels,geometry_view="complete shell",temperature_view="isosurfaces",geometry_triangles=int(geom.n_cells),field_stride=2,center_slice='original cell resolution',color_ranges=dict(temperature_K=[290,2300],speed_m_s=[0,150]),python_versions=dict(pyvista=pv.__version__,matplotlib=matplotlib.__version__,numpy=np.__version__),qa='Review image layout after each render'),indent=2)+'\n')
print('figures rendered',lines.n_points)
