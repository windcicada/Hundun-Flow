#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Read the production legacy output using VTK's independent binary reader."""
from pathlib import Path
import sys
from vtkmodules.vtkIOLegacy import vtkStructuredGridReader
from vtkmodules.vtkFiltersCore import vtkAppendFilter
from vtkmodules.vtkFiltersFlowPaths import vtkStreamTracer
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkCommonCore import vtkPoints
root=Path(sys.argv[1]);ranks=int(sys.argv[2])
index=(root/'solution.visit').read_text().splitlines()
assert index[0]=='!NBLOCKS '+str(ranks)
assert [float(s.split()[1]) for s in index if s.startswith('!TIME')]==[.125,.3]
for step in (1,2):
    combined=vtkAppendFilter();combined.MergePointsOn()
    owned=0;cells=0
    for rank in range(ranks):
        path=root/('s%dr%d.vtk'%(step,rank))
        assert b'BINARY\nDATASET STRUCTURED_GRID\n' in path.read_bytes()[:200]
        reader=vtkStructuredGridReader();reader.SetFileName(str(path));reader.ReadAllScalarsOn();reader.ReadAllVectorsOn();reader.Update()
        grid=reader.GetOutput();assert grid.GetNumberOfPoints()>0
        data=grid.GetPointData()
        u=data.GetArray('Velocity');p=data.GetArray('Pressure');pg=data.GetArray('PressureGauge');t=data.GetArray('Temperature');ghost=data.GetArray('vtkGhostType')
        assert all(a is not None for a in (u,p,pg,t,ghost))
        assert u.GetDataTypeAsString()=='double' and u.GetNumberOfComponents()==3
        for i in range(grid.GetNumberOfPoints()):
            x,y,z=grid.GetPoint(i)
            assert abs(u.GetTuple3(i)[0]-(1+x))<1e-12
            assert abs(t.GetTuple1(i)-(300+2*x+3*y+5*z))<1e-11
            assert abs(pg.GetTuple1(i)-(x+2*y+4*z))<1e-12
            assert abs(p.GetTuple1(i)-(101325+x+2*y+4*z))<3e-11
            owned+=ghost.GetTuple1(i)==0
        cells+=grid.GetNumberOfCells();combined.AddInputData(grid)
    assert owned==17*11*7 and cells==16*10*6,(owned,cells)
    combined.Update();merged=combined.GetOutput()
    merged.GetPointData().SetActiveVectors('Velocity')
    seeds=vtkPolyData();points=vtkPoints();points.InsertNextPoint(.05,.5,.5);seeds.SetPoints(points)
    tracer=vtkStreamTracer();tracer.SetInputData(merged);tracer.SetSourceData(seeds)
    tracer.SetMaximumPropagation(2);tracer.SetInitialIntegrationStep(.01);tracer.SetIntegrationDirectionToForward();tracer.Update()
    line=tracer.GetOutput();assert line.GetNumberOfPoints()>10
    assert line.GetBounds()[1]>.9,line.GetBounds()
print('VTK binary reader: stretched centres, affine fields, overlap, time index and cross-partition streamline PASS')
