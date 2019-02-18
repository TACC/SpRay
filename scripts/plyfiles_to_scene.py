# ========================================================================== #
# Copyright (c) 2017-2018 The University of Texas at Austin.                 #
# All rights reserved.                                                       #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# A copy of the License is included with this software in the file LICENSE.  #
# If your copy does not contain the License, you may obtain a copy of the    #
# License at:                                                                #
#                                                                            #
#     https://www.apache.org/licenses/LICENSE-2.0                            #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  #
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
# ========================================================================== #

#!/usr/bin/python

import os
import sys
import argparse
# import subprocess
import re
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# references
#   https://matplotlib.org/tutorials/index.html#tutorials-colors
#   https://matplotlib.org/tutorials/colors/colormapnorms.html#sphx-glr-tutorials-colors-colormapnorms-py

def RgbaToRgbString(rgba):
  rgb = rgba[:-1] 
  rgb_str = ' '.join(str(e) for e in rgb)
  return rgb_str

# def valueToColor(value, value_min, value_max, colormap_name, is_logarithm):
def valueToColor(value, value_min, value_max, colormap_name):
  norm = mpl.colors.Normalize(vmin=value_min, vmax=value_max, clip=True)
  # if is_logarithm:
  #   norm = mpl.colors.LogNorm(vmin=value_min, vmax=value_max)
  # mapper = cm.ScalarMappable(norm=norm, cmap=cm.coolwarm)
  # mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(colormap_name))
  mapper = cm.ScalarMappable(norm=norm, cmap=colormap_name)
  rgba = mapper.to_rgba(value)
  return rgba

def getMinPoint(a, b):
  x = min(a[0], b[0])
  y = min(a[1], b[1])
  z = min(a[2], b[2])
  return [x, y, z]

def getMaxPoint(a, b):
  x = max(a[0], b[0])
  y = max(a[1], b[1])
  z = max(a[2], b[2])
  return [x, y, z]

def mergeBounds(model_bounds_list):
  min_point = model_bounds_list[0][0:3]
  max_point = model_bounds_list[0][3:6]

  for bounds_list in model_bounds_list:
    min_point = getMinPoint(bounds_list[0:3], min_point)
    max_point = getMaxPoint(bounds_list[3:6], max_point)

  return min_point + max_point

def addLights(outfile):
  with open(outfile, "w") as fout:
    fout.write("light diffuse-sphere .1 .1 .1\n")

def getDomainId(filename):
  # dns_r0180_c000_d0002_iso0.8.ply
  found = re.search('_d([0-9]+)_', filename)
  return int(found.group(1))

def getIsoValue(filename):
  # dns_r0180_c000_d0002_iso0.8.ply
  found = re.search('_iso([0-9.]+).ply', filename)
  return float(found.group(1))

def loadPlyFilesAndGenerateScene(args, plyfiles, loader_executable, outfile, bounds_info_dict):
  total_num_faces = 0 
  total_num_vertices = 0 

  for domain_id in sorted(plyfiles):

    model_bounds_list = []
    domain_bounds_str = None # domain bounds

    if bounds_info_dict:

      if str(domain_id) not in bounds_info_dict:
        print("[error] unknown domain ID found:" + str(domain_id))
        sys.exit(1)

      coordinates_list = bounds_info_dict[str(domain_id)]
      domain_bounds_str = ' '.join(coordinates_list)

    with open(outfile, "a") as fout:
      fout.write("######################\n")
      fout.write("# " + str(domain_id) + "\n")
      fout.write("DomainBegin\n")

      if domain_bounds_str:
        fout.write("DomainBounds " + domain_bounds_str + "\n")

    domain_num_faces = 0 
    domain_num_vertices = 0

    num_ply_models = len(plyfiles[domain_id])

    for ply, iso in plyfiles[domain_id]:
      # rgba = valueToColor(iso, args.contour_range[0], args.contour_range[1], args.colormap[0], args.logarithmic_colormap)
      rgba = valueToColor(iso, args.contour_range[0], args.contour_range[1], args.colormap[0])
      rgb_str = RgbaToRgbString(rgba)

      cmd = loader_executable + " " + ply + " out.tmp"
      print("[python] " + cmd)
      os.system(cmd)
      
      lines=[] 
      with open("./out.tmp", "r") as ftemp:
        for line in ftemp:
          lines.append(line)
      
      os.system("rm -f ./out.tmp")

      num_faces = 0 
      num_vertices = 0
      ply_bounds_str = None

      for line in lines:
        if line.startswith("face"):
          num_faces = int(line.split()[1])
          total_num_faces += num_faces
          domain_num_faces += num_faces

        elif line.startswith("vertex"):
          num_vertices = int(line.split()[1])
          total_num_vertices += num_vertices
          domain_num_vertices += num_vertices

        elif line.startswith("bounds"):
          bounds = line.split()[1:]
          model_bounds_list.append([float(i) for i in bounds])
          ply_bounds_str = ' '.join(bounds)

        elif line.startswith('#'):
          pass

        elif line.startswith('file'):
          pass

        else:
          print("identifier not found in line: " + line)
          sys.exit(1)

        # if (not args.abspath) and ("file" in line):
        #   filename_only = os.path.basename(line.split()[1])
        #   ftemp.write("file " + filename_only + "\n")
        # else:
        #   ftemp.write(line)

      with open(outfile, "a") as fout:
        fout.write("ModelBegin\n")

        if args.abspath:
          fout.write("file " + ply + "\n")
        else:
          fout.write("file " + os.path.basename(ply) + "\n")
      
        if ply_bounds_str:
          fout.write("ModelBounds " + ply_bounds_str + "\n")
      
        fout.write("material matte " + rgb_str + "\n")
        fout.write("# number of vertices: " + str(num_vertices)  + "\n")
        fout.write("# number of faces: " + str(num_faces)  + "\n")
        fout.write("ModelEnd\n")

    with open(outfile, "a") as fout:
      if not domain_bounds_str: 
        if len(model_bounds_list) == num_ply_models:
          domain_bounds_list = mergeBounds(model_bounds_list) 
          f.write("DomainBounds " + ' '.join(domain_bounds_list) + "\n")

      fout.write("# [domain " + str(domain_id) + "] number of vertices: " + str(domain_num_vertices)  + "\n")
      fout.write("# [domain " + str(domain_id) + "] number of faces: " + str(domain_num_faces)  + "\n")
      fout.write("DomainEnd\n")

  with open(outfile, "a") as fout:
    fout.write("######################\n")
    fout.write("# total vertices " + str(total_num_vertices) + "\n")
    fout.write("# total faces " + str(total_num_faces) + "\n")

def generateScene(args, plyfiles, outfile, bounds_info_dict):
  for domain_id in sorted(plyfiles):
    
    domain_bounds_str = None # domain bounds
    if bounds_info_dict:

      if str(domain_id) not in bounds_info_dict:
        print("[error] unknown domain ID found:" + str(domain_id))
        sys.exit(1)

      coordinates_list = bounds_info_dict[str(domain_id)]
      domain_bounds_str = ' '.join(coordinates_list)

    with open(outfile, "a") as f:
      f.write("######################\n")
      f.write("# " + str(domain_id) + "\n")
      f.write("DomainBegin\n")

      if domain_bounds_str:
        f.write("DomainBounds " + domain_bounds_str + "\n")

    for ply, iso in plyfiles[domain_id]:
      # rgba = valueToColor(iso, args.contour_range[0], args.contour_range[1], args.colormap[0], args.logarithmic_colormap)
      rgba = valueToColor(iso, args.contour_range[0], args.contour_range[1], args.colormap[0])
      rgb_str = RgbaToRgbString(rgba)
      with open(outfile, "a") as f:
        f.write("ModelBegin\n")
      
        if args.abspath:
          f.write("file " + ply + "\n")
        else:
          f.write("file " + os.path.basename(ply) + "\n")
    
        f.write("material matte " + rgb_str + "\n")
        f.write("ModelEnd\n")

    with open(outfile, "a") as f:
      f.write("DomainEnd\n")
  
def parseBoundsLine(line):
  tokens = line.split() 

  if len(tokens) != 8:
    print("unknown format (not 8 elements in a line): " + line)
    sys.exit(1)

  cluster_id = tokens[0]
  domain_id = tokens[1]

  coordinates=[]

  for i, v in enumerate(tokens):
    if i > 1:
      coordinates.append(v)

  return cluster_id, domain_id, coordinates

def parseBoundsFile(bounds_file):
  # key: domain id, value: list of coordinates (x0,y0,z0,x1,y1,z1)
  bounds={}
  prev_cluster_id=None
  with open(bounds_file, "r") as f:
    for line in f:
      if not line.startswith('#'):
        cluster_id, domain_id, coordinates_list = parseBoundsLine(line)

        if prev_cluster_id == None:
          prev_cluster_id = cluster_id

        elif prev_cluster_id != cluster_id:
          print("unknown cluster ID found: ", line)
          sys.exit(1)

        if domain_id in bounds:
          print("not a unique domain ID: ", line)
          sys.exit(1)

        bounds[domain_id] = coordinates_list

  return bounds

def main():
  parser = argparse.ArgumentParser(description='SpRay scene generator')
  parser.add_argument('--bounds-file', nargs=1, required=False, help='.bounds file describing bounds of domains')
  parser.add_argument('--indir', nargs='+', required=True, help='a list of directories containing ply files')
  parser.add_argument('--loader', nargs=1, required=False, help='executable for ply loader')
  parser.add_argument('--out', nargs=1, default=['scene.domain'], help='output file')
  parser.add_argument('--abspath', action='store_true', help='use absoulte path for ply files')
  # parser.add_argument('--logarithmic-colormap', action='store_true', help='Use logarithmic normalization for color mapping')
  parser.add_argument('--contour-range', nargs=2, required=False, default=[0.0, 2.0], help='minimum and maximum contour values')
  parser.add_argument('--colormap', nargs=1, required=False, default=['coolwarm'], help='colormap name. See the list of colormap names in the link: https://matplotlib.org/examples/color/colormaps_reference.html')
  
  # parser.add_argument('--scale', nargs=3, type=float, help='scaling factor (x,y,z)')
  # parser.add_argument('--gap', nargs=1, type=float, default=[1.0], help='gap factor (between domains)')
  args = parser.parse_args()
  
  # print("program dump_ply_info indir [outfile]")
  
  # if len(sys.argv) < 3:
  #   print("[error] no input directory found")
  #   sys.exit(1)
 
  bounds_file=None 
  if args.bounds_file:
    bounds_file = args.bounds_file[0]
    if not os.path.isfile(bounds_file):
      print("[error] bounds file not found: " + bounds_file)
      sys.exit(1)

  # program to dump info
  if args.loader:
    loader_executable = args.loader[0]
  
  outfile = args.out[0]
  if os.path.isfile(outfile):
    print("[error] output file already exists: " + outfile)
    sys.exit(1)
  # os.system("rm -f " + outfile)

  print("[python] listing ply files")
  
  ply_dirs = args.indir
  plyfiles={}
  plyfiles_abspath=[]
  
  for ply_dir in ply_dirs:
    filenames = os.listdir(ply_dir)
    dir_name = os.path.abspath(ply_dir)
  
    c=0  
    for filename in filenames:
      if filename.endswith(".ply"):
        plyfile_name = os.path.join(dir_name, filename)
        domain_id = getDomainId(filename)
        iso_value = getIsoValue(filename)

        if domain_id in plyfiles:
          plyfiles[domain_id].append((plyfile_name, iso_value))
        else:
          plyfiles[domain_id] = [(plyfile_name, iso_value)]

        print("[python] " + plyfile_name)
        c += 1
    
    if c == 0: 
      print("[warning] empty directory " + ply_dir)
      # sys.exit(1)

  # add lights
  addLights(outfile)

  bounds_info_dict=None
  if bounds_file:
    # key: domain id, value: list of coordinates (x0,y0,z0,x1,y1,z1)
    bounds_info_dict = parseBoundsFile(bounds_file)
  
  if args.loader:
    loadPlyFilesAndGenerateScene(args, plyfiles, loader_executable, outfile, bounds_info_dict)
  else:
    generateScene(args, plyfiles, outfile, bounds_info_dict)
      
  print("[info] generated a scene file in " + outfile)

if __name__ == '__main__':
  main()

