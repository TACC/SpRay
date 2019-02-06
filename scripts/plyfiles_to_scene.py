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

def RgbaToRgbString(rgba):
  rgb = rgba[:-1] 
  rgb_str = ' '.join(str(e) for e in rgb)
  return rgb_str

def valueToColor(value, value_min, value_max):
  norm = mpl.colors.Normalize(vmin=value_min, vmax=value_max, clip=True)
  mapper = cm.ScalarMappable(norm=norm, cmap=cm.coolwarm)
  rgba = mapper.to_rgba(value)
  return rgba

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

def loadPlyFilesAndGenerateScene(args, plyfiles, loader_executable, outfile):
  total_num_faces = 0 
  total_num_vertices = 0 

  for domain_id in sorted(plyfiles):

    with open(outfile, "a") as fout:
      fout.write("######################\n")
      fout.write("# " + str(domain_id) + "\n")
      fout.write("DomainBegin\n")

    domain_num_faces = 0 
    domain_num_vertices = 0

    for ply, iso in plyfiles[domain_id]:
      rgba = valueToColor(iso, args.contour_range[0], args.contour_range[1])
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

      for line in lines:
        if "face" in line:
          num_faces = int(line.split()[1])
          total_num_faces += num_faces
          domain_num_faces += num_faces
        if "vertex" in line:
          num_vertices = int(line.split()[1])
          total_num_vertices += num_vertices
          domain_num_vertices += num_vertices
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
      
        fout.write("material matte " + rgb_str + "\n")
        fout.write("# number of vertices: " + str(num_vertices)  + "\n")
        fout.write("# number of faces: " + str(num_faces)  + "\n")
        fout.write("ModelEnd\n")

    with open(outfile, "a") as fout:
      fout.write("# [domain " + str(domain_id) + "] number of vertices: " + str(domain_num_vertices)  + "\n")
      fout.write("# [domain " + str(domain_id) + "] number of faces: " + str(domain_num_faces)  + "\n")
      fout.write("DomainEnd\n")

  with open(outfile, "a") as fout:
    fout.write("######################\n")
    fout.write("# total vertices " + str(total_num_vertices) + "\n")
    fout.write("# total faces " + str(total_num_faces) + "\n")

def generateScene(args, plyfiles, outfile):
  for domain_id in sorted(plyfiles):
    with open(outfile, "a") as f:
      f.write("######################\n")
      f.write("# " + str(domain_id) + "\n")
      f.write("DomainBegin\n")
    for ply, iso in plyfiles[domain_id]:
      rgba = valueToColor(iso, args.contour_range[0], args.contour_range[1])
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
  
def main():
  parser = argparse.ArgumentParser(description='Domains maker')
  parser.add_argument('--indir', nargs='+', required=True, help='a list of directories containing ply files')
  parser.add_argument('--loader', nargs=1, required=False, help='executable for ply loader')
  parser.add_argument('--out', nargs=1, default=['scene.domain'], help='output file')
  parser.add_argument('--abspath', action='store_true', help='use absoulte path for ply files')
  parser.add_argument('--contour-range', nargs=2, required=False, default=[0.0, 2.0], help='minimum and maximum contour values')
  
  
  # parser.add_argument('--scale', nargs=3, type=float, help='scaling factor (x,y,z)')
  # parser.add_argument('--gap', nargs=1, type=float, default=[1.0], help='gap factor (between domains)')
  args = parser.parse_args()
  
  # print("program dump_ply_info indir [outfile]")
  
  # if len(sys.argv) < 3:
  #   print("[error] no input directory found")
  #   sys.exit(1)
  
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
  
  if args.loader:
    loadPlyFilesAndGenerateScene(args, plyfiles, loader_executable, outfile)
  else:
    generateScene(args, plyfiles, outfile)
      
  print("[info] generated a scene file in " + outfile)

if __name__ == '__main__':
  main()

