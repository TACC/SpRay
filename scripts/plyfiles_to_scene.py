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

parser = argparse.ArgumentParser(description='Domains maker')
parser.add_argument('--indir', nargs='+', required=True, help='a list of directories containing ply files')
parser.add_argument('--loader', nargs=1, required=True, help='executable for ply loader')
parser.add_argument('--out', nargs=1, default=['scene.domain'], help='output file')


# parser.add_argument('--scale', nargs=3, type=float, help='scaling factor (x,y,z)')
# parser.add_argument('--gap', nargs=1, type=float, default=[1.0], help='gap factor (between domains)')
args = parser.parse_args()

# print("program dump_ply_info indir [outfile]")

# if len(sys.argv) < 3:
#   print("[error] no input directory found")
#   sys.exit(1)

# program to dump info
exe_name = args.loader[0]


print("[python] listing ply files")

ply_dirs = args.indir
plyfiles=[]

for ply_dir in ply_dirs:
  filenames = os.listdir(ply_dir)
  dir_name = os.path.abspath(ply_dir)

  c=0  
  for filename in filenames:
    if filename.endswith(".ply"):
      plyfile_name = os.path.join(dir_name, filename)
      plyfiles.append(plyfile_name)
      print("[python] " + plyfile_name)
      c += 1
  
  if c == 0: 
    print("empty directory " + ply_dir)
    sys.exit(1)


outfile = args.out[0]
os.system("rm -f " + outfile)

n = 0 
num_faces = 0 
num_vertices = 0 
for ply in plyfiles:
  cmd = exe_name + " " + ply + " out.tmp"
  print("[python] " + cmd)
  os.system(cmd)

  lines=[] 
  with open("./out.tmp", "r") as f:
    for line in f:
      lines.append(line)

  os.system("rm -f ./out.tmp")

  with open(outfile, "a") as f:
    f.write("######################\n")
    f.write("# " + str(n) + "\n")
    f.write("domain\n")
    for line in lines:
      f.write(line)
      if "face" in line:
        num_faces += int(line.split()[1])
      if "vertex" in line:
        num_vertices += int(line.split()[1])
    f.write("mtl diffuse 1 1 1\n")
  
  n += 1

with open(outfile, "a") as f:
  f.write("######################\n")
  f.write("# total vertices " + str(num_vertices) + "\n")
  f.write("# total faces " + str(num_faces) + "\n")
