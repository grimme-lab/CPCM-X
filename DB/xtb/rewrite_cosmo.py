#!/bin/python3

# Read .cosmo file and write a fortran file with the same data

import sys

# Break each line after every 120 characters
def break_lines(string):
    if (len(string) < 132):
        return string
    new_string = ""
    for i in range(0, len(string), 120):
        new_string += string[i:i+120] + "&\n&"
    return new_string[:-3]
# Read the file
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

solvent = sys.argv[1].split('.')[0]

# Read area
area=0
for line in lines:
    if "area" in line:
        area = float(line.split('=')[1])
        break
# Read volume (if it exists)
volume=0
for line in lines:
    if "volume" in line:
        volume = float(line.split('=')[1])
        break
# Read COSMO energy
energy=0
for line in lines:
    if "Total energy [a.u.]" in line:
        energy = float(line.split('=')[1])
        break

print(area, volume, energy)

# Read atom information
aid = []
axyz = []
aelement = []
i=0
while True:
    if "#atom" in lines[i]:
        break
    i+=1
i+=1
tmp_element=""
while 'coord_car' not in lines[i]:
    aid.append(int(lines[i].split()[0]))
    axyz.append([float(lines[i].split()[1]), float(lines[i].split()[2]), float(lines[i].split()[3])])
    tmp_element=lines[i].split()[4].lower()
    if len(tmp_element) == 1:
        tmp_element=tmp_element+" "
    aelement.append(tmp_element)
    i+=1

# Read segment information
su = []
sid = []
sarea = []
sxyz = []
spot = []
while True:
    if "segment_information" in lines[i]:
        break
    i+=1
i+=11
while i < len(lines):
    sid.append(int(lines[i].split()[1]))
    sxyz.append([float(lines[i].split()[2]), float(lines[i].split()[3]), float(lines[i].split()[4])])
    sarea.append(float(lines[i].split()[6]))
    su.append(float(lines[i].split()[7]))
    spot.append(float(lines[i].split()[8]))
    i+=1

# Write the file
with open(solvent+".fh", 'w') as f:
    f.write("!>COSMO file for "+solvent+" from xtb\n")
    f.write(break_lines("real(wp), parameter :: area_"+solvent+"_xtb = "+str(area)+", volume_"+solvent+"_xtb = "+str(volume)+", energy_"+solvent+"_xtb = "+str(energy)+"\n"))
    f.write(break_lines("integer, parameter, dimension("+str(len(aid))+") :: aid_"+solvent+"_xtb = (/ "+str(aid)[1:-1]+" /)\n"))
    f.write(break_lines("real(wp), parameter, dimension("+str(len(aid))+",3) :: axyz_"+solvent+"_xtb = reshape(btoa*(/ "+str(axyz)[1:-1]+" /), shape(axyz_"+solvent+"_xtb),order=(/2,1/))\n"))
    f.write(break_lines("character(len=2), parameter, dimension("+str(len(aid))+") :: aelement_"+solvent+"_xtb = (/ "+str(aelement)[1:-1]+" /)\n"))
    f.write(break_lines("integer, parameter, dimension("+str(len(sid))+") :: sid_"+solvent+"_xtb = (/ "+str(sid)[1:-1]+" /)\n"))
    f.write(break_lines("real(wp), parameter, dimension("+str(len(sid))+",3) :: sxyz_"+solvent+"_xtb = reshape(btoa*(/ "+str(sxyz)[1:-1]+" /), shape(sxyz_"+solvent+"_xtb),order=(/2,1/))\n"))
    f.write(break_lines("real(wp), parameter, dimension("+str(len(sid))+") :: sarea_"+solvent+"_xtb = (/ "+str(sarea)[1:-1]+" /)\n"))
    f.write(break_lines("real(wp), parameter, dimension("+str(len(sid))+") :: su_"+solvent+"_xtb = (/ "+str(su)[1:-1]+" /)\n"))
    f.write(break_lines("real(wp), parameter, dimension("+str(len(sid))+") :: spot_"+solvent+"_xtb = (/ "+str(spot)[1:-1]+" /)\n"))


