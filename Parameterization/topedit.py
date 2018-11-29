f = open("MOL_GMX.gro", "r")
a = []
no_lines = 0
for line in f:
    a.append(line)
    no_lines = no_lines + 1

no_atoms = int(a[1]) + 1
a[1] = " " + str(no_atoms) + "\n"

add = "    2  NA    NA  " + str(no_atoms) + "   0.000   0.000   0.000\n"
a.insert(no_lines - 1, add)
f.close()



f = open("MOL_GMX.gro", "w")
for line in a:
    f.write(line)
f.close()

f = open("MOL_GMX.top", "r")
top = []
itp = []
itp_temp = []
count = 0
start = 0
end = 0
no_lines_top = 0
for line in f:
    itp_temp.append(line)
    no_lines_top = no_lines_top + 1

i = 0
while i < no_lines_top:
    if itp_temp[i].startswith("[ moleculetype ]"):
        start = i
    if itp_temp[i].startswith("[ system ]"):
        end = i - 1
        break
    i = i + 1

i = start
while i < end:
    itp.append(itp_temp[i])
    i = i + 1

f.close()

f = open("MOL_GMX.itp", "w")
for line in itp:
    f.write(line)
f.close()

f = open("MOL_GMX.top", "w")
f.write(";Forcefield\n")
f.write('#include "amber99.ff/forcefield.itp"\n')
f.write('#include "/home/NormaKLangdon/2017/gaff/gaff.itp"\n')
f.write("\n")
f.write(";Monomer Topology\n")
f.write('#include "/home/NormaKLangdon/2017/MOL_GMX.itp"\n')
f.write("\n")
f.write(";Ion Topology\n")
f.write('#include "amber99.ff/ions.itp"\n')
f.write("\n")
f.write("[ system ]\n")
f.write("Vacuum Simulation\n")
f.write("\n")
f.write("[ molecules ]\n")
f.write("; Compound        nmols\n")
f.write(" MOL              1\n")
f.write(" NA               1\n")

f.close()