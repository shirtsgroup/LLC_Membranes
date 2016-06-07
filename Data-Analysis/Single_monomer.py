import matplotlib.pyplot as plt

f = open("P_wiggle_1ns.xvg")
a = []
for line in f:
    a.append(line)

x_values = []  # time (ps)
y_values = []  # energy (kJ/mol)

for i in range(23, len(a)):
    x_values.append(float(a[i][0:12]))
    y_values.append(float(a[i][14:28]))

print(len(x_values))
print min(y_values)

plt.plot(x_values, y_values)
# plt.ylabel('Energy (kJ/mol)')
plt.ylabel('Pressure (bar)')
plt.xlabel('Time (ps)')
plt.title('GROMACS Energies')
plt.show()
