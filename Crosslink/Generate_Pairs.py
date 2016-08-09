#!/usr/bin/python

# Script to generate a pairs list based on bond information in a Gromacs topology
# The pairs list lists atoms that are 1-4 neighbors
import numpy as np
import os
import sys
import time
import pprint as pp
from scipy.sparse import csr_matrix
import genpairs

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

f = open('%s/../Structure-Files/Monomer_Tops/HII_mon.itp' %location, 'r')

a = []
for line in f:
    a.append(line)

# locate [ atoms ] section

atoms_index = 0  # find index where [ atoms ] section begins
while a[atoms_index].count('[ atoms ]') == 0:
    atoms_index += 1

# count number of atoms
atoms_count = atoms_index + 2
nr = 0  # number of lines in 'atoms' section
while a[atoms_count] != '\n':
    atoms_count += 1  # increments the while loop
    nr += 1  # counts number of atoms

# locate [ bonds ] section

bonds_index = 0  # find index where [ bonds ] section begins
while a[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

nb = 0  # number of lines in the 'bonds' section
bond_count = bonds_index + 2
while a[bond_count] != '\n':
    bond_count += 1  # increments while loop
    nb += 1  # counting number of lines in 'bonds' section

bonds_list = []
for i in range(bonds_index + 2, bond_count):
    bonds_list.append(a[i])


# [ dihedrals ] ; propers

dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
while a[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
    dihedrals_p_index += 1

ndp = 0  # number of lines in the 'dihedrals ; proper' section
dihedrals_p_count = dihedrals_p_index + 3  # keep track of index of a
dihedrals_prev = []
while a[dihedrals_p_count] != '\n':
    dihedrals_prev.append([int(a[dihedrals_p_count][0:6]), int(a[dihedrals_p_count][6:13]), int(a[dihedrals_p_count][13:20]), int(a[dihedrals_p_count][20:27])])
    dihedrals_p_count += 1
    ndp += 1

# [ dihedrals ] ; impropers

dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
while a[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
    dihedrals_imp_index += 1

ndimp = 0  # number of lines in the 'dihedrals ; impropers' section
dihedrals_imp_count = dihedrals_imp_index + 2
for i in range(dihedrals_imp_count, len(a)):  # This is the last section in the input .itp file
    dihedrals_imp_count += 1
    ndimp += 1

# [ pairs ]

pairs_index = 0  # find index where [ pairs ] section begins
while a[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

npair = 0  # number of lines in the 'pairs' section
pairs_count = pairs_index + 2  # keep track of index of a
pairs_list_prev = []
while a[pairs_count] != '\n':
    pairs_list_prev.append([int(a[pairs_count][0:6]), int(a[pairs_count][6:13])])  # record pairs that already exist
    pairs_count += 1
    npair += 1

# [ angles ]

angles_index = 0  # find index where [ angles ] section begins
while a[angles_index].count('[ angles ]') == 0:
    angles_index += 1

na = 0  # number of lines in the 'angles' section
angle_count = angles_index + 2  # keep track of index of a
angles_prev = []
while a[angle_count] != '\n':
    angles_prev.append([int(a[angle_count][0:6]), int(a[angle_count][6:13]), int(a[angle_count][13:20])])
    angle_count += 1
    na += 1
print na
start = time.time()

# atoms_of_interest = []
# for i in range(0, 137):
#     atoms_of_interest.append(i + 1)
atoms_of_interest = [1328, 3345, 3759, 5189, 6062, 6604, 7048, 8463, 8907, 11338, 12196, 13769, 14484, 15628, 17359, 18217, 19474, 20919, 21649, 21920, 22221, 22650, 22793, 23621, 24908, 25623, 25781, 27083, 28227, 29341, 29356, 29770, 29785, 30801, 32231, 32630, 33089, 35520, 36077, 36363, 36506, 37078, 38666, 42226, 43513, 44243, 46787, 46817, 47502, 47645, 49504, 49948, 50535, 50806, 51378, 51664, 53109, 53222, 53380, 53523, 54795, 56940, 57083, 57098, 59529, 60515, 61087, 61245, 62833, 64692, 68538, 1329, 3346, 3760, 5190, 6063, 6605, 7049, 8464, 8908, 11339, 12197, 13770, 14485, 15629, 17360, 18218, 19475, 20920, 21650, 21921, 22222, 22651, 22794, 23622, 24909, 25624, 25782, 27084, 28228, 29342, 29357, 29771, 29786, 30802, 32232, 32631, 33090, 35521, 36078, 36364, 36507, 37079, 38667, 42227, 43514, 44244, 46788, 46818, 47503, 47646, 49505, 49949, 50536, 50807, 51379, 51665, 53110, 53223, 53381, 53524, 54796, 56941, 57084, 57099, 59530, 60516, 61088, 61246, 62834, 64693, 68539, 35949, 3330, 3774, 6047, 40668, 56955, 24509, 41669, 42098, 64248, 29627, 14627, 14326, 14785, 455, 35919, 36378, 54968, 40224, 55653, 38952, 22635, 6190, 22763, 8335, 60530, 60372, 26210, 29070, 28641, 11609, 28912, 28942, 49233, 67665, 31772, 33774, 18646, 35219, 2615, 36521, 37063, 39524, 42241, 60116, 61960, 64534, 47660, 12068, 64549, 32088, 49105, 33789, 50791, 15914, 34775, 53094, 37522, 35663, 53538, 19361, 22349, 56255, 23049, 24080, 43099, 10036, 26654, 45959, 64519, 16486, 35950, 3331, 3775, 6048, 40669, 56956, 24510, 41670, 42099, 64249, 29628, 14628, 14327, 14786, 456, 35920, 36379, 54969, 40225, 55654, 38953, 22636, 6191, 22764, 8336, 60531, 60373, 26211, 29071, 28642, 11610, 28913, 28943, 49234, 67666, 31773, 33775, 18647, 35220, 2616, 36522, 37064, 39525, 42242, 60117, 61961, 64535, 47661, 12069, 64550, 32089, 49106, 33790, 50792, 15915, 34776, 53095, 37523, 35664, 53539, 19362, 22350, 56256, 23050, 24081, 43100, 10037, 26655, 45960, 64520, 16487, 3172, 30214, 65836, 3173, 30215, 65837, 36028, 36029, 36030, 1398, 1399, 3399, 3400, 3401, 3425, 3426, 3853, 3854, 3855, 3829, 3830, 6116, 6117, 6118, 5259, 5260, 40747, 40748, 40749, 6142, 6143, 57024, 57025, 57026, 6664, 6665, 24588, 24589, 24590, 7118, 7119, 41748, 41749, 41750, 8523, 8524, 42177, 42178, 42179, 8977, 8978, 64317, 64318, 64319, 11408, 11409, 29686, 29687, 29688, 12266, 12267, 14696, 14697, 14698, 13839, 13840, 14385, 14386, 14387, 14554, 14555, 14864, 14865, 14866, 15698, 15699, 514, 515, 516, 17439, 17440, 35978, 35979, 35980, 18297, 18298, 36457, 36458, 36459, 19534, 19535, 55047, 55048, 55049, 20989, 20990, 40293, 40294, 40295, 21729, 21730, 55712, 55713, 55714, 21990, 21991, 39031, 39032, 39033, 22301, 22302, 22704, 22705, 22706, 22730, 22731, 6259, 6260, 6261, 22873, 22874, 22822, 22823, 22824, 23681, 23682, 8404, 8405, 8406, 24968, 24969, 60599, 60600, 60601, 25683, 25684, 60431, 60432, 60433, 25851, 25852, 26279, 26280, 26281, 27163, 27164, 29139, 29140, 29141, 28307, 28308, 28710, 28711, 28712, 29401, 29402, 11668, 11669, 11670, 29426, 29427, 28971, 28972, 28973, 29830, 29831, 29021, 29022, 29023, 29855, 29856, 49302, 49303, 49304, 30881, 30882, 67724, 67725, 67726, 32311, 32312, 31831, 31832, 31833, 32690, 32691, 33833, 33834, 33835, 33169, 33170, 18725, 18726, 18727, 35600, 35601, 35288, 35289, 35290, 36147, 36148, 2684, 2685, 2686, 36433, 36434, 36600, 36601, 36602, 36576, 36577, 37122, 37123, 37124, 37148, 37149, 39603, 39604, 39605, 38746, 38747, 42320, 42321, 42322, 42296, 42297, 60195, 60196, 60197, 43583, 43584, 62029, 62030, 62031, 44323, 44324, 64603, 64604, 64605, 46847, 46848, 47729, 47730, 47731, 46897, 46898, 12147, 12148, 12149, 47562, 47563, 64628, 64629, 64630, 47705, 47706, 32167, 32168, 32169, 49564, 49565, 49184, 49185, 49186, 50018, 50019, 33858, 33859, 33860, 50615, 50616, 50850, 50851, 50852, 50876, 50877, 15983, 15984, 15985, 51448, 51449, 34834, 34835, 34836, 51734, 51735, 53163, 53164, 53165, 53189, 53190, 37601, 37602, 37603, 53282, 53283, 35742, 35743, 35744, 53450, 53451, 53617, 53618, 53619, 53593, 53594, 19440, 19441, 19442, 54855, 54856, 22418, 22419, 22420, 57000, 57001, 56334, 56335, 56336, 57143, 57144, 23108, 23109, 23110, 57168, 57169, 24159, 24160, 24161, 59599, 59600, 43178, 43179, 43180, 60575, 60576, 10095, 10096, 10097, 61147, 61148, 26733, 26734, 26735, 61315, 61316, 46038, 46039, 46040, 62913, 62914, 64578, 64579, 64580, 64772, 64773, 16555, 16556, 16557, 68608, 68609, 3231, 3232, 3233, 30283, 30284, 30285, 65915, 65916, 65917, 36036, 3430, 3861, 6147, 40755, 57055, 24596, 41756, 42185, 64348, 29740, 14727, 14439, 14872, 568, 36032, 36465, 55055, 40324, 55766, 39039, 22735, 6290, 22876, 8435, 60630, 60485, 26310, 29170, 28741, 11722, 29025, 29029, 49333, 67778, 31885, 33887, 18733, 35319, 2715, 36608, 37176, 39611, 42328, 60203, 62060, 64634, 47760, 12155, 64636, 32175, 49192, 33889, 50904, 16014, 34888, 53194, 37609, 35750, 53625, 19448, 22449, 56342, 23162, 24167, 43186, 10149, 26741, 46046, 64632, 16586, 3284, 3285, 30313, 30314, 65922, 65923]


# Make a condensed list of bonds including only the bonds which we are interested in as well as a full list of bonds:

def read_lists(a, atoms_of_interest, start, end):
    list = []
    condensed = []
    for i in range(start, end):  # extract atom numbers in [ bonds ] section in order of appearance
        list.append([int(a[i][0:6]), int(a[i][6:13])])
        if int(a[i][0:6]) in atoms_of_interest:
            condensed.append([int(a[i][0:6]), int(a[i][6:13])])
        if int(a[i][6:13]) in atoms_of_interest:
            condensed.append([int(a[i][0:6]), int(a[i][6:13])])

    uniq_condensed = []
    for i in range(0, len(condensed)):
        a = condensed[i]
        b = [condensed[i][1], condensed[i][0]]
        if a not in uniq_condensed:
            if b not in uniq_condensed:
                uniq_condensed.append(a)

    return list, uniq_condensed

list, uniq_condensed = read_lists(a, atoms_of_interest, bonds_index + 2, bond_count)
print 'length of list: %s' %len(uniq_condensed)
start = time.time()


def gen_pairs_list(list, condensed):
    pairs = []
    angles = []
    dihedrals = []
    for i in range(0, len(condensed)):  # look at all bonds
        neighbor3 = []  # redefine lists for each loop since each loops looks at a different bond
        neighbor4 = []
        for j in range(0, len(list)):  # compare to all other bonds
            if list[j][0] == condensed[i][1]:  # look for 1-3 pairs by comparing first entry of each row of [ pairs ] section
                neighbor3.append(list[j][1])
                angles.append([condensed[i][0], list[j][0], list[j][1]])
            elif list[j][1] == condensed[i][1]:  # look for 1-3 pairs again
                if list[j][0] != condensed[i][0]:
                    neighbor3.append(list[j][0])
                    angles.append([condensed[i][0], list[j][1], list[j][0]])
        for h in neighbor3:
            for l in range(0, len(list)):
                if list[l][0] == h:
                    if list[l][1] != condensed[i][1]:
                        dihedrals.append([condensed[i][0], condensed[i][1], h, list[l][1]])
                        neighbor4.append(list[l][1])
                elif list[l][1] == h:
                    if list[l][0] != condensed[i][1]:
                        dihedrals.append([condensed[i][0], condensed[i][1], h, list[l][0]])
                        neighbor4.append(list[l][0])
        count = 0
        for m in neighbor4:
            pairs.append([condensed[i][0], neighbor4[count]])
            count += 1

        neighbor3 = []
        neighbor4 = []
        for j in range(0, len(list)):
            if list[j][1] == condensed[i][0]:
                angles.append([condensed[i][1], list[j][1], list[j][0]])
                neighbor3.append(list[j][0])
            elif list[j][0] == condensed[i][0]:
                if list[j][1] != condensed[i][1]:
                    angles.append([condensed[i][1], list[j][0], list[j][1]])
                    neighbor3.append(list[j][1])
        for h in neighbor3:
            for l in range(0, len(list)):
                if list[l][0] == h:
                    if list[l][1] != condensed[i][0]:
                        dihedrals.append([condensed[i][1], condensed[i][0], h, list[l][1]])
                        neighbor4.append(list[l][1])
                elif list[l][1] == h:
                    if list[l][0] != condensed[i][0]:
                        dihedrals.append([condensed[i][1], condensed[i][0], h, list[l][0]])
                        neighbor4.append(list[l][0])

        count = 0
        for m in neighbor4:
            pairs.append([condensed[i][1], neighbor4[count]])
            count += 1
    return pairs, dihedrals, angles


end1 = time.time()
print end1 - start

def uniq_pair(pairs, pairs_list_prev):
    for i in range(0, len(pairs)):
        a = pairs[i]  # the order it is read from the pairs list
        b = [pairs[i][1], pairs[i][0]]
        if a not in pairs_list_prev:
            if b not in pairs_list_prev:
                pairs_list_prev.append(a)
    return pairs_list_prev

def uniq_dihedral(dihedrals, dihedrals_prev):
    for i in range(0, len(dihedrals)):
        a = dihedrals[i]
        b = [dihedrals[i][3], dihedrals[i][2], dihedrals[i][1], dihedrals[i][0]]
        if a not in dihedrals_prev:
            if b not in dihedrals_prev:
                dihedrals_prev.append(dihedrals[i])
    return dihedrals_prev

def uniq_angles(angles, angles_prev):
    for i in range(0, len(angles)):
        a = angles[i]
        b = [angles[i][2], angles[i][1], angles[i][0]]
        if a not in angles_prev:
            if b not in angles_prev:
                angles_prev.append(angles[i])
    return angles_prev

pairs, dihedrals, angles = gen_pairs_list(list, uniq_condensed)

unique_pairs = uniq_pair(pairs, pairs_list_prev)
unique_dihedrals = uniq_dihedral(dihedrals, dihedrals_prev)
unique_angles = uniq_angles(angles, angles_prev)
print len(unique_pairs)
print len(unique_dihedrals)
print len(unique_angles)
# for i in range(0, len(unique_pairs)):
#     print '{:4d}{:4d}'.format(unique_pairs[i][0], unique_pairs[i][1])
#
# print len(unique_pairs)
#
# unique_pairs = uniq_pair(pairs1, pairs2, pairs_list_prev)

