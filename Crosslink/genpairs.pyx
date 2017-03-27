#!/usr/bin/python

# Script to generate a pairs list based on bond information in a Gromacs topology
# The pairs list lists atoms that are 1-4 neighbors
import cython

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

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def gen_pairs_list(list1, list2, xlinked_list1, xlinked_list2, no_bonds):
# 	pairs1 = []
# 	pairs2 = []
# 	dihedrals = []
# 	dc = 0  # dihedrals count
# 	for i in range(0, len(xlinked_list1)):  # look at all bonds
# 		neighbor3 = []  # redefine lists for each loop since each loops looks at a different bond
# 		neighbor4 = []
# 		for j in range(0, no_bonds):  # compare to all other bonds
# 			if list1[j] == xlinked_list2[i]:  # look for 1-3 pairs by comparing first entry of each row of [ pairs ] section
# 				neighbor3.append(list2[j])
# 			elif list2[j] == xlinked_list2[i]:  # look for 1-3 pairs again
# 				if list1[j] != xlinked_list1[i]:
# 					neighbor3.append(list1[j])
# 		for h in neighbor3:
# 			for l in range(0, len(list1)):
# 				if list1[l] == h:
# 					if list2[l] != xlinked_list2[i]:
# 						dihedrals.append([])
# 						dihedrals[dc].append(xlinked_list1[i])
# 						dihedrals[dc].append(xlinked_list2[i])
# 						dihedrals[dc].append(h)
# 						dihedrals[dc].append(list2[l])
# 						dc += 1
# 						neighbor4.append(list2[l])
# 				elif list2[l] == h:
# 					if list1[l] != xlinked_list2[i]:
# 						dihedrals.append([])
# 						dihedrals[dc].append(xlinked_list1[i])
# 						dihedrals[dc].append(xlinked_list2[i])
# 						dihedrals[dc].append(h)
# 						dihedrals[dc].append(list1[l])
# 						dc += 1
# 						neighbor4.append(list1[l])
# 		count = 0
# 		for m in neighbor4:
# 			pairs1.append(xlinked_list1[i])
# 			pairs2.append(neighbor4[count])
# 			#pairs.append('{:6d}{:7d}{:7d}'.format(list1[i], neighbor4[count], 1))
# 			count += 1
#
# 		neighbor3 = []
# 		neighbor4 = []
# 		for j in range(0, no_bonds):
# 			if list2[j] == xlinked_list1[i]:
# 				neighbor3.append(list1[j])
# 			elif list1[j] == xlinked_list1[i]:
# 				if list2[j] != xlinked_list2[i]:
# 					neighbor3.append(list2[j])
# 		for h in neighbor3:
# 			for l in range(0, len(list1)):
# 				if list1[l] == h:
# 					if list2[l] != xlinked_list1[i]:
# 						dihedrals.append([])
# 						dihedrals[dc].append(xlinked_list2[i])
# 						dihedrals[dc].append(xlinked_list1[i])
# 						dihedrals[dc].append(h)
# 						dihedrals[dc].append(list2[l])
# 						dc += 1
# 						neighbor4.append(list2[l])
# 				elif list2[l] == h:
# 					if list1[l] != xlinked_list1[i]:
# 						dihedrals.append([])
# 						dihedrals[dc].append(xlinked_list2[i])
# 						dihedrals[dc].append(xlinked_list1[i])
# 						dihedrals[dc].append(h)
# 						dihedrals[dc].append(list1[l])
# 						dc += 1
# 						neighbor4.append(list1[l])
#
# 		count = 0
# 		for m in neighbor4:
# 			pairs1.append(xlinked_list2[i])
# 			pairs2.append(neighbor4[count])
# 			#pairs.append('{:6d}{:7d}{:7d}'.format(list2[i], neighbor4[count], 1))
# 			count += 1
# 	return pairs1, pairs2, dihedrals

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

#eliminate duplicates


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

def calc_dist(C1x, C2x, C1y, C2y, C1z, C2z, exclude, dist):
	for i in range(0, len(C1x)):
		for k in range(0, len(C2x)):
			if exclude[k, i] == 1:  # make sure that C1 and C2 that are a part of the same monomer do not factor into this calculation
				dist[k, i] = 1000  # artificially high number to keep it from interfering
			else:
				dist[k, i] = ((C1x[i] - C2x[k])**2 + (C1y[i] - C2y[k])**2 + (C1z[i] - C2z[k])**2)**(0.5)
	return dist

def improper_dihedrals(b, start_imp, imp_of_interest, dihedrals_imp_count):
	dihedrals_imp = []
	for i in range(start_imp, len(b)):  # This is the last section in the input .itp file
		a_imp = [int(b[i][0:6]), int(b[i][6:13]), int(b[i][13:20]), int(b[i][20:27])]
		for k in range(0, len(imp_of_interest)):
			if set(a_imp) != set(imp_of_interest[k]):
				if a_imp not in dihedrals_imp:
					dihedrals_imp.append(a_imp)
		dihedrals_imp_count += 1
	return dihedrals_imp