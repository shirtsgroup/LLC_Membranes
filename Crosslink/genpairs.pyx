#!/usr/bin/python

# Script to generate a pairs list based on bond information in a Gromacs topology
# The pairs list lists atoms that are 1-4 neighbors
import cython

def read_lists(a, start, end):
	list1 = []  # list to hold atom numbers of 1st entry in each row of [ bonds ] section
	list2 = []  # list to hold atom numbers of 2nd entry in each row of [ bonds ] section
	for i in range(start, end):  # extract atom numbers in [ bonds ] section in order of appearance
		list1.append(int(a[i][0:6]))  # converts them to integers to save time later
		list2.append(int(a[i][6:13]))

	return list1, list2

@cython.boundscheck(False)
@cython.wraparound(False)
def gen_pairs_list(list1, list2, xlinked_list1, xlinked_list2, no_bonds):
	pairs1 = []
	pairs2 = []
	dihedrals = []
	dc = 0  # dihedrals count
	for i in range(0, len(xlinked_list1)):  # look at all bonds
		neighbor3 = []  # redefine lists for each loop since each loops looks at a different bond
		neighbor4 = []
		for j in range(0, no_bonds):  # compare to all other bonds
			if list1[j] == xlinked_list2[i]:  # look for 1-3 pairs by comparing first entry of each row of [ pairs ] section
				neighbor3.append(list2[j])
			elif list2[j] == xlinked_list2[i]:  # look for 1-3 pairs again
				if list1[j] != xlinked_list1[i]:
					neighbor3.append(list1[j])
		for h in neighbor3:
			for l in range(0, len(list1)):
				if list1[l] == h:
					if list2[l] != xlinked_list2[i]:
						dihedrals.append([])
						dihedrals[dc].append(xlinked_list1[i])
						dihedrals[dc].append(xlinked_list2[i])
						dihedrals[dc].append(h)
						dihedrals[dc].append(list2[l])
						dc += 1
						neighbor4.append(list2[l])
				elif list2[l] == h:
					if list1[l] != xlinked_list2[i]:
						dihedrals.append([])
						dihedrals[dc].append(xlinked_list1[i])
						dihedrals[dc].append(xlinked_list2[i])
						dihedrals[dc].append(h)
						dihedrals[dc].append(list1[l])
						dc += 1
						neighbor4.append(list1[l])
		count = 0
		for m in neighbor4:
			pairs1.append(xlinked_list1[i])
			pairs2.append(neighbor4[count])
			#pairs.append('{:6d}{:7d}{:7d}'.format(list1[i], neighbor4[count], 1))
			count += 1

		neighbor3 = []
		neighbor4 = []
		for j in range(0, no_bonds):
			if list2[j] == xlinked_list1[i]:
				neighbor3.append(list1[j])
			elif list1[j] == xlinked_list1[i]:
				if list2[j] != xlinked_list2[i]:
					neighbor3.append(list2[j])
		for h in neighbor3:
			for l in range(0, len(list1)):
				if list1[l] == h:
					if list2[l] != xlinked_list1[i]:
						dihedrals.append([])
						dihedrals[dc].append(xlinked_list2[i])
						dihedrals[dc].append(xlinked_list1[i])
						dihedrals[dc].append(h)
						dihedrals[dc].append(list2[l])
						dc += 1
						neighbor4.append(list2[l])
				elif list2[l] == h:
					if list1[l] != xlinked_list1[i]:
						dihedrals.append([])
						dihedrals[dc].append(xlinked_list2[i])
						dihedrals[dc].append(xlinked_list1[i])
						dihedrals[dc].append(h)
						dihedrals[dc].append(list1[l])
						dc += 1
						neighbor4.append(list1[l])

		count = 0
		for m in neighbor4:
			pairs1.append(xlinked_list2[i])
			pairs2.append(neighbor4[count])
			#pairs.append('{:6d}{:7d}{:7d}'.format(list2[i], neighbor4[count], 1))
			count += 1
	return pairs1, pairs2, dihedrals

#eliminate duplicates


def uniq_pair(pairs1, pairs2, pairs_list_prev):
	for i in range(0, len(pairs1)):
		a = [pairs1[i], pairs2[i]]  # the order it is read from the pairs list
		b = [pairs2[i], pairs1[i]]
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
