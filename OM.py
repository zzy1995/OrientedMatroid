import itertools 
import numpy as np
import networkx as nx
from math import factorial
from networkx.algorithms.distance_measures import diameter
import sys
import matplotlib.pyplot as plt
from scipy.linalg import null_space

EPS = np.finfo(np.float).eps
Todd_hyperplane = np.array([[6, 3, 0, -1, -1],
							[3, 6, -1, 0, -1],
							[35, 45, -6, -3, -8],
							[45, 35, -3, -6, -8],
							[1, 0, 0, 0, 0],
							[0, 1, 0, 0, 0],
							[0, 0, 1, 0, 0],
							[0, 0, 0, 1, 0]])

# s is the set of elements, r is the number of elements in a subset
def findsubsets(n, r): 
    if n == r:
    	return [list(range(n))]
    if r == 1:
    	return [[i] for i in range(n)]
    first = findsubsets(n - 1, r)
    second = findsubsets(n - 1, r - 1)
    for items in second:
    	items.append(n - 1)
    return first + second
      

def plusminus_to_num(string):
	if string == "+":
		return 1
	if string == "-":
		return -1
	if string == "0":
		return 0
	raise Exception("string is not +-0")

# construct a chirotope map from r-subset of [n], string is +-0 type
def chirotope(n, r, string):
	revlex = findsubsets(n, r)
	chirotope_map = {}
	for i in range(len(revlex)):
		chirotope_map[frozenset(revlex[i])] = plusminus_to_num(string[i])
	return chirotope_map

def circuits(n, r, chirotope_map):
	circuit = []
	circuit_nonzero = [i for i in itertools.combinations(range(n), r + 1)]
	for items in circuit_nonzero:
		signature = []
		for i in range(n):
			if i not in items:
				signature.append(0)
			elif i == min(items):
				signature.append(1)
			else:
				exchange = sum([i < x for x in items])
				new_set = list(items)
				del new_set[0]
				newlist = [x for x in new_set if x != i]
				newlist.append(items[0])
				signature.append(chirotope_map[frozenset(newlist)] 
					* chirotope_map[frozenset(new_set)] * ((-1) ** exchange))
		signature = np.array(signature)
		circuit.append(signature)
		circuit.append(-signature)
	return circuit


# return a list of np arrays of cocircuits signature
def cocircuits(n, r, chirotope_map):
	cocircuit = []
	cocircuit_zero = [i for i in itertools.combinations(range(n), r - 1)]
	for items in cocircuit_zero:
		signature = []
		for i in range(n):
			if i in items:
				signature.append(0)
			else:
				exchange = sum([i < x for x in items])
				new_set = list(items)
				new_set.append(i)
				signature.append(chirotope_map[frozenset(new_set)] * ((-1) ** exchange))
		signature = np.array(signature)
		if not True in [np.array_equal(signature, item) for item in cocircuit]:
			cocircuit.append(signature)
			cocircuit.append(-signature)
	return cocircuit

def make_graph(n, r, cocircuit):
	# This code does not work in non-uniform case
	#if factorial(n) // factorial(r - 1)// factorial(n - r + 1) * 2 != len(cocircuit):
		#raise Exception("number of cocircuits is not correct!")
	G = nx.Graph()
	G.add_nodes_from(range(len(cocircuit)))
	for i in range(len(cocircuit)):
		neighbor_count = 0
		x = cocircuit[i]
		for j in range(i + 1, len(cocircuit)):
			#if neighbor_count == 2 * (r - 1):
				#break
			y = cocircuit[j]
			sym = [z < 0 for z in x * y]
			common_zeros = np.logical_and(x == 0, y == 0)
			if True not in sym and sum(common_zeros) >= r - 2:
				G.add_edge(i, j)
				neighbor_count += 1
	return G

def construct_graph(string, n, r):
	f = chirotope(n, r, string)
	vertices = cocircuits(n, r, f)
	G = make_graph(n, r, vertices)
	return G


def find_diameter(string, n, r):
	G = construct_graph(string, n, r)
	return diameter(G)

def find_pairs(G):
	d = diameter(G)
	pairs = []
	num_nodes = nx.number_of_nodes(G)
	for i in range(num_nodes):
		for j in range(i + 1, num_nodes):
			if nx.algorithms.shortest_path_length(G, source = i, target = j) == d:
				pairs.append((i, j))
	return pairs

def gen_bin_comb(n):
	l = [-1 for _ in range(n)]
	yield l
	while True:
		i = n - 1
		while i >= 0:
			if l[i] == -1:
				l[i] = 1
				break
			else:
				l[i] = -1
			i -= 1
		else:
			return
		yield l

def check_tope(string, n, r):
	f = chirotope(n, r, string)
	vertices = cocircuits(n, r, f)
	G = make_graph(n, r, vertices)
	for vector in gen_bin_comb(n):
		multiplier = np.array(vector)
		non_negatives = []
		counting = np.ones(n)
		for i in range(len(vertices)):
			v = vertices[i]
			if (v * multiplier >= 0).all():
				non_negatives.append(i)
				counting = counting * (v * multiplier)
		if len(non_negatives) > 0:
			H = G.subgraph(non_negatives)
			diam = diameter(H)
			new_n = sum(counting == 0)
			if diam > new_n - r + 1:
				print(string)
				print(multiplier)
				print(non_negatives)
				print(new_n)
				print(diam)

def check_lps_conjecture(G, vertices, n, r):
	for vector in gen_bin_comb(n):
		multiplier = np.array(vector)
		non_negatives = []
		for i in range(len(vertices)):
			v = vertices[i]
			if (v * multiplier >= 0).all():
				non_negatives.append(i)
		if len(non_negatives) > 0:
			H = G.subgraph(non_negatives)
			for i in range(len(non_negatives) - 1):
				start_point = non_negatives[i]
				for j in range(i + 1, len(non_negatives)):
					end_point = non_negatives[j]
					tope_distance = nx.shortest_path_length(H, source = start_point, target = end_point)
					om_distance = nx.shortest_path_length(G, source = start_point, target = end_point)
					if tope_distance != om_distance:
						print(vector)
						print(i, j)
						print(tope_distance, om_distance)
						raise Exception("London-Paris-Shanghai Conjecture violated")
	print("London-Paris-Shanghai Conjecture Checked!")

def check_hirsch_tope(filename, n, r):
	with open(filename) as f:
		content = f.read().splitlines()
	for line in content:
		#print(line)
		if "0" in line:
			print("non-uniform cases now")
			break
		check_tope(line, n, r)


def check_conjecture_distance(G, string, i, j, n, r, distance, vertices):
	x = vertices[i]
	y = vertices[j]
	z = x * y
	same = np.sum(z == 1)
	if not (distance == n - r + 1 - same):
		print(string)
		path = nx.shortest_path(G, source = i, target = j)
		for items in path:
			print(vertices[items])
		print(distance)

def check_conjecture_steve(string, n, r):
	f = chirotope(n, r, string)
	vertices = cocircuits(n, r, f)
	G = make_graph(n, r, vertices)
	for i in range(len(vertices)):
		p = nx.shortest_path_length(G, source = i)
		if i % 2 == 0:
			for j in range(i + 2, len(vertices)):
				check_conjecture_distance(G, string, i, j, n, r, p[j], vertices)
		else:
			for j in range(i + 1, len(vertices)):
				check_conjecture_distance(G, string, i, j, n, r, p[j], vertices)


def makeEFM():
	vertices = []
	with open("EFM8.txt") as f:
		for line in f:
			inner = [int(elt.strip()) for elt in line.split(",")]
			v = np.array(inner)
			vertices.append(v)
			vertices.append(-v)
	return vertices

def EFMtopes():
	vertices = makeEFM()
	G = make_graph(8, 4, vertices)
	count = 0
	for vector in gen_bin_comb(8):
		multiplier = np.array(vector)
		non_negatives = []
		for i in range(len(vertices)):
			v = vertices[i]
			if (v * multiplier >= 0).all():
				non_negatives.append(i)
		if len(non_negatives) > 0:
			H = G.subgraph(non_negatives)
			print(count)
			print(nx.algorithms.planarity.check_planarity(H)[0])
			# nx.draw(H, with_labels = True)
			# plt.savefig("EFM/" + str(count) + ".png")
			# plt.close()
			count += 1


def makeRS():
	vertices = []
	with open("RS8.txt") as f:
		for line in f:
			inner = [int(elt.strip()) for elt in line.split(",")]
			v = np.array(inner)
			vertices.append(v)
			vertices.append(-v)
	return vertices

def RStopes():
	vertices = makeRS()
	G = make_graph(8, 4, vertices)
	count = 0
	for vector in gen_bin_comb(8):
		multiplier = np.array(vector)
		non_negatives = []
		for i in range(len(vertices)):
			v = vertices[i]
			if (v * multiplier >= 0).all():
				non_negatives.append(i)
		if len(non_negatives) > 0:
			H = G.subgraph(non_negatives)
			nx.draw(H, with_labels = True)
			print(count)
			print(nx.algorithms.planarity.check_planarity(H)[0])
			# plt.savefig("RS/" + str(count) + ".png")
			# plt.close()
			count += 1

def check95():
	with open("result_short.txt") as f:
		content = f.read().splitlines()
	counts = 0
	print(len(content))
	print("Starting to compute for n=9, r=5 using circuits!")
	for items in content:
		f = chirotope(9, 4, items)
		vertices = circuits(9, 4, f)
		G = make_graph(9, 5, vertices)
		d = diameter(G)
		if d != 6:
			print("Something wrong with line " + str(counts))
			print(items)
			print(d)
			break
		counts += 1
		if counts % 100 == 0:
			print("Now processing line " + str(counts))

# check if a graph's shortest paths are all crabbed
def check_crabbed(string, n, r):
	f = chirotope(n, r, string)
	vertices = cocircuits(n, r, f)
	G = make_graph(n, r, vertices)
	num_nodes = nx.number_of_nodes(G)
	for i in range(num_nodes):
		for j in range(i + 1, num_nodes):
			x = vertices[i]
			y = vertices[j]
			check_indices = np.where((x * y) >= 0)[0]
			is_crabbed = False
			all_paths = nx.all_shortest_paths(G, source = i, target = j)
			for p in all_paths:
				if is_crabbed:
					break
				for counter in range(1, len(p)):
					if counter == len(p) - 1:
						is_crabbed = True
					else:
						new1 = x[check_indices] * vertices[p[counter]][check_indices]
						new2 = y[check_indices] * vertices[p[counter]][check_indices]
						if ((-1) in new1) or ((-1) in new2):
							break
			if not is_crabbed:
				print("conjecture is wrong for")
				print(string)
				print(i, j)
				return False
	return True



def check_antipodal(string, n, r):
	G = construct_graph(string, n, r)
	pair = find_pairs(G)
	return (len(G) == 2 * len(pair))

def check_antipodal_file(filename, n, r):
	with open(filename) as f:
		content = f.read().splitlines()
	for items in content:
		if "0" in items:
			print("non-uniform cases now")
			break
		if not check_antipodal(items, n, r):
			print("counterexamples here")
			print(items)
			break


def write_ilan_file(string, n, r, filename):
	f = chirotope(n, r, string)
	v = cocircuits(n, r, f)
	file = open(filename, "w")
	for items in v:
		items[items == -1] = 2
		for ele in items:
			file.write(str(ele) + " ")
		file.write("\n")
	return

# read ilan txt file, returns an array of cocircuits
def convert_ilan(filename):
	with open("Ilan_mapped_OMs/" + filename) as f:
		content = f.read().splitlines()
	v = []
	for items in content:
		x = np.fromstring(items, dtype = int, sep = " ")
		x[x == 2] = -1
		v.append(x)
	return v


# check if the array of cocircuits can be achieved from chirotope f
def check_ilan_chirotope(f, cocircuits, n, r):
	viewed_set = set()
	negative_set = set()
	for v in cocircuits:
		nonzero_index = np.nonzero(v)[0]
		zeroset = np.where(v == 0)[0]
		if len(zeroset) != r - 1:
			raise Exception("Incorrect number of zeros")
		signature = []
		for i in range(n):
			if i in zeroset:
				signature.append(0)
			else:
				exchange = sum([i < x for x in zeroset])
				new_set = list(zeroset)
				new_set.append(i)
				signature.append(f[frozenset(new_set)] * ((-1) ** exchange))
		for i in negative_set:
			signature[i] = signature[i] * (-1)
		p_cocircuit = np.array(signature)
		n_cocircuit = -p_cocircuit
		if np.array_equal(v[list(viewed_set)], p_cocircuit[list(viewed_set)]):
			for i in nonzero_index:
				if i not in viewed_set:
					viewed_set.add(i)
					if v[i] != p_cocircuit[i]:
						negative_set.add(i)
		elif np.array_equal(v[list(viewed_set)], n_cocircuit[list(viewed_set)]):
			for i in nonzero_index:
				if i not in viewed_set:
					viewed_set.add(i)
					if v[i] != n_cocircuit[i]:
						negative_set.add(i)
		else:
			print(viewed_set)
			print(negative_set)
			print(v)
			print(p_cocircuit)
			return False
	print(negative_set)
	return True


# check if ilan's file is in one of the OM files 
def check_ilan_isOM(OM_file, ilan_file, n, r):
	v = convert_ilan(ilan_file)
	with open(OM_file) as f:
		content = f.read().splitlines()
	for item in content:
		f = chirotope(n, r, item)
		if check_ilan_chirotope(f, v, n, r):
			print(item)
			return True
	print("Not one of them!")
	return False

# A is the matrix Ax=0, n is number of hyperplanes, r is the rank of OM
# A should be of size n*r, every r-1 hyperplanes give a cocircuit
def cocircuits_from_arrangement(A, n, r):
	cocircuits = []
	for item in findsubsets(n, r - 1):
		ns = null_space(A[item])
		ns = ns.T[0]
		v = A.dot(ns)
		v[np.abs(v) < 1e-6] = 0
		if np.count_nonzero(v) != n - r + 1:
			print(item)
			raise Exception("Error in number of zeros")
		new = np.sign(v)
		cocircuits.append(new)
		cocircuits.append(-new)
	return cocircuits


def main():
	if len(sys.argv) == 1:
		return 
	filename = sys.argv[1]
	n = int(sys.argv[2])
	r = int(sys.argv[3])
	mainfile = open(filename+"_diameters.txt", "w")
	file_small = open(filename+"_smaller.txt", "w")
	file_counter = open(filename+"_counterexamples.txt", "w")
	with open(filename + ".txt") as f:
		content = f.read().splitlines()
	if len(content[0]) != factorial(n) // factorial(r) // factorial(n - r):
		raise Exception("n choose r is not equal to number of chirotopes!")
	for items in content:
		print(items)
		d = find_diameter(items, n, r)
		mainfile.write(items + ":"+str(d)+"\n")
		if d < n - r + 2:
			file_small.write(items+":"+str(d)+"\n")
		if d > n -r + 2:
			file_counter.write(items+":"+str(d)+"\n")
	mainfile.close()
	file_small.close()
	file_counter.close()

if __name__ == "__main__":
	main()
