︠6934c1f5-7b05-4515-ac67-91609890ea44︠
import re
BaseRing = IntegerRing()
#BaseRing = RationalField()
#BaseRing = FiniteField(2)

def extract_numbers_from_string(text):
    numbers = re.findall(r'-?\d+', text)
    return [int(num) for num in numbers]

def get_key_from_value(dictionary, target_value):
    for key, value in dictionary.items():
        if value == target_value:
            return key
    return None  # Or raise an exception, depending on desired behavior

def magnitude_homology(g, lmax=7):
    kmax = lmax + 1
    d = g.distance_all_pairs()
    vertices = g.vertices(sort=True)

    # populate the generators recursively
    generators = dict(((s, t, k, l),[]) for s in vertices for t in vertices
                                     for k in range(kmax + 2) for l in range(lmax + 1))


    def add_generators(a, l, x):
        k = len(a) - 1
        if k <= kmax and l <= lmax:
            generators[(a[0], a[len(a) - 1], k, l)].append(a)
            for y in vertices:
                #if x != y and y not in a:        # uncomment this line for Eulerian magnitude homology
                #if x != y:                       # uncomment this line for ordinary magnitude homology
                    add_generators(a + [y], l + d[x][y], y)


    for x in vertices:
        add_generators([x], 0, x)

    # number the generators, so as to produce differentials rapidly
    for s in vertices:
        for t in vertices:
            for l in range(lmax + 1):
                for k in range(kmax + 1):
                    generators[(s, t, k, l)] = dict((tuple(a), i)
                                                    for (i, a) in enumerate(generators[(s, t, k, l)]))


    def differential(s, t, k, l):
        m = {}
        h = generators[(s, t, k-1, l)]
        for (a, i) in generators[(s, t, k, l)].items():
            for z in range(len(a) - 2):
                if d[a[z]][a[z + 1]] + d[a[z + 1]][a[z + 2]] == d[a[z]][a[z + 2]]:
                    j = h[a[:z + 1]+a[z + 2:]]
                    if z % 2:
                        m[(j, i)] = m.get((j, i), 0) + 1
                    else:
                        m[(j, i)] = m.get((j, i), 0) - 1
        return matrix(BaseRing, len(h), len(generators[(s, t, k, l)]), m)


    def chains(s, t, l):
        differentials = dict((k, differential(s, t, k, l)) for k in range(1, kmax + 1)
                             if generators[(s, t, k, l)] or generators[(s, t, k-1, l)])
        return ChainComplex(differentials, base_ring=BaseRing, degree=-1)


    def homology(s, t, l):
        return chains(s, t, l).homology(generators=False)


    return dict(((s, t, l), homology(s, t, l)) for s in vertices for t in vertices
                                               for l in range(lmax+1))

def magnitude_homology_generators(g, lmax, k_input, vertex_1, vertex_2):
    d = g.distance_all_pairs()
    vertices = g.vertices(sort=True)

    # populate the generators recursively
    generators = dict(((s, t, k, lmax),[]) for s in vertices for t in vertices for k in range(k_input+2))

    def add_generators(a, l, x):
        k = len(a) - 1
        if k < k_input+2 and l <= lmax:
            for y in vertices:
                #if x != y and y not in a:           # uncomment this line for Eulerian magnitude homology
                #if x != y:                          # uncomment this line for ordinary magnitude homology
                    add_generators(a + [y], l + d[x][y], y)
        if k < k_input+2 and l == lmax:
            generators[(a[0], a[len(a) - 1], k, l)].append(a)

    for x in vertices:
        add_generators([x], 0, x)

    # number the generators, so as to produce differentials rapidly
    for s in vertices:
        for t in vertices:
            for k in range(k_input + 2):
                generators[(s, t, k, lmax)] = dict((tuple(a), i)
                    for (i, a) in enumerate(generators[(s, t, k, l)]))

    #print(generators[(vertex_1, vertex_2, k_input, lmax)])

    def differential(s, t, k, l):
        m = {}
        h = generators[(s, t, k-1, l)]
        q = generators[(s, t, k, l)]
        for (a, i) in q.items():
            for z in range(len(a) - 2):
                if d[a[z]][a[z + 1]] + d[a[z + 1]][a[z + 2]] == d[a[z]][a[z + 2]]:
                    j = h[a[:z + 1]+a[z + 2:]]
                    if z % 2:
                        m[(j, i)] = m.get((j, i), 0) + 1
                    else:
                        m[(j, i)] = m.get((j, i), 0) - 1
        return matrix(BaseRing, len(h), len(generators[(s, t, k, l)]), m)

    def chains(s, t, l):
        differentials = dict((k, differential(s, t, k, l)) for k in range(k_input-1, k_input+2)
                             if generators[(s, t, k, l)] or generators[(s, t, k-1, l)])
        return ChainComplex(differentials, base_ring=BaseRing, degree=-1)

    def homology(s, t, l):
        return chains(s, t, l).homology(generators=True)

    gens = homology(vertex_1, vertex_2, lmax)[k_input]
    gen_word = str(gens[0][1])[8:]
    gen_array = extract_numbers_from_string(gen_word)
    for i in range(len(gen_array)):
        if gen_array[i] != 0:
            a = get_key_from_value(generators[(vertex_1, vertex_2, k_input, lmax)], i)
            print(gen_array[i], a)

#RP2
#g = Graph([
#   (0,1),(0,2),(0,3),(0,4),(0,5),(0,6),
#   (1,7),(1,8),(1,10),(1,13),(1,17),
#   (2,7),(2,9),(2,11),(2,14),(2,18),
#   (3,8),(3,9),(3,12),(3,15),(3,19),
#   (4,10),(4,11),(4,12),(4,16),(4,20),
#   (5,13),(5,14),(5,15),(5,16),(5,21),
#   (6,17),(6,18),(6,19),(6,20),(6,21),
#   (7,22),(7,26),
#   (8,23),(8,24),
#   (9,25),(9,27),
#   (10,22),(10,23),
#   (11,22),(11,29),
#   (12,23),(12,30),
#   (13,24),(13,28),
#   (14,25),(14,29),
#   (15,24),(15,25),
#   (16,29),(16,31),
#   (17,26),(17,28),
#   (18,26),(18,27),
#   (19,27),(19,30),
#   (20,30),(20,31),
#   (21,28),(21,31),
#   (22,32),(23,32),(24,32),(25,32),(26,32),(27,32),(28,32),(29,32),(30,32),(31,32)
#]); graph_name = 'RP^2'

#g = Graph([
#    (100,0),
#    (0,1),(0,2),(0,3),(0,4),
#    (1,5),(1,6),(1,8),#(1,11), keeps torsion
#    (2,5),(2,7),(2,9),
#    (3,6),(3,7),(3,10),
#    (4,8),(4,9),(4,10),
#    (1,14),(2,14),(3,14),(4,14),
#    (5,11),(6,12),(7,13),(8,13),(9,12),(10,11),
#    (5,13),(6,11),(7,12),(8,12),(9,11),(10,13),
#    (11,14),(12,14),(13,14)
#]); graph_name = 'P_4^(1,2,3)'

# Altered Boolean Poset for B_4
#g = Graph([
#    (0,1),(0,2),(0,3),(0,4),
#    (1,5),(1,6),(1,7),
#    (2,5),(2,8),(2,9),
#    (3,6),(3,8),(3,10),
#    (4,7),(4,9),(4,10),
#    (5,11),(6,11),(7,12),(8,13),(9,11),(10,12),
#    (5,12),(6,13),(7,13),(8,12),(9,13),(10,11),
#    (11,14),(12,14),(13,14)
#]); graph_name= 'graph'

#g = Graph([
#    (0,1),(0,2),(0,3),(0,4),(0,5),
#    (1,6),(1,7),(1,8),(1,9),
#    (2,6),(2,10),(2,11),(2,12),
#    (3,7),(3,10),(3,13),(3,14),
#    (4,8),(4,11),(4,13),(4,15),
#    (5,9),(5,12),(5,14),(5,15),
#    (6,16),(6,17),(6,18),
#    (7,16),(7,19),(7,20),
#    (8,17),(8,19),(8,21),
#    (9,18),(9,20),(9,21),
#    (10,16),(10,22),(10,23),
#    (11,17),(11,22),(11,24),
#    (12,18),(12,23),(12,24),
#    (13,19),(13,22),(13,25),
#    (14,20),(14,23),(14,25),
#    (15,21),(15,24),(15,25),
#    (16,26),(17,26),(19,26),(23,26),(24,26),
#    (16,27),(18,27),(20,27),(25,27),(22,27),
#    (17,28),(18,28),(21,28),(24,28),(25,28),
#    (19,29),(20,29),(21,29),(22,29),(23,29),
#    (26,30),(27,30),(28,30),(29,30)
#])

# Regular CW structure on (2,1)-Lens space
#g = Graph([
#    (-1,0),(-1,1),(-1,2),(-1,5),
#    (0,10),(0,20),(0,30),(0,40),(0,50),(0,60),
#    (1,10),(1,30),(1,15),(1,35),(1,12),(1,23),
#    (2,20),(2,40),(2,12),(2,23),(2,25),(2,45),
#    (5,50),(5,60),(5,15),(5,35),(5,25),(5,45),
#    (10,105),(10,106),(10,102),(10,104),
#    (20,205),(20,206),(20,102),(20,203),
#    (30,305),(30,306),(30,203),(30,304),
#    (40,405),(40,406),(40,104),(40,304),
#    (50,105),(50,305),(50,205),(50,405),
#    (60,106),(60,306),(60,206),(60,406),
#    (15,125),(15,236),(15,105),(15,306),
#    (35,235),(35,126),(35,305),(35,106),
#    (12,125),(12,126),(12,102),(12,304),
#    (23,235),(23,236),(23,203),(23,104),
#    (25,125),(25,235),(25,205),(25,406),
#    (45,126),(45,236),(45,405),(45,206),
#    (1025,125),(1025,102),(1025,105),(1025,205),
#    (1045,236),(1045,104),(1045,105),(1045,405),
#    (2035,235),(2035,203),(2035,205),(2035,305),
#    (3045,126),(3045,304),(3045,305),(3045,405),
#    (1026,126),(1026,102),(1026,106),(1026,206),
#    (1046,235),(1046,104),(1046,106),(1046,406),
#    (2036,236),(2036,203),(2036,206),(2036,306),
#    (3046,125),(3046,304),(3046,306),(3046,406),
#    (1025,-2),(1045,-2),(2035,-2),(3045,-2),(1026,-2),(1046,-2),(2036,-2),(3046,-2)
#]); graph_name = 'graph'

# Graph with a = 5; b = 10; c = 4
#No Torsion
#g = Graph([
#    (0,1),(0,2),(0,3),(0,4),(0,5),
#    (1,6),(1,7),(1,9),(1,11),
#    (2,6),(2,8),(2,10),(2,12),
#    (3,7),(3,8),(3,13),(3,14),
#    (4,10),(4,11),(4,13),(4,15),
#    (5,9),(5,12),(5,14),(5,15),
#    (6,16),(7,17),(8,18),(9,18),(10,17),(11,19),(12,19),(13,16),(14,16),(15,18),
#    (6,17),(7,19),(8,16),(9,16),(10,19),(11,18),(12,17),(13,17),(14,18),(15,19),
#    (16,20),(17,20),(18,20),(19,20)
#]); graph_name = 'this graph'

# a = 6; b = 15; c = 5
#g = Graph([
#    (0,1),(0,2),(0,3),(0,4),(0,5),(0,6),
#    (1,7),(1,8),(1,10),(1,11),(1,14),
#    (2,7),(2,9),(2,13),(2,15),(2,17),
#    (3,8),(3,9),(3,12),(3,16),(3,18),
#    (4,10),(4,12),(4,13),(4,19),(4,20),
#    (5,11),(5,15),(5,16),(5,19),(5,21),
#    (6,14),(6,17),(6,18),(6,20),(6,21),
#    (7,26),(8,23),(9,25),(10,22),(11,25),(12,26),(13,24),(14,24),(15,22),(16,24),(17,23),(18,22),(19,23),(20,25),(21,26),
#    (7,22),(8,24),(9,23),(10,23),(11,26),(12,25),(13,26),(14,25),(15,25),(16,22),(17,24),(18,26),(19,24),(20,22),(21,23),
#    (22,27),(23,27),(24,27),(25,27),(26,27)
#])

# Regular CW structure on 3-Moore space
#g = Graph([
#    (-1,0),(-1,1),(-1,2),
#    (0,10),(0,20),(0,30),(0,40),(0,50),(0,60),
#    (1,10),(1,30),(1,50),(1,12),(1,23),
#    (2,20),(2,40),(2,60),(2,12),(2,23),
#    (102,10),(102,20),(102,12),
#    (203,20),(203,30),(203,23),
#    (304,30),(304,40),(304,12),
#    (405,40),(405,50),(405,23),
#    (506,50),(506,60),(506,12),
#    (106,10),(106,60),(106,23),
#    (102,-2),(203,-2),(304,-2),(405,-2),(506,-2),(106,-2)
#]); graph_name = 'this graph'

# Regular CW struture on L(3,1)
#g = Graph([
#    (-1,0),(-1,1),(-1,2),(-1,7),
#    (0,10),(0,20),(0,30),(0,40),(0,50),(0,60),(0,70),(0,80),
#    (1,10),(1,30),(1,50),(1,12),(1,23),(1,17),(1,37),(1,57),
#    (2,20),(2,40),(2,60),(2,12),(2,23),(2,27),(2,47),(2,67),
#    (7,70),(7,80),(7,17),(7,27),(7,37),(7,47),(7,57),(7,67),
#    (102,10),(102,20),(102,12),
#    (203,20),(203,30),(203,23),
#    (304,30),(304,40),(304,12),
#    (405,40),(405,50),(405,23),
#    (506,50),(506,60),(506,12),
#    (106,10),(106,60),(106,23),
#    (107,10),(107,70),(107,17),
#    (207,20),(207,70),(207,27),
#    (307,30),(307,70),(307,37),
#    (407,40),(407,70),(407,47),
#    (507,50),(507,70),(507,57),
#    (607,60),(607,70),(607,67),
#    (108,10),(108,80),(108,57),
#    (208,20),(208,80),(208,67),
#    (308,30),(308,80),(308,17),
#    (408,40),(408,80),(408,27),
#    (508,50),(508,80),(508,37),
#    (608,60),(608,80),(608,47),
#    (127,12),(127,17),(127,27),
#    (237,23),(237,27),(237,37),
#    (347,12),(347,37),(347,47),
#    (457,23),(457,47),(457,57),
#    (567,12),(567,57),(567,67),
#    (167,23),(167,17),(167,67),
#    (1027,102),(1027,107),(1027,207),(1027,127),
#    (2037,203),(2037,207),(2037,307),(2037,237),
#    (3047,304),(3047,307),(3047,407),(3047,347),
#    (4057,405),(4057,407),(4057,507),(4057,457),
#    (5067,506),(5067,507),(5067,607),(5067,567),
#    (1067,106),(1067,107),(1067,607),(1067,167),
#    (1028,102),(1028,108),(1028,208),(1028,567),
#    (2038,203),(2038,208),(2038,308),(2038,167),
#    (3048,304),(3048,308),(3048,408),(3048,127),
#    (4058,405),(4058,408),(4058,508),(4058,237),
#    (5068,506),(5068,508),(5068,608),(5068,347),
#    (1068,106),(1068,108),(1068,608),(1068,457),
#    (1027,-2),(1028,-2),(2037,-2),(2038,-2),(3047,-2),(3048,-2),(4057,-2),(4058,-2),
#    (5067,-2),(5068,-2),(1067,-2),(1068,-2)
#]); graph_name = 'L(3,1)'

# Whitney Twist graph 1
#g = Graph([
#    (1,2),(2,3),(3,1),(1,4),(1,5)
#]); graph_name = 'this graph'

# Whitney Twist graph 2
#g = Graph([
#    (1,2),(2,3),(3,1),(1,4),(3,5)
#]); graph_name = 'this graph'

#g = graphs.CycleGraph(13); graph_name='this graph'

g.show('this graph')
lmax = 10 #what do you want the maximum length to be

print(graph_name)
print('lmax = {0}'.format(lmax))

homology = magnitude_homology(g, lmax)

total_rank = dict(((k, l), 0) for k in range(0, lmax + 1) for l in range(0, lmax + 1))


vertices = g.vertices(sort=True)

for s in vertices:
    for t in vertices:
        for l in range(lmax + 1):
            for degree, group in sorted(homology[s, t, l].items()):
                for i in range(len(group.invariants())):
                    if group.invariants()[i] == 0:
                        total_rank[degree, l] += 1
                    else:
                        group
                        print (s,t,degree,l)
                        magnitude_homology_generators(g, l, degree, s, t)

for l in range(0, lmax + 1):
    print ('\n\r{0:2d}: '.format(l), end="|")
    for k in range(0, lmax + 1):
        if total_rank[k, l] != 0:
            print ('{0:2d}'.format(total_rank[k, l]), end="  |")
        else:
            print ('    ', end="|")
︡e2d9867d-8d25-4697-b5e5-aaa7f51296ce︡{"file":{"filename":"/tmp/tmpd3dq29xw/tmp_hzmuo_qi.svg","show":true,"text":null,"uuid":"a4d3de25-45b6-485f-8d58-c7a81b0947f1"},"once":false}︡{"stdout":"this graph\n"}︡{"stdout":"lmax = 10\n"}︡{"stdout":"\n\r 0: |13  |    |    |    |    |    |    |    |    |    |    |\n\r 1: |    |26  |    |    |    |    |    |    |    |    |    |\n\r 2: |    |    |    |    |    |    |    |    |    |    |    |\n\r 3: |    |    |52  |    |    |    |    |    |    |    |    |\n\r 4: |    |    |    |    |    |    |    |    |    |    |    |\n\r 5: |    |    |    |104  |    |    |    |    |    |    |    |\n\r 6: |    |    |    |52  |    |    |    |    |    |    |    |\n\r 7: |    |    |26  |26  |156  |    |    |    |    |    |    |\n\r 8: |    |    |    |78  |156  |    |    |    |    |    |    |\n\r 9: |    |    |    |    |286  |234  |    |    |    |    |    |\n\r10: |    |    |    |    |390  |312  |    |    |    |    |    |"}︡{"done":true}
︠ff978ed0-0935-4bef-88de-fb09d9c1460d︠











