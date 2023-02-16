# Split the graph into two subsets such that:
#   1.The subsets are the same size.
#   2.There are as few edges as possible between the two subsets




from pyqubo import Binary, Constraint, Array
import neal
import networkx as nx
import matplotlib.pyplot as plt


G = nx.Graph()

G.add_edges_from([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])

x = Array.create('x',shape = 4, vartype='BINARY')

obj_fun = 0

for (i,j) in G.edges():
    obj_fun = x[i] + x[j] - 2*x[i]*x[j]

cons = (sum(x)-2)**2

M = 50

H = obj_fun + M*cons

model = H.compile()

sampler = neal.SimulatedAnnealingSampler()

bqm = model.to_bqm()

sampleset = sampler.sample(bqm, num_reads=10)

decoded_samples = model.decode_sampleset(sampleset)

best_sample = min(decoded_samples, key=lambda x: x.energy)


S = [best_sample.sample[k] for k in sorted(best_sample.sample.keys())]
sol_nodes = [i for i in range(4) if S[i]]
print(sol_nodes)
print({k:best_sample.sample[k] for k in sorted(best_sample.sample.keys())})

k = G.subgraph(sol_nodes)
notS = list(set(G.nodes()) - set(sol_nodes))
othersubgraph = G.subgraph(notS)
pos = nx.spring_layout(G)
plt.figure()

nx.draw_networkx(G, pos=pos, with_labels=True)
solution_name = "graph_part_plot_solution.png"
nx.draw_networkx(k, pos=pos, with_labels=True, node_color='r', font_color='k')
nx.draw_networkx(othersubgraph, pos=pos, with_labels=True, node_color='b', font_color='w')
plt.savefig(solution_name, bbox_inches='tight')