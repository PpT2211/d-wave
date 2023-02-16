from pyqubo import Binary, Constraint, Array
import neal

x = Array.create('x',shape = (5), vartype='BINARY')

M = 50.0

values = [1,2,3,4,5]
wt = [2,4,6,8,10]

eq = 0
cons = 0
for i in range(len(values)):
    eq -= values[i]*x[i]

for i in range(len(wt)):
    cons += wt[i]*x[i]


H = eq + M*Constraint((-cons+24)**2, label='sum(x)=3')

model = H.compile()

sampler = neal.SimulatedAnnealingSampler()

bqm = model.to_bqm()

sampleset = sampler.sample(bqm, num_reads=10)

decoded_samples = model.decode_sampleset(sampleset)

best_sample = min(decoded_samples, key=lambda x: x.energy)

sol = {k:best_sample.sample[k] for k in sorted(best_sample.sample.keys())}
print(sol)
