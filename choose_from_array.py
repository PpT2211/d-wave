# You have the set A = {1, 2, 3, 4, 5}. Find 3 numbers in A that add up to 8.



from pyqubo import Binary, Constraint, Array
import neal

x = Array.create('x',shape = (5), vartype='BINARY')

M = 50.0

arr = [1, 2, 3, 4, 5]

eq = [arr[i]*x[i] for i in range(5)]

H = (sum(eq) - 8)**2 + M*Constraint((sum(x)-3)**2,label='sum(x)=3')

model = H.compile()

sampler = neal.SimulatedAnnealingSampler()

bqm = model.to_bqm()

sampleset = sampler.sample(bqm, num_reads=10)

decoded_samples = model.decode_sampleset(sampleset)

best_sample = min(decoded_samples, key=lambda x: x.energy)

sol = [best_sample.sample[k] for k in sorted(best_sample.sample.keys())]

final = []

for s in range(len(sol)):
    if sol[s]:
        final.append(arr[s])

print(final)