# problem : there are three boxes each with a number. select exactly 2 boxes such that the sum of the numbers in the two boxes is minimum.



from pyqubo import Binary, Constraint
import neal

a, b, c = Binary('a'), Binary('b'), Binary('c')

M = 50.0

H = (17*a + 21*b + 19*c) + M*Constraint((a+b+c-2)**2,label='a+b+c=2')

model = H.compile()

sampler = neal.SimulatedAnnealingSampler()

bqm = model.to_bqm()

sampleset = sampler.sample(bqm, num_reads=10)

decoded_samples = model.decode_sampleset(sampleset)

best_sample = min(decoded_samples, key=lambda x: x.energy)

sol = {k:best_sample.sample[k] for k in sorted(best_sample.sample.keys())}

print(sol)