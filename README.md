# OM2M
## Calculation of heterozygosity trajectories
Function that returns (dHs/dt, dHt/dt)
Assumptions: infinite alleles per site.
asserions: two lineages can coalesce only if they are in the same deme.
Derivation: backwards-coalescence frame of thought. 

Within-demes:
H(t+1) = Ht(1-1/2N-2m-2Miu) + 2mHgen(t) + 2miu(*), whose reordering gives dH/timestep,

As apears in the function below. [note that the case of coalescence doesn't have a
term of its own, because its contribution to H is zero]
* : If we assume an infinite allele model, any mutation always leads to a contribution 
to H. If we assume a 2-allele model with back mutation, we need to add the terms:
2miu(1-Ht) - 2miu*Ht, where the first is a mutation in a pair that were identical, 
and the second a back-mutation in an initially non-identical pair.
Hgeneral (HG): we assume no two events occur in a single timestep, so ignore the option
of migration and coalescence or mig+mutation, and we even ignore the option that both
individuals happen to be from the same deme and one of them had a mutation (although
this last thing isn't completely justified, I think).
HG(t+1) = HG(t)(1-1/L-2mu) + 1/L * H(t) + 2mu
HG(t+1)-HG(t) = 1/L(H(t)-HG(t)) + 2mu(1-HG(t))
Which looks, admittedly, quite different from Mike's derivation which is based on a
forward-in-time-in-his-head intuition/derivation-that-I-didn't-follow.
