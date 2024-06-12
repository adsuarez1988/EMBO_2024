import msprime
import demesdraw
from matplotlib import pyplot as plt

## creating both populations ##
demography = msprime.Demography()
demography.add_population(name="N1", initial_size=200_000)
demography.add_population(name="N2", initial_size=6_600)

## creating initial population ##
demography.add_population(name="Nacn", initial_size=7_000_000)

## migracion ##

demography.set.migration_rate(source="N1", dest="N2", rate=0.1)
demography.set.migration_rate(source="N2", dest="N1", rate=0.1)

## setting split time ##
demography.add_population_split(time=1000, derived=["N1", "N2"], ancestral="Nacn") # split 1k generations ago

# Plot a schematic of the model
demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
plt.show()


ts = msprime.sim_ancestry(
        {"N1": 50, "N2": 50}, 
        demography=demography, 
        recombination_rate=3.5e-9, 
        sequence_length=1_000,
        random_seed=1234)
print(ts)

#### mutation rate ####

mts = msprime.sim_mutations(ts, rate=8.4e-9, random_seed=1234)
print(mts.tables.sites)




