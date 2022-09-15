from cobra import Model, Reaction, Metabolite
import cobra
from cobra.io import load_model


model = load_model("iMM904")

reaction = Reaction('iMM904')
reaction.name = 'Saccharomyces cerevisiae S288C'
reaction.subsystem = ''
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

#accoa_m = Metabolite(
#    'accoa_m',
#    formula='C23H34N7O17P3S',
#    name='Acetyl-CoA',
#    compartment='c'
#)

print(model.metabolites.accoa_m.summary())