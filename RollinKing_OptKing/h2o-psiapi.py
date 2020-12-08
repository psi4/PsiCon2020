# Demonstration of python optking using psi api.
# See optking/tests for lots of psi api examples.
import psi4
import optking

psi4.core.set_output_file('psi-out.dat')
h2o = psi4.geometry("""
     O
     H 1 1.0
     H 1 1.0 2 104.5
""")
psi4_options = { 'basis': 'cc-pvdz' }
psi4.set_options(psi4_options)

json_output = optking.optimize_psi4('hf')

print('Contents of QCSchema output:')
for k in json_output.keys():
    print('\t' + k)

print('Properties of (next to) last geom:')
for k in json_output['trajectory'][-1]['properties']:
    print('\t' + k)

print('Extras of (next to) last geom:')
for k in json_output['trajectory'][-1]['extras']['qcvars']:
    print('\t' + k)

print('Energy {:15.10f}'.format(json_output['energies'][-1]))
print('NRE    {:15.10f}'.format(json_output['trajectory'][-1]['properties']['nuclear_repulsion_energy']))
