import psi4
import optking

# Test optimization with explicit loop using externally provided gradients
# In this case, using Leonard-Jones potential
h2o = psi4.geometry("""
     O
     H 1 1.0
     H 1 1.0 2 104.5
""")

psi4_options = {
  'basis': 'sto-3g',
  'g_convergence': 'gau_verytight',
  'intrafrag_step_limit_max': 0.2
}
psi4.set_options(psi4_options)

# Create an OptHelper object.  The comp_type == 'user' implies that the
# user will provide the gradients directly.
opt = optking.OptHelper('hf', comp_type='user')
opt.build_coordinates()

for step in range(30):
    # Compute one's own energy and gradient
    E, gX = optking.lj_functions.calc_energy_and_gradient(opt.geom, 2.5, 0.01, True)

    # Insert these values into the 'user' computer.
    opt.E = E
    opt.gX = gX
    opt.energy_gradient_hessian()

    # Take step
    opt.step()
    conv = opt.test_convergence()
    if conv is True:
        print("Optimization SUCCESS:")
        break
else:
    print("Optimization FAILURE:\n")

json_output = opt.close()

assert conv is True
E = json_output['energies'][-1]  # TEST
RefEnergy = -0.03  # - epsilon * 3, where -epsilon is depth of each Vij well
assert psi4.compare_values(RefEnergy, E, 6, "L-J Energy upon optimization")

