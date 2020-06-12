import sys
import time
import random
sys.path.append('/home/eg475/molpro_stuff/driver')
from ase.io import read
from molpro import Molpro



# example of how to call molpro calculator

template_path='/home/eg475/molpro_stuff/driver/example/template_e_f.inp'
# template_path = '/opt/project/example/template_e_f.inp'

calc_args = { \
              'template': template_path,
              'molpro' : '/opt/molpro/bin/molpro',
              'energy_from' : 'RKS',
              # 'append_lines' : None,
              #  'test_mode' : True,
              # 'working_dir' : {/scratch-ssd/eg475/tmp},
              'extract_forces' : True}

with open(template_path, 'r') as f:
    print('Molpro template:')
    for line in f.readlines():
        print(line.rstrip())


molpro=Molpro(calc_args=calc_args)
methane = read('methane.xyz')
methane.set_calculator(Molpro(calc_args=calc_args))
energy = methane.get_potential_energy()
forces = methane.get_forces()
print(f'Methane energy: {energy}')
print(f'And forces: \n{forces}')



print('timing more atoms')
start = time.time()
no_atoms = 4
energies = []
for _ in range(no_atoms):
    at = methane.copy()
    seed = random.randint(0, 100000000)
    print('seed', seed)
    at.rattle(stdev=0.01, seed=seed)
    at.set_calculator(Molpro(calc_args=calc_args))
    energy=at.get_potential_energy()
    energies.append(energy)
    print('Energy:', energy)

print(f'Time: {(time.time()-start)/no_atoms} s/structure; {no_atoms} structures') 


