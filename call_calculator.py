from ase.io import read
from quippy.potential import Potential
from ase.io.extxyz import key_val_dict_to_str
import pdb

# example of how to call molpro calculator

template_path='/home/eg475/molpro_stuff/driver/template_e_f.inp'
calc_args = { 'template' : f'{{{template_path}}}',
              'molpro' : '{/opt/molpro/bin/molprop}',
              'energy_from' : 'RKS',
              # 'append_lines' : None,
              'test_mode' : True,
              # 'working_dir' : {/scratch-ssd/eg475/tmp},
              'extract_forces' : True}

with open(template_path, 'r') as f:
    print('Molpro template:')
    for line in f.readlines():
        print(line.rstrip())

from molpro import Molpro
print('Molpro imported')

methane = read('methane.xyz')
methane.set_calculator(Molpro)
print('calculator set')
energy = methane.get_potential_energy()
print(f'finally calculated energy: {energy}')
