from ase.io import read
from quippy.potential import Potential
from ase.io.extxyz import key_val_dict_to_str
import pdb

# example of how to call molpro calculator

calc_args = { 'template' : '{/home/eg475/molpro_stuff/driver/quippy_test/template_e_f.inp}',
              'molpro' : '{/opt/molpro/bin/molprop}',
              'energy_from' : 'RKS',
              # 'append_lines' : None,
              # 'test_mode' : False,
              # 'working_dir' : {/scratch-ssd/eg475/tmp},
              'extract_forces' : True}

molpro = Potential(args_str='FilePot command=/home/eg475/molpro_stuff/driver/molpro_driver.py', calc_args=calc_args)

methane = read('methane.xyz')
methane.set_calculator(molpro)
print('finally calculated energy:', methane.get_potential_energy())
