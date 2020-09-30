from molpro import Molpro
from ase.build import molecule
from ase.io import read, write
from pathlib import Path
import os
import shutil
from ase.units import Ha
import unittest
import math
import pytest
from ase.calculators.calculator import SCFError
import time

def test_write_input():
    '''creates a molecule and prepares molpro template for it.
    Does not cover every clause in the function though. '''
    at = molecule('CH4')
    mp = Molpro()
    mp.write_input(atoms=at)

    if not os.path.isfile('MOLPRO/molpro.xyz'):
        raise FileNotFoundError('.xyz file was not created')

    if not os.path.isfile('MOLPRO/molpro.ase'):
        raise FileNotFoundError('.ase parameter file was not found')

    expected_input = 'geomtyp=xyz\n' \
                     'geom=MOLPRO/molpro.xyz\n' \
                     'basis=6-31G*\n' \
                     'rks, b3lyp\n' \
                     'force\n'

    read_input = Path('MOLPRO/molpro.inp').read_text()
    assert expected_input == read_input
    # shutil.rmtree('MOLPRO')

    # test setting max no iterations
    mp = Molpro(maxit=300)
    mp.write_input(atoms=at)
    expected_input = 'geomtyp=xyz\n' \
                     'geom=MOLPRO/molpro.xyz\n' \
                     'basis=6-31G*\n' \
                     'rks, b3lyp, maxit=300\n' \
                     'force\n'

    read_input = Path('MOLPRO/molpro.inp').read_text()
    assert expected_input == read_input
    # shutil.rmtree('MOLPRO')

    # test giving a command block
    command_block = '{rks, b3lyp\n' \
                    'wf, nelec=2, symmetry=1, spin=0\n'\
                    'occ, 1, 0, 0, 0, 0, 0, 0, 0}'

    mp = Molpro(command_block=command_block)
    mp.write_input(atoms=at)
    expected_input = 'geomtyp=xyz\n' \
                     'geom=MOLPRO/molpro.xyz\n' \
                     'basis=6-31G*\n' \
                     '{rks, b3lyp\n'\
                     'wf, nelec=2, symmetry=1, spin=0\n'\
                     'occ, 1, 0, 0, 0, 0, 0, 0, 0}\n'\
                     'force\n'

    read_input = Path('MOLPRO/molpro.inp').read_text()
    assert expected_input == read_input
    # shutil.rmtree('MOLPRO')

    # test making a template
    expected_input =  'geomtyp=xyz\n' \
                    'geom=MOLPRO/molpro.xyz\n' \
                     'basis=6-31G*\n' \
                     'rks, b3lyp, maxit=300\n' \
                     'force\n'

    old_template = 'geomtyp=xyz\n' \
                     'geom=\n' \
                     'basis=6-31G*\n' \
                     'rks, b3lyp, maxit=300\n'

    if not os.path.isdir('MOLPRO'):
        os.makedirs('MOLPRO')
    with open('MOLPRO/template.txt', 'w') as f:
        f.write(old_template)

    mp = Molpro(template_path='MOLPRO/template.txt')
    mp.write_input(atoms=at)
    read_input = Path('MOLPRO/molpro.inp').read_text()
    assert expected_input == read_input
    # shutil.rmtree('MOLPRO')



def test_getting_energies():
    at = molecule('CH4')
    mp = Molpro()
    # mp.run_molpro=False
    at.set_calculator(mp)
    e_start = time.time()
    energy = at.get_potential_energy()
    e_time = time.time() - e_start

    assert math.isclose(energy, -40.48163242 * Ha)

    f_start = time.time()
    at.get_forces()
    f_time = time.time() - f_start

    # getting forces should be just reading them out if task=='gradient' (default)
    assert e_time/2 > f_time




def test_catching_maxing_out_of_scf_iterations():
    at = molecule('CH4')
    mp = Molpro(maxit=3)
    # mp.run_molpro=False
    at.set_calculator(mp)
    with pytest.raises(SCFError) as execinfo:
        # Not Sure this is the right way to do this
        energy = at.get_potential_energy()








