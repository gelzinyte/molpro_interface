"""The module defines an interface to MOLPRO 2012 (?).
    Based on Python 2 molpro driver for quippy by Alan Nichol, James Kermode
    and others (????). Adapted to ASE by Elena Gelžinytė.
    Loosely based on castep ase calculator and Tamas' Orca (based on NWChem?).
"""

# TODO list:
#   Raise appropriate calculator errors
#   Clean up imports to follow python conventions
#   can I leave scratch directory to its own devices if not specified?

import os
import pdb
import ase
from lxml import etree as et
from ase.calculators.calculator import FileIOCalculator, Parameters, \
                                       CalculatorSetupError, Calculator, \
                                       ReadError, PropertyNotImplementedError,\
                                    SCFError, PropertyNotPresent

__all__ = ['Molpro']



class Molpro(FileIOCalculator):
    r"""
Molpro Interface Documentation


Introduction
============

Something something


Getting Started
===============

something something 
something something


Running the Calculator
======================


Arguments:
==========

=========================  ===================================================
Keyword                    Description
=========================  ===================================================
``directory``              Relative path where all Molpro-related items will
                           be placed. Will be created if doesn't exist.
                           TODO: add behaviour with how to deal with existing
                           directory.

``label``                  The prefix for .inp, .xyz, .out, .xml, etc files.

``molpro_command``         Command to run molpro. Can also be set via the bash
                           environment variable ``MOLPRO_COMMAND`. If none is
                           given or found, will default to ``molpro``.

``check_molpro_version``   Should implement???

``keyword_tolerance``      Should implement?

=========================  ===================================================


Arguments:
==========

=========================  ===================================================
Keyword                    Description
=========================  ===================================================


=========================  ===================================================


Additional Settings
===================

=========================  ===================================================
Internal Setting           Description
=========================  ===================================================
``_molpro_command``        (``=molpro``): the actual shell command used to
                           call MOLPRO

``_prepare_input_only``     (``=False``): If set to True, the calculator will
                           create \*.inp and \*.xyz files but not start the
                           calculation itself.

``_rename_existing_dir``   (``=True``) : when using a new instance
                           of the calculator, this will move directories out of
                           the way that would be overwritten otherwise,
                           appending a date string.

``_set_atoms``             (``=False``) : setting this to True will overwrite
                           any atoms object previously attached to the
                           calculator when reading the output file.  By de-
                           fault, the read() function will only create a new
                           atoms object if none has been attached and other-
                           wise try to assign forces etc. based on the atom's
                           positions.  ``_set_atoms=True`` could be necessary
                           if one uses internal geometry optimization
                           (``blah``)
                           because then the positions get out of sync.
                           *Warning*: this option is generally not recommended
                           unless one knows one really needs it. There should
                           never be any need, if MOLPRO is used as a
                           single-point calculator.

=========================  ===================================================

Special features
================


``.keyword.clear()``
  ??maybe implement?

``.initialize()``
  Creates all needed input in the ``_directory``. This could then be copied to
  and ran in a place without ASE or Python.

``print(calc)``
  Prints a short summary of the calculator settings and atoms.


Notes/Issues
============

blah blah


End Molpro Interface Documentation
"""

    # Class attributes
    # TODO implement/check how this is used/called, now only what I think is ok
    atoms_keys = ['charges']

    atoms_obj_keys = [
        'dipole',
        'energy_free',
        'forces',
        'positions']

    internal_keys = [
        '_molpro_command',
        '_prepare_input_only',
        '_rename_existing_dir',
        '_set_atoms'
    ]
    # TODO check this!!!
    implemented_properties = ['energy', 'forces']
    # TODO maybe implement default_parameters?
    # default_parameters = dict(
    #     memory='500, m',
    #     basis='6-31G',
    #     program='hf',
    #     task='forces')

    def __init__(self, restart=None, ignore_bad_restart_file=False, label='MOLPRO/molpro',
               atoms=None, molpro_command=None, **kwargs):

        self.__name__='Molpro'

        # initialise the FileIOCalculator. This does:
        #   Initialises ase.calculators.calculator.Calculator
        #   sets self.command to be the given command.
        #
        # Initialising Calculator does:
        #    sets atoms, reads results from file if restart, sets label ~directory/prefix
        #    reads parameters or gets default parameters
        #    TODO set defaut parameters somewhere somehow
        #    and sets name.
        # pdb.set_trace()

        Calculator.__init__(self, restart=restart,
                ignore_bad_restart_file=ignore_bad_restart_file, label=label,
                atoms=atoms,  **kwargs)



        # Maybe I need this, maybe I don't
        # Calculator state variables
        self._calls = 0
        self._warnings = []
        self._error = None
        self._interface_warnings = []

        # Internal Keys to allow to tweak behaviour
        self.molpro_command = get_molpro_command(molpro_command)
        # self._force_write = True
        self._prepare_input_only = False
        self._rename.existing_dir = True
        self._set_atoms = False
        self._default_forces = False   #call energies and forces when calling calculator
        # default. Gets overwritten if set in parameters (TODO! Also does
        #  maxit only apply to scf or other optimisations?) or
        self._scf_maxit = 60


        # Physical result variables
        # check if I need this.
        self._point_group = None

        # TODO cycle through parameters and set appropriate ones here.
        # TODO capitalise appropriate parameters (?) 





    def calculate(self, atoms=None, properties=['energy'], systm_changes=\
                  ['positions', 'numbers', 'cell', 'pbc', 'initial_charges',
                   'initial_magmoms']):
        # TODO check what to do with properties and system_changes
        self.write_input(self, atoms, properties, system_changes)
        if not self._prepare_input_only:
            self.run()
            self.read()
            self.push_old_state()


    def run(self):
        self._calls += 1

        # add scratch directory
        # deal with setting workind directory and giving labels, rather than anything else
        stdout, stderr = shell_stdouterr(f'{self.molpro_command} {self.input_fname} -o {self.output_fname}')

        if stdout:
            print(f'molpro call stdout:\n{stdout}')

        if stderr:
            print(f'molpro call stderr:\n{stderr}')

    def push_old_state(self):
        pass


    # TODO difference between read and read_results?
    def read(self):
        # maybe introduce function to read open file or just generic file as with castep.
        # deal with reading from .out not only .xml, also from a fileobject
        if not os.path.isfile(self.label + '.xml'):
            # would this be ok speed and memory-wise?
            # xml = parse_xml()
            raise ReadError('xml output file does not exist')
        elif os.path.isfile(self.label + '.out'):
            raise RuntimeError('.xml is not present and reading from .out files not implemented yet')
        # else:


        self.parameters = Parameters.read(self.label + '.ase')
        # decide what to do with _set_atoms and
        # Check if 'overriding previous calculator-attached results' works like this
        if not self._set_atoms:
            self.atoms = ase.io.read(self.label + '.xyz')
        else:
            # Means used non-single point energy calculation and need to update atom positions.
            raise NotImplementedError('reading atoms from .xml or .out not yet implemented')

        self.read_results()

    # def parse_xml(self):
    #     check convergence
    #     if calc_type not given, get the last one set in the input file
    #       I guess either in the .inp or from the xml.
    def check_warnings_errors(self):
        root = et.parse(slef.label + '.xml').getroot()

    def check_SCF_convergence(self):
        # and read other stuff too
        program = self.parameters.program
        root = et.parse(slef.label + '.xml').getroot()
        job = root.find('{http://www.molpro.net/schema/molpro-output}job')
        jobsteps = job.findall(
            '{http://www.molpro.net/schema/molpro-output}jobstep')
        for jobstep in jobsteps:
            if program in jobstep.attrib['command']:
                for child in jobstep:
                    text = child.text
                    if text is not None and 'ITERATION' in text:
                        text_lines = text.splitlines()
                        #                 print(text)
                        for idx, line in enumerate(text_lines):
                            if 'MAX. NUMBER OF ITERATIONS' in line:
                                maxit = int(re.search(r'\d+', line).group())
                                print('maxit:', maxit)
                            if 'Final occupancy' in line:
                                last_scf_line = text_lines[idx - 2].strip()
                                # ^ and search are probably doubling up
                                match = re.search(r'^\d+',
                                                  last_scf_line).group()
                                last_iteration = int(match)
                                print(last_iteration)
                                if last_iteration == 11:
                                    print('yay')

                            # TODO Spot messages for not converging and
                            #  stopping after a few iterations

    def read_other_stuff(self):
        pass

    def read_results(self):
        # TODO do this somehow better? Definitely for multiple cases
        # check how to check for 'program' TODO
        if 'program' not in self.parameters.keys():
        #     TODO read it out from the input somehow
            raise ReadError('\'program\' not in calculator\'s parameters, not sure which energy to extract')
        #

        # TODO check for errors and warnings
        self.check_warnings_errors()

        # Check for SCF
        # TODO or pass something, not root?
        self.check_SCF_convergence()



        self.read_energy()
        # TODO think about how to deal with forces
        # TODO use 'self.parameters.task.fin('gradient') > -1 maybe?
        if self.parameters.task.find('forces') > -1:
            self.read_forces()

        self.read_other_stuff()



    @property
    def name(self):
        return self.__name__

    @property
    def input_fname(self):
        return self.label + '.inp'

    @property
    def output_fname(self):
        return self._label + '.out'

    def write_input(self, atoms, properties=None, system_changes=None):
        # make so that either read from .xml, .out, .inp or something else.
        # creates a directory
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        if 'program' not in p.keys():
            raise CalculatorSetupError('must give a name of the program')
        p.write(self.label + '.ase')
        ase.io.write(self.label+'.xyz', atoms) # TODO check why not self.atoms
        with open(self.input_fname, 'w') as f:
            if 'memory' in p.keys():
                # assume p.memory is given as a string 'amount, units'
                f.write(f'memory, {p.memory}\n')
            f.write('geomtyp=xyz\n')
            f.write(f'geom={self.label}.xyz\n')
            if 'basis' in p.keys():
                f.write(f'basis={p.basis}\n')
            f.write(f'{p.program}\n')

        pass

    # Optional:
    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
    #
    # def check_state(self):
    #     pass

    # def set_label(self):
    #     pass


    def read_parameters_from_input_file(self):
        pass


def get_molpro_command(self, molpro_command=''):
    # TODO documment this somewhere
    if molpro_command:
        return molpro_command
    elif 'MOLPRO_COMMAND' in os.environ:
        return os.environ['MOLPRO_COMMAND']
    else:
        return 'molpro'

def shell_stdouterr(raw_command, cwd=None):
    """Abstracts the standard call of the commandline, when
    we are only interested in the stdout and stderr
    """
    stdout, stderr = subprocess.Popen(raw_command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True,
                                      shell=True, cwd=cwd).communicate()
    return stdout.strip(), stderr.strip()



'''
from molpro2 import Molpro
from ase.build import molecule
at = molecule('CH4')
mp = Molpro(program='hf')
'''

