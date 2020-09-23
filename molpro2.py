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
from ase.calculators.calculator import FileIOCalculator, \
                                       CalculatorSetupError, Calculator, \
                                       ReadError, PropertyNotImplementedError,\
                                       SCFError, PropertyNotPresent

from ase.calculators.calculator import Parameters

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

    implemented_properties = ['energy', 'forces']

    default_parameters = dict(task='gradient',
                              basis='6-31G*',
                              functional='b3lyp')


    def __init__(self, restart=None, ignore_bad_restart_file=False, label='MOLPRO/molpro',
               atoms=None, molpro_command=None, **kwargs):

        # TODO add an option to read out a template file to read out all the parameters
        # TODO function to make a template file out of everything
        # TODO should I input and link .xyz file or coordinates directly?
        # TODO deal with the directories - scratch (what's the default?), working
        # dir to write in and out all the tmp files, an option to not delete that
        # and how to deal with duplicate files that molpro just moves over
        # Does creating and removing directories with every calculation slow it down?
        # Maybe add a cleanup function that once everything is done removes all the tmp dirs
        # TODO add a function to use the calculator parallely.
        # TODO have a problem if .xyz file that is passed in has a number of at.info stuff.
        # Set it up so that the molpro input has no info, but when it's read out, at.info from original at file is returned.
        # actually probably taken care of by ase itself.
        # TODO when to make functions part of Molpro class and when just include them in this .py file?

        self.__name__ = 'Molpro'
        self.molpro_command = get_molpro_command(molpro_command)

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

        # TODO sort out how command=molpro_command works out in FileIOCalculator
        FileIOCalculator.__init__(self, restart=restart,
                ignore_bad_restart_file=ignore_bad_restart_file, label=label,
                atoms=atoms, command=molpro_command, **kwargs)


        self._no_scf_iterations = None
        self._prepare_input_only = False
        self._rename.existing_dir = True
        self._scf_maxit = 60
        self._geomtyp = 'xyz'


        #######################
        # scrape this and put tests somewhere else
        ###########################

        # TODO get default parameters here maybe and then overwrite?
        # TODO document these in the docummenttion
        # TODO how to enforce compulsory arguments? Or have defaults? 
        # set keyword arguments
        for keyword, item in kwargs.items():
            # TODO maybe make everything be uppercase so capitalisation miscommunication is avoided??
            # set and check for specific allowed arguments related to molpro
            if keyword == 'task':
                supported_tasks = ['gradient', 'energy_only', 'optg']
                if item not in supported_tasks:
                    raise CalculatorSetupError(
                        f'{item} unrecognised, {keyword} must be one of {supported_tasks}')
                self.task = item
                
            elif keyword == 'basis':
                # if keyword not in []  # TODO scrape a list out of the basis set library?
                self.basis = item
                
            elif keyword == 'command':
                supported_commands = ['RKS', 'UKS', 'HF'] # TODO extend/scrape/transfer from previous calculator
                if item not in supported_commands:
                    raise CalculatorSetupError(f'{item} unrecognised, {keyword} must be one of {supported_commands}')
                self.command = item
                
            elif keyword == 'functional':
                supported_functionals = ['b-lyp', 'pbe', 'b3lyp']
                if item not in supported_functionals: # TODO extend/scrape/transfer
                    raise CalculatorSetupError(f'{item} unrecognised, {keyword} must be one of {supported_functionals}')

            elif keyword == 'extended_command':
                if self.command not in item:
                    # TODO how to know what is the appropriate error to raise?
                    raise Exception(f"Keyword 'command'={self.command} is not found in the {keyword}")
                self.extended_commnad = item




    def calculate(self, atoms=None, properties=['energy'], systm_changes=\
                  ['positions', 'numbers', 'cell', 'pbc', 'initial_charges',
                   'initial_magmoms']):
        # TODO check what to do with properties and system_changes
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self, atoms, properties, system_changes)
        if not self._prepare_input_only:
            self.run()
            self.read()
            self.push_old_state()

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
        pass

    def check_SCF_convergence(self):
        # and read other stuff too
        program = self.parameters.program
        root = et.parse(slef.label + '.xml').getroot()
        job = root.find('{http://www.molpro.net/schema/molpro-output}job')
        jobsteps = job.findall(
            '{http://www.molpro.net/schema/molpro-output}jobstep')
        for jobstep in jobsteps:
            if program.upper() in jobstep.attrib['command']:
                for child in jobstep:
                    text = child.text
                    if text is not None and 'ITERATION' in text:
                        text_lines = text.splitlines()
                        for idx, line in enumerate(text_lines):
                            if 'MAX. NUMBER OF ITERATIONS' in line:
                                maxit = int(re.search(r'\d+', line).group())
                                self._scf_maxit = maxit
                            if 'Final occupancy' in line:
                                last_scf_line = text_lines[idx - 2].strip()
                                # ^ and search are probably doubling up
                                match = re.search(r'^\d+',
                                                  last_scf_line).group()
                                last_iteration = int(match)
                                if last_iteration == self._scf_maxit:
                                    raise SCFError('Maxed out of SCF iterations')
                                elif last_iteration > self._scf_maxit:
                                    raise RuntimeError('last_iteration > self._scf_maxit, something is worng with the code')

                            # TODO Spot messages for not converging and
                            #  stopping after a few iterations
                        # only interested in this (TODO: unless more than one program), so can return (?).
                        return

    def read_other_stuff(self):
        pass





    @property
    def name(self):
        return self.__name__

    @property
    def input_fname(self):
        return self.label + '.inp'

    @property
    def output_fname(self):
        return self._label + '.out'



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

