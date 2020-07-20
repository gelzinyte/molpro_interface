"""The module defines an interface to MOLPRO 2012 (?).
    Based on Python 2 molpro driver for quippy by Alan Nichol, James Kermode
    and others (????). Adapted to ASE by Elena Gelžinytė.
    Loosely based on castep ase calculator.
"""

# TODO list:
#   Raise appropriate calculator errors
#   Clean up imports to follow python conventions

import os
import pdb
import ase
from ase.calculators.calculator import FileIOCalculator, Parameters, \
                                       CalculatorSetupError, Calculator

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
    #     program='hf')

    def __init__(self, restart=None, ignore_bad_restart_file=False, label='molpro',
               atoms=None, molpro_command=None, directory='MOLPRO', **kwargs):

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
                atoms=atoms, directory=directory, **kwargs)

        self.molpro_command = self.get_molpro_command(molpro_command)






    def calculate(self, atoms=None, properties=['energy'], systm_changes=\
                  ['positions', 'numbers', 'cell', 'pbc', 'initial_charges',
                   'initial_magmoms']):
        pass

    def read(self):
        pass

    def read_results(self):
        pass

    def write_input(self, atoms, properties=None, system_changes=None):
        # creates a directory
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        if 'program' not in p.keys():
            raise CalculatorSetupError('must give a name of the program')
        p.write(self.label + '.ase')
        # write xyz file to be included by molpro.
        # TODO workout how to preserve all at.info entries
        # alternative is to write .xyz directly to template, but I think this
        # is more convenient to inspect (e.g. with .xyz viewer). Maybe add
        # an option?
        ase.io.write(self.label+'.xyz', atoms) # TODO check why not self.atoms
        with open(self.label+'.inp', 'w') as f:
            if 'memory' in p.keys():
                # assume p.memory is given as a string 'amount, units'
                f.write(f'memory, {p.memory}')
            f.write('geomtyp=xyz\n')
            f.write(f'geom={self.label}.xyz')
            if 'basis' in p.keys():
                f.write(f'basis={p.basis}')
            f.write(p.program) # TODO add check and raise error if command is not given



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

    def get_molpro_command(molpro_command=''):
        # TODO documment this somewhere
        if molpro_command:
            return molpro_command
        elif 'MOLPRO_COMMAND' in os.environ:
            return os.environ['MOLPRO_COMMAND']
        else:
            return 'molpro'
