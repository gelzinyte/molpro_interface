#!/usr/bin/env python2
# vim: set encoding=utf-8
"""High-level drivers for electronic structure and related codes

Public abstract classes:
    ESCalculator        Interface to an electronic structure program
    ChargeExtractor     Get partial charges from an electronic structure

Public classes:
    MolproCalculator    Interface to the MOLPRO quantum chemistry code
    QTAIMCharges        Driver for the 'bader' QCT charge extractor

Public functions:
    compute_qm_single   Run a single electronic structure calculation
    compute_qm_traj     Run several electronic structure calculations
    compute_qm_traj_mproc
                        Run many calculations in parallel
"""

from __future__ import print_function
import abc
import errno
import logging
import os
import shutil
import subprocess
import tempfile

import numpy as np
import quippy
from quippy import Potential


logger = logging.getLogger(__name__)


def _quote_escape(name):
    """Escape names for use in FilePot arg strings"""
    if name is not None:
        return r'\"' + name + r'\"'
    else:
        return ''


class PersistentResult(object):

    """Interface to a result that can be retrieved at a later time"""

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def result_exists(self):
        """Check whether a result is available"""
        return False

    @abc.abstractmethod
    def get_result(self):
        """Get the result, assuming it is available

        It is up to the caller to make sure the result is available
        before calling this method.
        """
        raise NotImplementedError


class ESCalculator(PersistentResult):

    """General interface to an electronic structure program."""

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def calc(self, geom, work_dir):
        """Run the calculation on the given geometry

        Returns the result, which is also accessible using the
        get_result() method.  The return value of the latter contains
        the same data as the return value of the last invocation of this
        method.

        Any files created over the course of the calculation are placed
        in the given working directory.  The final result is stored in a
        more persistent location.
        """
        raise NotImplementedError


class ChargeExtractor(PersistentResult):

    """Interface to programs for getting partial charges from electronic
    structure calculation results."""

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def calc(self, density_file):
        """Extract the charges from the given density file

        Returns the result, which is also accessible using the
        get_result() method.  The return value of the latter contains
        the same data as the return value of the last invocation of this
        method.
        """
        raise NotImplementedError


class Workspace(object):

    """Interface to a directory to hold intermediate files.

    The object created should be used in a with-statement; this
    creates the scratch directory and deletes it when the statement
    exits.
    """

    def __init__(self, location):
        """Register a workspace in the directory given by 'location'"""
        self.location = location

    def __enter__(self):
        """Create the workspace directory

        Returns the name of the directory created"""
        self.work_dir = tempfile.mkdtemp(prefix='qm-',
                                         dir=self.location)
        return self.work_dir

    def __exit__(self, exc_type, exc_val, traceback):
        """Clean up the workspace (i.e. delete the directory)"""
        #print("TEST VERSION NOT REMOVING OUTPUTS")
        try:
            shutil.rmtree(self.work_dir)
        except OSError as ose:
            logger.warning("Got error while removing working directory %s",
                           ose.strerror)


class MolproCalculator(ESCalculator):

    """Class to run MOLPRO calculations

    Mainly wraps the MOLPRO FilePot interface.  Uses a context manager
    to manage the working directory; use a with-statement to set up the
    directory and gain access to a Potential object that will do the
    calculation.

    Public attributes:
        atoms_fname     Name of the file where the intermediate geometry
                        is written.  This also determines where the
                        calculation results are stored.
        work_dir        The directory where intermediate files for the
                        calculation are stored.  This directory is
                        deleted when this class's context manager exits.
        density_fname   Filename where the electron density is stored.
    """

    def __init__(self, template_fname, write_density=False, energy_from='RKS',
                 do_forces=True, nprocs=None, cmdl_args=''):
        """Set up a FilePot interface to MOLPRO

        Parameters:
            template_fname  Name of the MOLPRO template file to use
            write_density   Whether to write the electron density to a
                            file for later analysis (default False)
            energy_from     Keyword to use to extract the energy in the
                            MOLPRO output XML file (default 'RKS')
            do_forces       Whether the calculation will produce forces
                            that should be extracted (default True)
            nprocs          Number of processors to use for MOLPRO's
                            OpenMP parallelism (default 1, i.e. no
                            parallelism)
            cmdl_args       Additional arguments to be passed to molpro
                            on the command line

        If 'write_density' is True, the density will be written to a
        file named 'total_density.cube' in the working directory.
        """
        #TODO support passing in a calculation title?
        #TODO support putting atoms output file in working dir
        #TODO the cmdl_args string is a potential security loophole -
        #     sanitize first? Or just support a few known molpro args?
        pot = Potential('FilePot command=molpro_driver.py')
        pot.set(energy_from=energy_from)
        template_fname_abs = os.path.abspath(template_fname)
        pot.set(template=_quote_escape(template_fname_abs))
        pot.set(extract_forces=True)
        self.write_density = write_density
        self._pot = pot
        self._cur_atoms_fname = None
        self.density_fname = None
        self._noconverge_fname = None
        self.do_forces = do_forces
        self.cmdl_args = cmdl_args
        if nprocs is not None:
            self._set_nprocs(nprocs)

    def _set_nprocs(self, nprocs):
        """Set the number of processors for a parallel calculation"""
        if nprocs != 1:
            molpro_cmd = _quote_escape('molprop -n {} {:s} %s'.format(
                nprocs, self.cmdl_args))
        else:
            molpro_cmd = _quote_escape('molpros {:s} %s'.format(self.cmdl_args))
        self._pot.set(molpro=molpro_cmd)

    def set_atoms_fname(self, atoms_fname):
        self._cur_atoms_fname = atoms_fname
        atoms_fname_abs = os.path.splitext(os.path.abspath(atoms_fname))[0]
        # And _of course_, the atoms filename has to stay unescaped...
        self._pot.set(filename=atoms_fname_abs)
        self._out_fname = atoms_fname_abs + '.out'
        self._noconverge_fname = os.path.join(
            os.path.dirname(atoms_fname_abs), 'not-converged')

    def get_atoms_fname(self):
        """Filename where the intermediate geometry is stored."""
        return self._cur_atoms_fname

    atoms_fname = property(get_atoms_fname, set_atoms_fname)

    def calc(self, geom, work_dir, nprocs=None):
        """Do the calculation on the given geometry

        Other parameters:
            work_dir        Directory to hold temporary files used
                            during the calculation, also holds output
                            files that may be used by other calculations
            nprocs          Number of processors for a parallel
                            calculation

        Returns an Atoms object with forces and energies stored inside.

        Note: Modifies the input Atoms object; make a copy if you don't
        want this to happen.
        """
        if work_dir is None:
            work_dir = os.getcwd()
        self._pot.set(working_dir=_quote_escape(work_dir))
        if nprocs is not None:
            self._set_nprocs(nprocs)
        elif self._pot.get('molpro') is None:
            self._set_nprocs(1)
        atoms_fname_stem = os.path.splitext(
            os.path.basename(self.atoms_fname))[0]
        if self.write_density:
            self.density_fname = os.path.join(work_dir, atoms_fname_stem,
                                              'total_density.cube')
            cube_lines = 'cube,total.cube,-1 title,molpro'
            self._pot.set(append_lines=_quote_escape(cube_lines))
        molpro_out_xml = os.path.join(work_dir, atoms_fname_stem,
                                      atoms_fname_stem + '.xml')
        # Associate the working dir with the atoms_fname for crash diagnosis
        logger.info("Running MOLPRO on %s in working directory %s",
                    self.atoms_fname, work_dir)
        try:
            if self.do_forces:
                self._pot.calc(geom, energy=True, force=True)
            else:
                self._pot.calc(geom, energy=True)
        except RuntimeError as rte:
            # TODO Careful - this catch clause can also catch simple
            # programming errors. A way to distinguish non-converged
            # calculations from other screw-ups would be nice.
            # TODO TODO this is becoming exasperating. Need to look up
            # the error codes that molpro_driver.py actually spits out.
            logger.warning("Calculation did not finish on %s",
                           self.atoms_fname)
            logger.warning("Exception details:")
            logger.warning(rte)
            geom.energy = np.nan
            if self.do_forces:
                geom.force[:] = np.nan
            out_dir = os.path.dirname(self.atoms_fname)
            open(self._noconverge_fname, 'w').close()
            return geom
        finally:
            # Save the output files for later diagnosis
            shutil.copy2(molpro_out_xml,
                         os.path.splitext(self.atoms_fname)[0] + '.xml')
        # (Re)read the output file to get all the energies
        return self.get_result()

    def result_exists(self):
        """Determine whether a calculation has already been run

        Looks for the result associated with the currently stored
        'atoms_fname'
        """
        #TODO Check to make sure the file contains a valid Atoms?
        if os.path.isfile(self._noconverge_fname):
            # Don't want to redo the whole non-converging calculation
            #TODO would be better to have a separate method to check for
            #     non-convergence
            return True
        else:
            return os.path.isfile(self._out_fname)

    def get_result(self, fail_noconverge=False):
        """Retrieve this calculation's result as an Atoms object"""
        if os.path.isfile(self._noconverge_fname):
            logger.warning('Calculation on %s apparently failed',
                           self.atoms_fname)
            if fail_noconverge or not os.path.isfile(self._out_fname):
                atoms_out = quippy.Atoms(self.atoms_fname)
                atoms_out.energy = np.nan
                if self.do_forces:
                    atoms_out.properties['force'] = (np.nan *
                                                     np.ones((3, atoms_out.n)))
            else:
                logger.warning('Reading anyway.')
                atoms_out = quippy.Atoms(self._out_fname, format='xyz')
        else:
            #TODO support for rereading xml output file?
            atoms_out = quippy.Atoms(self._out_fname, format='xyz')
        return atoms_out


class QTAIMCharges(ChargeExtractor):

    """Extract Bader's topological/QTAIM charges from a density file

    Public attributes:
        atoms_fname     Filename specifying the geometry used for the
                        calculation
    """

    def __init__(self):
        """Construct a QTAIM charge calculator"""
        self._result_dir = None
        self._out_fname = None
        self._atoms_fname = None
        self._abs_atoms_fname = None
        self._noconverge_fname = None

    def set_atoms_fname(self, fname):
        self._atoms_fname = fname
        # Store an absolute filename because calc() changes directories
        self._abs_atoms_fname = os.path.abspath(fname)
        abs_directory = os.path.abspath(os.path.dirname(fname))
        self._result_dir = abs_directory
        self._out_fname = os.path.join(abs_directory, 'ACF.dat')
        self._noconverge_fname = os.path.join(
            os.path.dirname(self._abs_atoms_fname), 'not-converged')

    def get_atoms_fname(self):
        """Filename where the geometry for this calculation is stored"""
        return self._atoms_fname

    atoms_fname = property(get_atoms_fname, set_atoms_fname)

    def calc(self, density_fname):
        """Run the calculation and return the charges

        Parameters:
            density_fname   Where to find the density file (if relative,
                            must be relative to the directory containing
                            'atoms_fname')

        Result files are written into the same directory as
        'atoms_fname'
        """
        if os.path.isfile(self._noconverge_fname):
            return self.get_result()
        try:
            orig_dir = os.getcwd()
            os.chdir(self._result_dir)
            # The bader program doesn't return an error if the file
            # doesn't exist, so check explicitly
            with open(density_fname, 'r'):
                pass
            subprocess.check_call(['bader', density_fname])
            charges = self.get_result()
            return charges
        except subprocess.CalledProcessError as cpe:
            logger.warning("bader failed on file %s with return code %s",
                           density_fname, cpe.returncode)
        finally:
            os.chdir(orig_dir)

    def result_exists(self):
        """Determine whether a calculation has already been run

        Looks for the result associated with the currently stored
        'atoms_fname'
        """
        if not os.path.isfile(self._noconverge_fname):
            return os.path.isfile(self._out_fname)
        else:
            return True

    def get_result(self):
        """Return this calculation's result as an array of charges"""
        acf_charge_col = 4
        atoms_ref = quippy.io.AtomsReaders['xyz'](self._abs_atoms_fname)[0]
        if os.path.isfile(self._noconverge_fname):
            return np.nan * np.ones((atoms_ref.n, ))
        atoms_pcharge = np.genfromtxt(self._out_fname,
                                      skip_header=2, skip_footer=4)
        charges = atoms_ref.properties['Z'] - atoms_pcharge[:, acf_charge_col]
        # The np.array conversion is necessary because otherwise, the
        # result goes out of scope (grrr...)
        return np.array(charges)


def compute_qm_single(geom, el_calc, atoms_fname, wd_loc, pcharge_calc=None):
    """Run an electronic structure calculation on a single geometry

    Parameters:
        geom        Atoms object specifying the geometry
        el_calc     ESCalculator object to do the calculations
        atoms_fname Filename to hold the intermediate geometry
        wd_loc      Where to create the working directory
        pcharge_calc
                    ChargeExtractor object to use to extract partial
                    charges.  If None (the default), no partial charges
                    will be extracted.

    If the results of all calculations requested are already present,
    they are read and returned.  Otherwise, all calculations are run in
    full.

    Returns an Atoms object with the relevant quantities attached.
    """
    el_calc.atoms_fname = atoms_fname
    if pcharge_calc is not None:
        pcharge_calc.atoms_fname = atoms_fname
    if el_calc.result_exists():
        result = el_calc.get_result()
        if pcharge_calc is None:
            return result
        elif pcharge_calc.result_exists():
            charge = pcharge_calc.get_result()
            result.properties['partial_charge'] = charge
            return result
    with Workspace(wd_loc) as work_dir:
        result = el_calc.calc(geom, work_dir)
        if pcharge_calc is not None:
            charge = pcharge_calc.calc(el_calc.density_fname)
            result.properties['partial_charge'] = charge
    return result


def compute_qm_traj(traj, out_pref, el_calc, step=1,
                    write_xyz=True, pcharge_calc=None):
    u"""Run electronic structure calculations over a trajectory

    Parameters:
        traj        Iterable that sequentially returns Atoms objects in
                    the trajectory, can be e.g. a quippy.io.AtomsReader
                    object
        out_pref    Prefix for the output dirs (extension will be
                    stripped)
        el_calc     ESCalculator object to use for the calculations
        step        Run DFT only every (step)th snapshot of the
                    trajectory. Default 1, meaning use every step.
        write_xyz   Whether to write an XYZ file with the same name
                    as the main output directory (dirname(out_pref))
                    Default True
        pcharge_calc
                    ChargeExtractor object to use to extract partial
                    charges. If None (the default), no partial charges
                    are extracted.
    Not yet supported:
        vacuum_pad  Length by which to pad each dimension of the
                    periodic box (default 10 Å). Assumes an orthogonal
                    box.
        atom_id     Atom for which to return atom-local properties. If
                    None (the default), return properties for all atoms.

    Note: The vacuum_pad option is unnecesary for codes that don't do
    periodic calculations (e.g. MOLPRO).

    Returns a list of Atoms snapshots containing the desired properties.
    """
    atoms_snaps = []
    scratch_dir = os.path.expandvars('/scratch/$LOGNAME')
    if write_xyz:
        base_dir = os.path.dirname(out_pref)
        out_writer = quippy.AtomsWriter(base_dir + '.xyz')
    for step_idx, atoms in enumerate(traj):
        if (step_idx % step) != 0:
            continue
        output_dir_fmt = os.path.splitext(out_pref)[0] + "-step{:04d}"
        output_dir = output_dir_fmt.format(step_idx)
        try:
            os.makedirs(output_dir)
        except OSError as ose:
            if ose.errno != errno.EEXIST:
                raise ose
        atoms_fname = os.path.join(output_dir, "atoms-qm.xyz")
        result = compute_qm_single(atoms, el_calc, atoms_fname, scratch_dir,
                                   pcharge_calc)
        if write_xyz:
            out_writer.write(result)
        atoms_snaps.append(result)
    if write_xyz:
        out_writer.close()
    #TODO convert to AtomsList?
    return atoms_snaps


def _compute_qm_single_mproc(args):
    """Wrapper around compute_qm_single() for multiprocessing use

    Ports the tuple of args passed in from the multiprocessing map()
    invocation to the interface of compute_qm_single().  Signature is
    basically the same as compute_qm_single(), but see
    compute_qm_traj_mproc() for differences.

    TODO: It would have been nice to be able to just pass around
    calculator objects instead of constructing them here, but that
    causes scary memory errors - the fix for that will have to wait.
    """
    el_calc = MolproCalculator(*args[1])
    if args[4] is not None:
        pcharge_calc = QTAIMCharges(*args[4])
    else:
        pcharge_calc = None
    return compute_qm_single(args[0], el_calc, args[2], args[3],
                             pcharge_calc)


def compute_qm_traj_mproc(traj, out_pref, el_init_args, pool, step=1,
                          pcharge_init_args=None, map_chunksize=4):
    u"""Run electronic structure calculations over a trajectory

    Parallel asynchronous version using multiprocessing

    Parameters:
        traj        Iterable that sequentially returns Atoms objects in
                    the trajectory, can be e.g. a quippy.io.AtomsReader
                    object
        out_pref    Prefix for the output dirs (extension will be
                    stripped)
        el_init_args
                    Arguments to initialize the ESCalculator object for
                    the calculations
        pool        multiprocessing.Pool object to take the calculations
                    as independent processes
        step        Run DFT only every (step)th snapshot of the
                    trajectory. Default 1, meaning use every step.
        pcharge_init_args
                    Arguments to initialize the ChargeExtractor object
                    for extracting partial charges. If None (the
                    default), no partial charges are extracted.
        map_chunksize
                    Size of the chunks into which the atoms list is
                    split for iteration. Passed to the pool.map_async()
                    function as the chunksize argument.

    Returns an AsyncResult object that provides a list of the result
    Atoms objects when its get() method is called.
    """
    scratch_dir = os.path.expandvars('/scratch/$LOGNAME')
    def args_gen():
        for step_idx, atoms in enumerate(traj):
            if (step_idx % step) != 0:
                continue
            output_dir_fmt = os.path.splitext(out_pref)[0] + "-step{:04d}"
            output_dir = output_dir_fmt.format(step_idx)
            try:
                os.makedirs(output_dir)
            except OSError as ose:
                if ose.errno != errno.EEXIST:
                    raise ose
            atoms_fname = os.path.join(output_dir, "atoms-qm.xyz")
            yield (atoms, el_init_args, atoms_fname, scratch_dir,
                   pcharge_init_args)
    return pool.map_async(_compute_qm_single_mproc, args_gen(),
                          chunksize=map_chunksize)

