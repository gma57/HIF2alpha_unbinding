#!/usr/bin/env python
# 
# wcrawl_functions.py
#
# Written 15.09.29; updated 18.04.10 by Alex DeGrave
from __future__ import print_function, division; __metaclass__ = type
import numpy
from scipy.spatial.distance import pdist
import h5py
import os
import sys
from westpa.core import h5io
from westpa.cli.tools.w_crawl import WESTPACrawler
import mdtraj

class IterationProcessor(object):
    '''
    This class performs analysis on each iteration.  It should contain a method
    ``process_iteration``, which may be called as 
    ``process_iteration(self, n_iter, iter_group)``, where ``n_iter`` refers to
    the weighted ensemble iteration index, and ``iter_group`` is the HDF5 group 
    for the given iteration. The method should return an array or other values, 
    which the ``process_iter_result`` method of the ``Crawler`` class recieves 
    as the argument ``result``. 
    '''
    # Store the location of the PDB file to be used as the topology
    #topology_filename = '/ocean/projects/mcb180038p/jml230/bdpa_wsh2045_p2_3/common_files/system.prmtop'
    #strip_topology_filename = '/ocean/projects/mcb180038p/jml230/bdpa_wsh2045_p3/common_files/bdpa_wsh2045_nowat.prmtop'
    #ionless_topology_filename = '/ocean/projects/mcb180038p/jml230/bdpa_wsh2045_p3/common_files/bdpa_wsh2045_noions.prmtop'
    # Define the pattern used for finding each segment's traj file
    #iter_pattern = 'iterations/iter_{n_iter:08d}'
    traj_pattern = '/ix/lchong/ris130/WESTPA/new_simulations_unb/unb1_MULTIMAB_rates/traj_segs/{n_iter:06d}/{seg_id:06d}'
    # In this example, there are three frames saved for each segment.
    nframes = 11

    def __init__(self):
        '''
        Initialize the IterationProcessor class
        '''
    def process_iteration(self, n_iter, iter_group):
        '''
        The main analysis function that w_crawl calls for each iteration.
        This should be changed based on your analysis. This method could
        contain all the code for your analysis, or it could call an outside
        function. 

        ----------
        Parameters
        ----------
        n_iter: (int) The index of the weighted ensemble iteration for which
          analysis should be performed.
        iter_group: (H5py group) The hdf5 group corresponding to iteration
          n_iter, from the the main WESTPA data file (typically west.h5)

        -------
        Returns
        -------
        result: (numpy.ndarray) In general this could be an object, which is
          later processed by Crawler.process_iter_result. Here, it is an array
          of the center of mass of the protein. The array has shape 
          (n_segments, n_timepoints, 3), where dimension 0 indexes the segment, 
          dimension 1 indexes the frame number, and dimension 2 indexes the 
          x/y/z coordinate of the center of mass.
        '''
        # Find the number of segments in the iteration at hand
        #print("starting",n_iter)
        nsegs=iter_group['seg_index'].shape[0]
        # parent_iter = n_iter-1
        
        # The dimensionality of the data you wish to store
        #data_dims = 1
        
        # Create an array to hold your data
        #iteration_data_array = []
        #for i in range(0,6):
        #    iteration_data_array.append(numpy.zeros((nsegs, self.nframes)))
            
        #iteration_data_array=numpy.zeros((nsegs, self.nframes))

        #iteration_data_array.append(numpy.zeros((nsegs, self.nframes, 923, 3)))
        data_dims=6
        iteration_data_array = numpy.zeros((nsegs, self.nframes, data_dims))

        
        #print(iteration_data_array.shape)
        home = '/ix/lchong/ris130/WESTPA/new_simulations_unb/unb1_MULTIMAB_rates'
        surf_file = home + "/crawl_keydistances/cpptraj_dist.in"
	#wat_prmtop = home + "/common_files/bound_new.prmtop"
        #nowat_prmtop = home + "/common_files/bdpa_wsh2045_nowat.prmtop"
        #noions_prmtop = home + "/common_files/bdpa_wsh2045_noions.prmtop"
        #ref_folded = home + "/reference/wsh2045_eq3.pdb"
        #ref_avg = home + "/reference/wsh2045_avg.pdb"
        #ref_helix = home + "/reference/wsh2045_avg_helix.pdb"
        
        #refPDBfile = '/ocean/projects/mcb180038p/jml230/bdpa_wsh2045_p3_r7_1/reference/wsh2045_eq3.pdb'
        #ref_topology = '/ocean/projects/mcb180038p/jml230/bdpa_wsh2045_p3_r7_1/common_files/bdpa_wsh2045_nowat.prmtop'
        

        #atom_slice = numpy.asarray(range(0,923))

        #parentTraj = 'parent.ncrst'
        #childTraj = 'seg.nc'


        # Iterate over each segment
        for iseg in range(nsegs):
            #print(iseg)

            pathy = self.traj_pattern.format(n_iter=n_iter, seg_id=iseg)
            #print(n_ter,iseg)
            os.chdir(pathy)
            
            try:
                os.symlink(surf_file, pathy+'/cpptraj_dist.in')
            except FileExistsError:
                os.unlink(pathy+'/cpptraj_dist.in')
                os.symlink(surf_file, pathy+'/cpptraj_dist.in')
            #try:
            #    os.symlink(wat_prmtop, pathy+'/system.prmtop')
            #except FileExistsError:
            #    pass
            #try:
            #    os.symlink(nowat_prmtop, pathy+'/bdpa_wsh2045_nowat.prmtop')
            #except FileExistsError:
            #    pass
            #try:
            #    os.symlink(ref_folded, pathy+'/wsh2045_eq3.pdb')
            #except FileExistsError:
            #    os.unlink(pathy+'/wsh2045_eq3.pdb')
            #    os.symlink(ref_folded, pathy+'/wsh2045_eq3.pdb')
            #try:
            #    os.symlink(noions_prmtop, pathy+'/bdpa_wsh2045_noions.prmtop')
            #except FileExistsError:
            #    pass

            command = 'cpptraj -i cpptraj_dist.in > /dev/null '
            #command = 'cpptraj -i percent_calc.in > /dev/null'
            os.system(command)

            


            x1,x2,x3,x4,x5,x6 = numpy.loadtxt(f"{pathy}/dist.dat", skiprows=1, usecols=(3,6,9,12,15,18), unpack=True)
            #x1 is complex, x2 is protein, x3 is ligand
            #ene_int = x1 - (x2+x3)
            arr1 = numpy.column_stack(numpy.array([x1]))
            arr2 = numpy.column_stack(numpy.array([x2]))
            arr3 = numpy.column_stack(numpy.array([x3]))
            arr4 = numpy.column_stack(numpy.array([x4]))
            arr5 = numpy.column_stack(numpy.array([x5]))
            arr6 = numpy.column_stack(numpy.array([x6]))

            #arrN = numpy.column_stack(numpy.array([-ene_int]))


            #y1,y2,y3 = numpy.loadtxt(f"{pathy}/new_contacts_intra5.dat", skiprows=1, usecols=(1,3,5), unpack=True)
            #c = numpy.loadtxt(f"{pathy}/drms-avg.dat", skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10))
            #d = numpy.loadtxt(f"{pathy}/drms-helix.dat", skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10))

            #b = x1 + x2 + x3
            #c = y1 + y2 + y3
            #e = b + c
            iteration_data_array[iseg, :, 0] = arr1[:,0]
            iteration_data_array[iseg, :, 1] = arr2[:,0] 
            iteration_data_array[iseg, :, 2] = arr3[:,0] 
            iteration_data_array[iseg, :, 3] = arr4[:,0] 
            iteration_data_array[iseg, :, 4] = arr5[:,0] 
            iteration_data_array[iseg, :, 5] = arr6[:,0] 

            #iteration_data_array[2][iseg, :] = x3
            #iteration_data_array[3][iseg, :] = y1
            #iteration_data_array[4][iseg, :] = y2
            #iteration_data_array[5][iseg, :] = y3
            #iteration_data_array[iseg, :, 1] = b[1]

            #ref = mdtraj.load(refPDBfile)
            #ref_slice = ref.topology.select('resid 5 to 16 or resid 23 to 35 or resid 40 to 53')
            #atom_slice = ref.topology.select('resid 0 to 58')
            #ref = ref.atom_slice(atom_slice)
            #try:
            #    struct0 = mdtraj.load(f'{parentTraj}', top=refPDBfile, atom_indices=atom_slice)
            #except OSError:
            #    log.warning("Parent traj file doesn't exist, loading reference structure coords")
            #    struct0 = numpy.squeeze(model.reference_structure._xyz)
    
            #struct1 = mdtraj.load_netcdf(f'{childTraj}', top=ref_topology, atom_indices=atom_slice)
 
            #struct0 = struct0.superpose(ref, atom_indices=ref_slice)
            #struct1 = struct1.superpose(ref, atom_indices=ref_slice)
            
            #all_coords = numpy.squeeze(numpy.concatenate((struct0._xyz, struct1._xyz)))
            
            #iteration_data_array[3][:, :, :] = all_coords

            #for num, val in enumerate(b):
            #    iteration_data_array[iseg, num, 0] = b[num]
            #    iteration_data_array[iseg, num, 1] = c[num]

            #os.unlink(pathy + '/dmatrix_calc.in')
            #os.unlink(pathy + '/system.prmtop')
            #os.unlink(pathy + '/bdpa_wsh2045_nowat.prmtop')
            #os.unlink(pathy + '/wsh2045_eq3.pdb')
            #os.unlink(pathy + '/bdpa_wsh2045_noions.prmtop')
           
            #with open('SASA_calc.dat', 'w') as f_open:
            #    f_open.write("frame\tHC_SASA\tAMIDE_SASA\tm-value\tdeltaCP\n")
            #    for idx in range(0,b.shape):
            #        d = 
            #        f_open.write(print(a[idx]+"\t"+b[idx]+))
            
 
            with open('/ix/lchong/ris130/WESTPA/new_simulations_unb/unb1_MULTIMAB_rates/crawl_keydistances/update.txt', 'a') as fo:
                fo.write("finished iteration "+str(n_iter)+" and segment "+str(iseg))
                fo.write('\n')

            os.chdir(home)

        return iteration_data_array

class Crawler(WESTPACrawler):
    '''
    In this example, w_crawl works as follows:

    We supply the ``Crawler`` class, which handles writing data. The 
    Crawler specifies 3 methods: initialize, finalize, and process_iter_result.

    ``initialize`` is called only once--when w_crawl starts up. The job of 
    initialize is to create the output file (and HDF5 file).

    Like ``initialize``, ``finalize`` is also called only once--when w_crawl
    finishes calculations for all iterations. The job of ``finalize`` is to
    gracefully close the output file, preventing data corruption.

    The method ``process_iter_result`` is called once per weighted ensemble
    iteration. It takes the weighted ensemble iteration (n_iter) and the result
    of the calculations for an iteration (result) as arguments, and stores the
    data in the output file.

    The actual calculations are performed by the IterationProcessor class 
    defined above. In particular, the IterationProcessor.process_iteration 
    method performs the calculations; the return value of this method is passed
    to Crawler.process_iter_result.
    '''

    def initialize(self, iter_start, iter_stop):
        '''
        Create an HDF5 file for saving the data.  Change the file path to
        a location that is available to you. 
        '''
        self.output_file = h5io.WESTPAH5File('./crawl.h5', 'w')
        h5io.stamp_iter_range(self.output_file, iter_start, iter_stop)
        pass

    def finalize(self):
        self.output_file.close()
        pass

    def process_iter_result(self, n_iter, result):
        '''
        Save the result of the calculation in the output file.

        ----------
        Parameters
        ----------
        n_iter: (int) The index of the weighted ensemble iteration to which
          the data in ``result`` corresponds.
        result: (numpy.ndarray) In general this could be an arbitrary object
          returned by IterationProcessor.process_iteration; here it is a numpy
          array of the center of geometry.
        '''
        pass
        # Initialize/create the group for the specific iteration
        iter_group = self.output_file.require_iter_group(n_iter)

        iteration_data_array = result
        
        # Save datasets
        dataset = iter_group.create_dataset('dist_lig', 
                                            data=iteration_data_array, 
                                            scaleoffset=6, 
                                            compression=4,
                                            chunks=h5io.calc_chunksize(
                                                    iteration_data_array.shape,
                                                    iteration_data_array.dtype
                                                                       )
                                            )
#        dataset1 = iter_group.create_dataset('iter_energy_N', 
#                                            data=iteration_data_array[1], 
#                                            scaleoffset=6, 
#                                            compression=4,
#                                            chunks=h5io.calc_chunksize(
#                                                    iteration_data_array[1].shape,
#                                                    iteration_data_array[1].dtype
#                                                                       )
#                                            )
#        dataset2 = iter_group.create_dataset('H1_H3', 
#                                            data=iteration_data_array[2], 
#                                            scaleoffset=6, 
#                                            compression=4,
#                                            chunks=h5io.calc_chunksize(
#                                                    iteration_data_array[2].shape,
#                                                    iteration_data_array[2].dtype
#                                                                       )
#                                            )
#        dataset = iter_group.create_dataset('H1_intra', 
#                                            data=iteration_data_array[3], 
#                                            scaleoffset=6, 
#                                            compression=4,
#                                            chunks=h5io.calc_chunksize(
#                                                    iteration_data_array[3].shape,
#                                                    iteration_data_array[3].dtype
#                                                                       )
#                                            )
#        dataset1 = iter_group.create_dataset('H2_intra', 
#                                            data=iteration_data_array[4], 
#                                            scaleoffset=6, 
#                                            compression=4,
#                                            chunks=h5io.calc_chunksize(
#                                                    iteration_data_array[4].shape,
#                                                    iteration_data_array[4].dtype
#                                                                       )
#                                            )
#        dataset2 = iter_group.create_dataset('H3_intra', 
#                                            data=iteration_data_array[5], 
#                                            scaleoffset=6, 
#                                            compression=4,
#                                            chunks=h5io.calc_chunksize(
#                                                    iteration_data_array[5].shape,
#                                                    iteration_data_array[5].dtype
#                                                                       )
#                                            )


# Entry point for w_crawl
iteration_processor = IterationProcessor()
def calculate(n_iter, iter_group):
    '''Picklable shim for iteration_processor.process_iteration()'''
    global iteration_processor 
    return iteration_processor.process_iteration(n_iter, iter_group)

crawler = Crawler()
