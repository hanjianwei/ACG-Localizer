ACG Localizer Version 1.2.1
www.rwth-graphics.de/localization/
tsattler@cs.rwth-aachen.de

------------
What is this
------------
ACG Localizer is the C++ implementation of the localization frameworks proposed
in the papers 
  
  T. Sattler, B. Leibe, L. Kobbelt. Fast Image-Based Localization using Direct
2D-to-3D Matching. ICCV 2011.

  T. Sattler, B. Leibe, L. Kobbelt. Improving Image-Based Localization by Active
Correspondence Search. ECCV 2011.

This is a short readme file to help you compiling and using ACG Localizer.
It also states the conditions under which you might use the software.

------------------
About the software
------------------
This software is intended to give researchers a reference implementation of the
localization methods presented in the papers

  T. Sattler, B. Leibe, L. Kobbelt. Fast Image-Based Localization using Direct
2D-to-3D Matching. ICCV 2011.

  T. Sattler, B. Leibe, L. Kobbelt. Improving Image-Based Localization by Active
Correspondence Search. ECCV 2011.

It is by no means a polished product and neither is it intended to be, although
we are planning to release a more user-friendly version of the code by the end
of 2012. The program has been developed and tested under Debian Linux with the
GCC compiler.
Please understand that we do not plan to port it to Windows or Mac OS X and that
we will not offer support for porting the software.


-------------------------------
License and conditions of usage
-------------------------------
The ACG Localizer software (except for the 3rd party libraries which have their
own licenses) is released under the GNU General Public License V3 (see COPYING).
For commercial use of the software, please contact Torsten Sattler
(tsattler@cs.rwth-aachen.de).

Any research or publication using this software must explicitly acknowledge
using this software. In case of a publication, citing the corresponding paper is
sufficient. If you are using acg_localizer, then you should cite

Torsten Sattler, Bastian Leibe and Leif Kobbelt.  Fast Image-Based Localization
using Direct 2D-to-3D Matching.
 13th IEEE International Conference on Computer Vision (ICCV 2011). 2011.
 
If you are using acg_localizer_active_search, then you should cite 

Torsten Sattler, Bastian Leibe and Leif Kobbelt.  Improving Image-Based
Localization by Active
Correspondence Search.
 12th European Conference on Computer Vision (ECCV 2012). 2012.
 
If you are just using the RANSAC functionality provided in this software
release, you should cite the appropriate papers (see RANSAC.hh).


------------
Requirements
------------
ACG Localizer uses the following libraries. If not specified otherwise, you have
to install them by yourself.
* ANN. 
   The Approximate Nearest Neighbor Library written by David M. Mount and Sunil
   Arya, available at http://www.cs.umd.edu/~mount/ANN/ . ACG Localizer has been
   tested with version 1.1.2. Also, we recommend to set the data type of
   ANNcoord to float to reduce memory requirements. 
* FLANN.
   The Fast Library for Approximate Nearest Neighbors written Marius Muja and
   David G. Lowe, available at
   http://people.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN . ACG Localizer has
   been tested with version 1.6.11. Due to an API breaking change in version
   1.7.0, we suggest that you use 1.6.11, which (as of 02/29/2012) can still be
   found at http://people.cs.ubc.ca/~mariusm/uploads/FLANN/flann-1.6.11-src.zip
   We plan to update our source code to be compliant with newer versions of
   FLANN.   
* GMM
   Available at http://download.gna.org/getfem/html/homepage/gmm/index.html .
* LAPACK
* OpenMesh
   The vector class of OpenMesh is used for representing 3D points. OpenMesh is
   available at http://openmesh.org . However, since only the vector class is
   required, the required files from OpenMesh are already included in the source
   files of ACG Localizer.
* jhead
   A tool to read exif tags from jpg files, available at
   http://www.sentex.net/~mwandel/jhead/index.html . It is already included in
   this release of ACG Localizer.
* SFMT
  The SIMD oriented Fast Mersenne Twister implementation by Mutsuo Saito and 
  Makoto Matsumoto. The software is used to generate random numbers for RANSAC
  sampling. It is also already included in this release of ACG Localizer.

We use CMake for compiling the project. Please make sure that you have CMake
installed (http://www.cmake.org/).


------------
Installation
------------

IMPORTANT: Before compiling the ACG Localizer, please copy the modified version
of hkmeans_index.h which can be found in the directory flann_modification/ into
the FLANN_DIR/src/cpp/flann/algorithms directory and replace the original
hkmeans_index.h file. This modified file contains some additional functions
needed to read out the coarser vocabularies defined by different levels in a
vocabulary tree. After copying the file, you have to re-compile FLANN.

Assuming that you unpacked ACG Localizer into the directory
SOME_DIRECTORY/ACG_Localizer, you first have to edit the CMake files that are
responsible for finding the required libraries. These files are located in
SOME_DIRECTORY/ACG_Localizer/cmake. You have to edit the FindANN.cmake and
FindFLANN.cmake files to fill in the directories where the include and library
files of both libraries are located. You probably don't have to edit
FindGMM.cmake and FindLAPACK.cmake while FindOpenMesh.cmake does not have to be
edited. After editing the files, you can start building ACG Localizer using
cmake:
* cd SOME_DIRECTORY/ACG_Localizer
* mkdir build
* cd build
* cmake ..
If no library is missing and if CMake threw no errors you can build the software
by simply typing
* make

To get the best performance, make sure that you set the build type to "Release".

Typing 'make install' will install the generated executables in
SOME_DIRECTORY/ACG_Localizer/build/bin . Notice that to achieve timing results
similar to those in the paper, you have to
set the build type of CMake to "Release". You can do this as follows (assuming
that you have installed the gui-version of cmake):
* cd SOME_DIRECTORY/ACG_Localizer/build
* cmake-gui . (or "ccmake ." when you want to use the command line version of
cmake)
In the field "CMAKE_BUILD_TYPE" type "Release" (without quotation marks), then
click "Configure" and "Generate". Afterwards, compile the software with make.


-----
Usage
-----
During the compilation, four executables are generated:
* Bundle2Info
* compute_desc_assignments
* acg_localizer
* acg_localizer_knn
* acg_localizer_active_search

The last three executables are the actual localization methods. One for the
vocabulary-based prioritized search proposed in the ICCV 2012 paper
(acg_localizer), one for for kd-tree based localization (acg_localizer_knn) as
also described in the ICCV paper, and one using the active correspondence search
mechanism described in the ECCV 2012 paper. The first two executables can be
used to extract the information needed by the localization methods from a
structure-from-motion reconstruction computed by Bundler
(http://phototour.cs.washington.edu/bundler/). If you are using any other
software to compute structure-from-motion reconstructions you will have to write
your own converter that converts to the file format needed by the localization
methods.

Suppose that you have a structure-from-motion reconstruction computed by
Bundler, located in the directory BUNDLE_DATA, where the output generated by
Bundler can be found at BUNDLE_DATA/bundle/bundle.out and the list of used jpg
images is stored at BUNDLE_DATA/list.txt. Then you can generate the data needed
by acg_localizer as follows:
* cd BUNDLE_DATA

* Bundle2Info bundle/bundle.out list.txt bundle.info
  This will create a binary file called bundle.info which contains information
about the number of cameras and 3D points contained in the reconstruction as
well as information about the 3D points and their SIFT descriptors. Bundle2Info
will try to open the keyfiles containing the SIFT descriptors. It assumes that
the keypoint files are text files that have the same format as the ones
generated by David Lowe's SIFT binary (http://www.cs.ubc.ca/~lowe/keypoints/).
Bundle2Info further assumes that if list.txt contains a jpg image called a.jpg
that the corresponding keypoint file is called a.key. Notice that Bundler uses
gzip to compress the keypoint files, so you need to decompress them before
running Bundle2Info. Run Bundle2Info without any parameters for information
about its parameters.
  
* compute_desc_assignments bundle.info 1 100000 clusters.txt bundle.desc_assignments.integer_mean.kdtree.clusters.100k.bin 6 0 0
  The executable compute_desc_assignments will create a binary file called
bundle.desc_assignments.integer_mean.kdtree.clusters.100k.bin that contains all
information needed by acg_localizer. 
  The first parameter is the binary file we have just created. 
  The second parameter specifies the number of trees in a random kd-tree forest
that should be used to assign each descriptor of every 3D point in the
reconstruction to a visual word (cluster centers obtained from k-means). 
An example cluster file containing 100k visual words can be found at
www.rwth-graphics.de/localization/.
Note that this cluster file is compatible with SiftGPU
(http://cs.unc.edu/~ccwu/siftgpu/), but might not work well with other SIFT
implementations due to a different ordering of the descriptors.
  The next parameter specifies the number of visual words / cluster centers,
while "clusters.txt" is a text file containing those cluster centers. Any
clustering computed by k-means can be used as long as you make sure that you
compute the clustering on unnormalized SIFT descriptors similar to those
contained in the keypoint files. The computed cluster centers can then be
written into a text file, where the entries are separated by blanks and every
line contains one cluster center. 
  The second to last parameter describes the type of 3D point representation
that should be used. For example, 5 is the "all descriptor" representation and 6
is the "integer mean per vw" representation (cf. the paper). If you want to use
the kd-tree localization method (acg_localizer_knn) you should use 2 or 3.
  The second to last parameter defines the search structure that is used to
compute the assignments of the descriptors of the 3D points to the visual words.
For using the ICCV localization method (acg_localizer), this parameter has to be
set to 0 to use a single kd-tree. When the resulting assignments file should be
used with the ECCV method (acg_localizer_active_search), the parameter has to be
set to 1 to use a vocabulary tree. It is crucial that the search structure that
is used to compute the assignments of the points is the same that is used in the
localization method!
The last parameter defines the file type of the bundle.info file. All .info
files generated by Bundle2Info have type 0.
But the aachen.info file that is part of the Aachen dataset from the paper 

Torsten Sattler, Tobias Weyand, Bastian Leibe, Leif Kobbelt. Image Retrieval for
Image-Based Localization Revisited. BMVC 2012.

has a slightly different format, so you have to set the last parameter to 1. If
you are interested in obtaining the dataset, please contact Torsten Sattler
(tsattler@cs.rwth-aachen.de).

  Run compute_desc_assignments without parameters for a description of the
program.
  
You can now run the localization methods, for example by typing

* acg_localizer list.query.keys.txt 1 100000 clusters.txt bundle.desc_assignments.integer_mean.kdtree.clusters.100k.bin 0 0.2 100 results.txt

The first parameter specifies a text file containing the filenames of all
keypoint files that should be used as query images. For example, if you want to
register the jpg files a.jpg and b.jpg, then list.query.keys.txt should look
like this:
a.key
b.key
  acg_localizer assumes that given the filename of a keypoint file, the
corresponding jpg file can be found by replacing the .key ending with .jpg.
Notice that the jpg files have to contain the width and height of the image in
their exif tags as it is required for the pose estimation.
  The next three parameters specify the number of trees to use, the number of
cluster centers / visual words and the text file containing the centers,
similarly to compute_desc_assignments. We strongly advice you to use the same
settings as for compute_desc_assignments since we do a very approximate nearest
neighbor search when assigning features in the query images to visual words (a
kd-tree based search visiting at most 10 leaf nodes). If you use different
settings, it is very unlikely to get similar assignments as those from
compute_desc_assignments, which will lead to poor performance. 
  The fifth parameter specifies the filename of an assignment file generated by
compute_desc_assignments. 
  The sixth parameter tells the method how descriptors are stored in that file.
If it is set to 0, descriptors are assumed to be saved as unsigned char values,
while setting it to 1 assumes that descriptors are stored using floating point
values. If you want to use the "all descriptor" representation you have to set
it to 2 to signal that multiple descriptors stored in one visual word might
belong to the same 3D point. 
  The next parameter is the assumed inlier ratio for the RANSAC-based pose
estimation (parameter R from the paper). Setting it to 0.2 forces RANSAC to
assume an inlier ratio of 20%, i.e., RANSAC will take at most ceil( log( 0.05) /
log( 1 - 0.2^6 ) ) samples (RANSAC stops if the probability of missing the best
hypothesis is below 5%). The recommended value for this parameter is 0.2 .
  The next parameter corresponds to the parameter N_t from the paper and signals
that the prioritized search should be stopped after finding 100
correspondences. 
The recommended value for this parameter is 100.
  The last parameter is the filename of a text file in which acg_localizer will
dump information about the localization process. For every query image it
contains one line containing the six following numbers:
1. The number of inliers found by RANSAC.
2. The number of correspondences found using prioritized search.
3. The time needed to assign all 2D features in the query image to their visual
words (in seconds).
4. The time needed to compute the correspondences (in seconds).
5. The time needed by RANSAC to compute the pose (in seconds).
6. The total time needed for localization (in seconds), excluding loading of 
query features, adjusting the coordinate system of the 2D features such that
the origin coincides with the center of the image, and the computation of the
visual word assignments.

Details on the parameters of acg_localizer can be found by running it without
any parameters.
Using the method for kd-tree based localization (acg_localizer_knn) is very
similar, please type acg_localizer_knn for a description of its parameters.

Running the ECCV localization method is similar. Here is an example:

After computing the descriptor assignments using a vocabulary tree with the
command 
* compute_desc_assignments bundle.info 1 100000 clusters.txt bundle.desc_assignments.integer_mean.voctree.clusters.100k.bin 6 1

the acg_localizer_active_search can be used as follows
* acg_localizer_active_search list.query.keys.txt bundle.out 100000 clusters.txt bundle.desc_assignments.integer_mean.voctree.clusters.100k.bin 0 results.txt 200 1 1 1 10

The parameters are similar to acg_localizer:
  The first parameter is again the list of query images, similar to
acg_localizer.
  The second parameter is the filename of the Bundler reconstruction which
should be used as the 3D model for the localization method, i.e., the same
Bundler file that you used to compute the .info file with Bundle2Info.
  The third parameter is the number of visual words to use.
  The fourth parameter is the filename of a text file containing the cluster
centers.
  The fifth parameter is the assignment file previously computed using
compute_desc_assignments.
  The sixth parameter specifies the prioritization scheme to use to combine
2D-to-3D and 3D-to-2D 
matching. Setting the parameter to 0 chooses the "Combined" strategy
(recommended), setting it to 1 chooses the "Direct" strategy, and setting it to
2 chooses the "Afterwards" strategy. Please see the paper for an explanation of
the different strategies.
  The seventh parameter is the filename of a text file to which
acg_localizer_active_search will report timings and statistics. It has the same
format as the output file of acg_localizer.
  Parameter eight defines the number N_3D of nearest neighbors in 3D that should
be found during active search. 200 is the recommended value for that parameter.
  The last parameters specify various filtering methods:
  Parameter nine defines whether the RANSAC Pre-Filter described in the paper
should be used (1) or not (0).
We recommed using the filer, i.e., setting the parameter to 1.
  Parameter ten defines whether to use the point filter proposed in the paper.
Set to 1 (recommended) to enable the filter and to 0 otherwise.
  If parameter eleven is set to 0, the original images included in the
reconstruction are used for visibility filtering. If the parameter is set to 1,
then k cameras are clustered together, where k is defined by the last 
parameter. Again, more details can be found in the paper.


------------
Change Log
------------
version 1.0   - initial release
version 1.1   - updated this readme file since newer versions of FLANN break 
               the old API.
version 1.2   - added acg_localizer_active_search, which implements the 
                localization framework proposed in the ECCV 2012 paper.
                Fixed error that would compute the wrong descriptor indices 
                for datasets with more 33 million descriptors.
                Some performance optimization, resulting in faster localization
                times.  
version 1.2.1 - added support for the file format of the aachen.info file from 
                the Aachen dataset, which is slightly different to the file 
                format generated by Bundle2Info.
version 1.2.2 - Bugfix: ANN will throw an error if the number of 3D points in one 
                connected component is smaller than N_3D, i.e., the number of  
                3D points we are searching for. This bug should now be fixed.
                If you still experience it, please contact Torsten Sattler
                (tsattler@cs.rwth-aachen.de)

--------
Feedback
--------
We appreciate your feedback! If you have found bugs, have comments on the
software, questions regarding certain lines in the source code (or just any
basic question) please send them to tsattler@cs.rwth-aachen.de !


----------------
Acknowledgements
----------------
The authors want to thank Martin Habbecke for contributing his math and pose
estimation functions to this project. We also thank Jan Moebius for his help
with the CMake build system  and Darko Pavic for help with the timer class. 

 
