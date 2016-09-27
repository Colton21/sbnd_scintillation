# sbnd_scintillation
Low energy light simulations

To use this project you will likely need access to C++ 11 features, such as the custom sorting function, where a pair
is sorted based on the first element, rather than a simple array of doubles. It also makes use of features inherent to 
CERN ROOT, but this is out of familiarity and that the library files (not included in the project) are in the .root format.

As the code stands, it can be compiled simply by typing "make -B", to recompile everything. This only takes a few seconds.
As many of the configurable parameters are in 'libraryanalyze_light_histo.h', when you change a parameter in the .h file,
I prefer to recompile everything in the project.

The Makefile genrates and executable that can be run with "./libraryanalyze_light_histo", or whatever you change the name to.
If you happen to be missing the data file, a segmentation violation will occur. Before the crash readout, you will find
that the requested file could not be found. Change your path, and it should then run fine.

At the moment the code creates two root files - where the event_file.root should contain the information needed to perform any
analysis. The event_tree has data on an event-by-event basis, and data_tree has the information based on detected photons.

Notes on number of events and memory consumption:

At the time of writing, the memory consumption of the code is rather significant. A suspicion has to do with how the timing 
parameterisation distributions are rather finely sampled, but this is untested. Efforts have been made already to attempt
to reduce overall memory consumption, but it is possible that a more efficient method could be developed.

A current benchmark is that for the full foils configuration (most light), running for 10 readout windows, 12 ms, gives
around 670 events (decays), and at the end of the file generation requires less than 4GB of memory, despite the output files
themselves being rather small.

To do an isolated study, this may be sufficient, however statistically small samples run the risk of undersampling our
simulation space, such that the fluctuation of the events is significant. For this reason, I have been generating between
4500-5000 data files (12ms) for a timing between 54-60s, giving over 3 million events. This means that each voxel in the
SBND detector (5cm^3 voxelisation) has around 10 decays each.

The variation in detection efficiency based on position is important, making higher statistics samples very useful.

In the case of supernova events, as each event produces significantly more light, I can typically generate around 150 events
per file, with memory consumption under 4GB.
