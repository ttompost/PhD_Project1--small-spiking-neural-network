# (Ching et al., 2010) TC, RE, and synaptic DynaSim mechanism files.

 DynaSim-compatible mechanism files for simulation of the thalamus of (Ching et al., 2010).

Adding these mechanism files and associated functions into where you keep your mechanism files for [DynaSim](https://github.com/DynaSim/DynaSim), e.g. `/your/path/to/dynasim/models`, should enable you to simulate the computational thalamus from:

Ching, S., Cimenser, A., Purdon, P. L., Brown, E. N., & Kopell, N. J. (2010). Thalamocortical model for a propofol-induced alpha-rhythm associated with loss of consciousness. Proceedings of the National Academy of Sciences, 107(52), 22665–22670. http://doi.org/10.1073/pnas.1017069108

Note that this code diverges from the given equations due to typos and errors in the original equations. The code contained here has been ground-truthed and cross-checked across both the original code run for the paper and the source material that THAT code is based on, several times. In other words, the code contained here is **more** correct than the equations in the paper.

## Install

The easiest way to get started with this is

1. Install DynaSim
    1. `git clone` [DynaSim](https://github.com/DynaSim/DynaSim) to where you want to put the code
    2. Add `addpath(genpath('/your/path/to/dynasim'))` to your
         `~/Documents/MATLAB/startup.m` file. Create the necessary folders and files
         if need be.
3. `git clone` this repo into '/your/path/to/dynasim/models', i.e. the 'models'
     subdirectory of your copy of the DynaSim repo.
4. Set your own data directory inside the script, make your own changes, etc.
5. Believe it or not...that should be it! You should be able to start MATLAB in
     any directory and run this script successfully! Let me know if there are
     problems, at austin.soplata 'at symbol' gmail 'dot' com

Note that there are extra synaptic mechanism files unused in the `database` that do things like simulate a persistent AMPAergic spike train down to the cells, etc., in order to mimic the cortex's input.
