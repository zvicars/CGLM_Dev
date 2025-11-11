#Project Description
This project contains the software necessary to carry out and analyze lattice-based
Monte-Carlo simulations of water near surfaces. This includes reading input files, evolving trajectories,
outputting files that contain trajectory information (binary) or information regarding one of several 
"observation volumes" for calculating free energies. 

#Installation
    Compiling and installing this package follows the standard cmake approach.
    Via terminal in the top directory of the package (same folder as readme):
        "mkdir build"
        "cd build"
        "cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH"
        "make -j 1"
        "make test"
        "make install"
    The binaries will be present in $INSTALL_PATH/bin.

#Binaries and their functionality
    CGLM: Primary driver for carrying out simulations.
        Inputs: <string>, filepath for input JSON file
        Outputs: various, depends on choice of parameters in input JSON file. typically a binary trajectory and 
        additional outputs depending on applied biasing potentials
   
    analysis: Primary driver for post-processing simulation data, contains a bunch of algorithms for
    analyzing lattice-based simulation results. 
        Inputs: <string>, filepath for the input JSON file
        Outputs: Various, depends on the type and number of post-procecssing calculations performed
    
    b2xyz: A utility for converting binary trajectory files produced via CGLM to and xyz
    molecule file format (https://en.wikipedia.org/wiki/XYZ_file_format) that can be visualized 
    easily. Preferred visualization software is Ovito for this data (3.7.12 or older preferred 
    due to paywalled features in later versions)
        Inputs: Hyphenated argument list
            -b <int>, what frame in the binary trajectory to start on
            -e <int>, what frame in the binary trajectory to end on
            -s <int>, frequency, write every s'th frame
            -f <string>, input file name (binary trajectory file from CGLM)
            -o <string>, output file name (xyz file)
        Outputs: Plaintext xyz file at the path specified using -o.

    genphi: A utility for constructing a "phi" file, which is used as an an input for the CGLM program.
    A phi file contains information about solid surfaces, represented as a per-voxel external potential, 
    phi(x,y,z), which penalizes or promotes the occupancy of a single voxel.
        Inputs: A hyphenated argument list
            -box <int> <int> <int>, number of voxels in x, y, z. Total number of voxels is x*y*z
            -gs <float>, side-length for each voxel, actual box size is (gs*x, gs*y, gs*z)
            -i <string>, input file containing contributions to phi field (rectangles, spheres, cylinders with
            associated potential energy functions)
            -o <string>, output file (binary file containing external potential information)
            -int <int>, optional input to average potential functions at int*int*int different points per voxel
            to obtain smoother potential fields
        Output: Binary file containing external potential information, serves as an input to CGLM.
    
    phi2xyz: A utility for converting binary phi files to xyz files that can be readily visualized (ovito preferred).
        Inputs: Filepaths pointing to the input location and desired output location, -f <input file> -o <output file>  
        Output: An xyz formatted plaintext file containing x, y, z, phi, where phi is a float value that specifies the 
        value of the field at coordinates x, y, z.

    gro2phifile: A utility for converting gromacs GRO files to genphi input files.
        Inputs: A single argument containing the path to a JSON-formatted file containing additional input arguments (samples provided)
        Outputs: A text file that can then serve as an input to "genphi"

#additional notes and example files
A folder containing example calculations will be present in the "examples" directory in the folder specified using "DCMAKE_INSTALL_PREFIX". 
This allows the accompanying scripts to function without the user having to specify the location of the binaries. These examples will use 
python for data analysis, so an install of python3 will be necesary (with some other packages like numpy). Each example will contain a short
readme specifying what it does and how to visualize the data produced.