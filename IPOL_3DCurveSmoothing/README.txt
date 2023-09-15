================================================================================
############################## 3D CURVE SMOOTHING ##############################
================================================================================

# ABOUT
-------

* Author    : Luis Alvarez (lalvarez@ulpgc.es) ...
* Copyright : (C) 2009-2019 IPOL Image Processing On Line http://www.ipol.im/
* License   : CC Creative Commons "Attribution-NonCommercial-ShareAlike" 
              see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en

================================================================================

# CONTENTS
----------

  - Overview
  - Requirements
  - Compilation
  - Known Issues
  - Usage
  - Source Code Organization
  - Example
  - Thanks

================================================================================

# OVERVIEW
----------

This source code provides an implementation of a method to automatically smooth
3D curves.

 http://www.ipol.im/pub/art/2016/130/ <==== MODIFICAR ENLACE CUANDO TENGAMOS LA DEMO!!!!

The program produces 2 outputs:
 (1) A text file with the points of the smoothed curve.
 (2) An obj file with a comparison of the 3D curves (original and smoothed). If
     we provide a segmentation as input (optional), it will be also included in
     the obj file. To visualize it you can use any compatible viewer or the 
     online viewer available at http://3dviewer.net/

================================================================================

# REQUIREMENTS
--------------

The code is written in ANSI C++, and should compile on any system with
an ANSI C compiler.

================================================================================

# COMPILATION
-------------

We have checked our code on:
	- Windows 10 using  MinGW with GCC version 7.2.0 and Code::Blocks
	- Ubuntu 16.04.3 LTS with GCC version 5.4.0 and makefile
	- MacOSX 10.14.6 (Mojave) with Homebrew GCC version 9.2.0 and makefile (please
    see section Known Issues for more information)

To compile the code, simply use the provided makefile, with the command `make`. 
If you want to use the OpenMP library, please use `make OMP=1`. Optionally, you 
can obtain coverage information by including `COV=1` to the `make` command. Then 
you can generate an html document (after running the program) with the 
information by using the commands (or alternatively, use the script
generateCoverageInfo.sh included with the code):

	lcov --capture --directory . --output-file coverage.info
	lcov --remove coverage.info '/usr/include/*' -o coverage_filtered.info
	genhtml coverage_filtered.info --output-directory outcov

(NOTE: the example commands have been executed in the same directory in which 
the demo is placed. Moreover, the directory outcov must be created to locate the
output from the genhtml command. Bear in mind that, in order to cover most of
the code, you can perform several runs with the different options before
generating the coverage information)

Alternatively, you can manually compile by using the following command:

	g++ -Wall -Wextra -O3 -std=c++98 main/main3DCurveSmoothing.cpp 
                                   main/curve3D.cpp auxiliary/obj3D.cpp -lm -o 
                                   3DCurveSmoothing

	Or the following alternative compilation line for using the OpenMP library:

	g++ -Wall -Wextra -O3 -fopenmp -std=c++98 main/main3DCurveSmoothing.cpp 
                                            main/curve3D.cpp auxiliary/obj3D.cpp
                                            -lm -o 3DCurveSmoothing

(NOTE: to generate coverage information, remember to include the option
--coverage)

================================================================================

# KNOWN ISSUES
--------------

By default, macOS Mojave provides the frontend `clang` to compile `C`, `C++`,
`Objective-C`, ... code. However, clang is not able to deal with the use of the
OpenMP library (if you try to compile you receive a message of `clang: error:
unsupported option -fopenmp`). For this reason, we need to install `GCC` by
means of Homebrew (https://brew.sh/) following these steps:

  1. If you do not have already installed XCode, install it from App Store.
  2. Install the command line tools: open a terminal and write `xcode-select
     --install`, press enter, and follow the installation steps.
  3. Install Homebrew by writing in a terminal `/usr/bin/ruby -e "$(curl -fsSL
     https://raw.githubusercontent.com/Homebrew/install/master/install)"`.
  4. If you have installed XCode for the first time and never opened it, you
     will get a message in the terminal when you try to run Homebrew (for
     instance when running `brew install ...`). In this message, you will be
     asked to agree on the XCode terms or write a command to avoid it. Simply
     open XCode, accept the terms, and let it install some extra packages. Once
     the installation has finished, you can close XCode and pass to the next
     step.
  5. Once Homebrew is installed and the XCode terms have been accepted, you can
     install GCC by writing `brew install gcc` in your terminal.

The previous steps will install GCC in your computer, but instead of being in
the path `/usr/bin` it will be located at `/usr/local/bin`.

We have included a new option in the makefile called APPLE. When it is
activated, we use the GCC compiler located at `/usr/local/bin`. Please, change
such a path and `g++` version if it is needed. You can compile the code with
make as `make OMP=1 APPLE=1`. As an alternative, you can also create links to
such compiler in `/usr/bin`, however, bear in mind that in recent macOS versions
the writing on system paths is disabled (System Integrity Protection), which
requires to disable `csrutil` and reboot (under your responsibility).

As described in the previous section, you can alternatively compile by using the
following line:

  /usr/local/bin/g++-9 -Wall -Wextra -O3 -fopenmp -std=c++98 
              main/main3DCurveSmoothing.cpp main/curve3D.cpp auxiliary/obj3D.cpp
              -lm -o 3DCurveSmoothing

in which we indicate the path to the `gcc` compiler that we have installed
through Homebrew.

(NOTE: `clang` is able to compile without using OpenMP)

================================================================================

# USAGE
-------

This program takes 7 parameters (one of them is [optional]):

  exe_file input_curve.txt w tol max_iter output_smoothed_curve.txt 
                                             3Dcomparison.obj [segmentation.hdr]


  1. exe_file:                    name of the executable file (called 
                                  ./3DCurveSmoothing)
  2. input_curve.txt:             name of the ASCII file with the 3D curve point
                                  coordinates
  3. w:                           weight parameter to balance the energy
  4. tol:                         tolerance to stop iterations
  5. max_iter:                    maximum number of iterations
  6. output_smoothed_curve.txt:   filename for the output ASCII file with the
                                  smoothed 3D curve
  7. 3Dcomparison.obj:            output obj file to compare the original and
                                  smoothed curves
  8. [OPTIONAL] segmentation.hdr: 3D image (unsigned char in Analyze format)
                                  with a segmentation of a 3D object for which
                                  the given 3D curve is the centerline. If the
                                  image value is equal to zero the point is
                                  outside the 3D object.
                                
================================================================================

# SOURCE CODE ORGANIZATION
--------------------------

The source code is organized in the following folders:

  - auxiliary: inside this folder there is auxiliary code, namely:
    - color.h:            class to manage color values
    - dbh.h:              structs for the datatypes of Analyze 3D images (based
                          on Analyze format from Mayo Clinic (https://web.archive.org/web/20121116093304/http://wideman-one.com/gw/brain/analyze/formatdoc.htm))
    - image3D.h:          class to store 3D images
    - obj3D.h, obj3D.cpp: class to generate Wavefront .obj files
    - point3d.h:          class to store and manage 3D points

  - main: the main code which includes the algorithms described in the paper:
    - curve3D.h, curve3D.cpp:   methods to manage 3D curves
    - main3DCurveSmoothing.cpp: main program

  - example: example input file and results:
    - input_curve.txt: ASCII file containing the 3D points of the original curve
                       (contains a first line with the number of points and
                       afterward one point per line separated by spaces in float
                       type (x. y. z.)):
        494
        178.246400 194.540300 359.649600
        178.312631 194.715571 358.667310
        178.414299 194.916885 357.693074
        178.564228 195.167875 356.736766
        ... 
    - Segmentation.hdr and Segmentation.img: segmentation for which
      the curve is the centerline in Analyze format (unsigned char)
    - curve_smoothed.txt:     output ASCII file with the result of smoothing the
                              input curve without including the segmentation
    - curve_smoothed_seg.txt: output ASCII file with the result of smoothing the
                              input curve including the segmentation as input
    - 3DCurves_w1.obj:        .obj file with the comparison between input (red)
                              and output (blue) curves
    - 3DCurves_w1_seg.obj:    .obj file with the comparison between input (red)
                              and output (blue) curves including the
                              segmentation volume in gray
    - phantom_curve.txt: input ASCII file with the curve described by the
                         equation C(t) = ((6Pi-t)cos(t),(6Pi;-t)sin(t),3t)
                         for t in the interval [0,6Pi]
    - phantom_curve_with_noise.txt: the same curve as the previous file but with
                         noise added
    - curve_smoothed_phantom.txt: output ASCII file withe the result of
                         smoothing the phantom curve with noise
    - 3DPhantomCurves.obj: output file with the comparison between the phantom
                         curves (noisy vs smoothed).
  - documentation: Doxygen documentation in HTML format for further information

The files of c++ code inside the folder main will be subjected to the IPOL
peer-review process.

================================================================================

# EXAMPLE
---------

You can test the program with the provided test file (input_curve.txt) in the 
following way:

  ./3DCurveSmoothing example/input_curve.txt 1. 0.00001 1000 
    example/output_curve_smoothed.txt example/3DCurvesComparison_w1.obj

or alternatively, you can also indicate the segmentation (in Segmentation.hdr)
associated with the curve:

  ./3DCurveSmoothing example/input_curve.txt 1. 0.00001 1000
    example/output_curve_smoothed.txt example/3DCurvesComparison_w1.obj
    example/Segmentation.hdr

Furthermore, you can compare your results with the results present inside the
folder example.

================================================================================

# THANKS
--------

The authors would be grateful to receive any comment, especially about
portability issues, errors, bugs or strange results.
