# Hydrogen in Blender
Generates blender scenes of various hydrogen orbitals in the DeBroglie-Bohm picture.

## Installation
- Make sure the python blender uses has the packages ``numpy``, ``scipy`` and ``sympy`` installed.
- For Linux:
  - Go to ``path/to/your/blender/installation/python/bin``
  - Open a terminal and type: ``./python3.7m -m ensurepip``
  - Install the packages with: ``./python3.7m -m pip install numpy scipy sympy``
- Start blender and make a new script file
- Copy the contents of ``script.py`` to your script
- Change the parameters to your liking and run the script

## Parameters
-``n``,``l`` and ``m`` are the quantum numbers of the orbital you want to visualize
- ``num_points`` is the number of electron trajectories that will be calculated (also the number of blobs you can later see)
- ``num_frames`` is the number of frames that will be animated
- ``end_time`` is the total physical time that passes in ``num_frames`` frames
- ``framestep`` is the intervall at which keyframes will be taken
- ``r``, ``theta``, ``phi`` is the discretization of space. Here you can mainly change number of points along each axis, but consider that increasing the number of points increases the memory requirement of the program drastically
