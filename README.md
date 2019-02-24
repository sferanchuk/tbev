# tbev_compare

Algorithms and utilities to support a comparative study of proteins for substrains of tick-borne encephalitis virus.

## List of folders

md - scripts for molecular dynamics simulations

comppdb - software for a comparison of MD trajectories

pkfs - pieces of C++ code for protein modeling

drawalign - several utilities to represent alignments/sequences and per-position distributions generated by a comparison of MD trajectories

alignmut - utilities to represent alignments and statistics of per-position substitutions

vtkview - ulitilty to draw protein structure together with a distribution of NMA vector, implemented within VTK 3D graphics environment 

bubbleview - utility to represent a comparative distribution of NMA vectors, implemented with the use of d3.js functionality

## Dependencies and comments

Several of provided pieces of software have unresolved dependencies, due to ethical and/or legal restrictions which prevent to open the missing parts of code.

The parties involved into private relations with the author are not listed here. But, in addition, following dependencies are required to compile and use the utilities listed above:

md: amber11 package
comppdb: u3b.f - optimal superposition of 3D objects
pkfc: u3b.f, libpdb++, libSAS
vtkview: VTK library, Bio3D
bubbleview: d3b.js, Bio3D



