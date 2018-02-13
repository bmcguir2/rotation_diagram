#####################################################
#                                                   #
#               rotation_diagram.py README          #
#                                                   #
#####################################################

Last Updated: Sept 29, 2017

This program takes an input file containing data on a molecule and observed transitions and does a rotation diagram fit to the data.

#####################################################
#                                                   #
#                   Getting Started                 #
#                                                   #
#####################################################

Actually running the program is trivial once the input file is set up:

> make_rot_diagram(input_file)

An example input file is provided.  It follows the following format:

First Line: Molecule Name (LaTeX formatting accepted)
Second Line: Rotational Constants (space or tab delimited).  If linear: give one. If symmetric (B=C), give 2 (A B).  If asymmetric, give three (A B C)
Third Line: Optional column labels to keep things straight: #Freq		#dT		#t_err	#Eup	#g	#Aij	#dV	#V_err	#tbg
Fourth+ lines: The data as outlined above, space or tab delimited.

The data are:

Frequency (MHz)
Intensity (K)
Intensity Uncertainty (K)
Upper State Energy (K)
Degeneracy
Aij (s-1)
Linewidth (km/s)
Linewidth Uncertainty (km/s)
Background Continuum Temperature (K)

The program will look at all lines that don't start with '#', so you can comment out lines as needed.  They do not need to be in frequency order.

That's it!  If you find any issues, please let me know.
