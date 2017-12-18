#!/usr/bin/env python

#############################################################
#						Revision History					#
#############################################################

# reads in line parameters and creates a rotation diagram

# 0.0 - Project start
# 0.5 - switching to weighted error analysis
# 1.0 - adds ability to scale partition function 'manually'
# 1.1 - adds ability to check calculated partition function

#############################################################
#							Preamble						#
#############################################################

import sys

#Python version check

if sys.version_info.major != 3:

	print("This code is written in Python 3.  It will not execute in Python 2.7.  Exiting (sorry).")
	
	quit()

import numpy as np
from numpy import exp as exp
import time as tm
import warnings
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import itertools
import pylab as pl
from datetime import datetime
from scipy.optimize import curve_fit
#warnings.filterwarnings('error')

version = 1.1

h = 6.626*10**(-34) #Planck's constant in J*s
k = 1.381*10**(-23) #Boltzmann's constant in J/K
kcm = 0.69503476 #Boltzmann's constant in cm-1/K
ckm = 2.998*10**5 #speed of light in km/s
ccm = 2.998*10**10 #speed of light in cm/s

molecule = ''

#############################################################
#							Warning 						#
#############################################################

print('\nWarning! This code is in beta. I am fairly sure everything works, but, well, check your results, please.') 

#############################################################
#							Calls	 						#
#############################################################

pl.rc('text',usetex=True)

#############################################################
#							Functions						#
#############################################################
	
#read_cat reads the catalog file in

def read_molecule(molecule_file):

	'''
	Reads in a molecule file line by line
	'''

	my_array = []

	with open(molecule_file) as input:
	
		for line in input:
		
			my_array.append(line)		
			
	return my_array	
	
#parse_molfile takes the read in raw array from the molecule_file and splits it out into the appropriate variables

def parse_molfile(raw_array):

	#get the molecule name

	name = raw_array[0].strip('\n')
	
	#load in the rotation constants into a string to then split out and turn into floats
	
	temp = raw_array[1].strip('\n').split()
	
	rot_cons = []
	
	for i in range(len(temp)):
	
		rot_cons.append(float(temp[i]))
		
	frequency = []
	dT = []
	dT_err = []
	Eup = []
	g = []
	Aij = []
	dV = []
	dV_err = []
	tbg = []
	
	for i in range(2,len(raw_array)):

		if raw_array[i] == '' or raw_array[i] == '\n' or raw_array[i][0] == '#':
		
			continue
			
		else:
	
			temp = raw_array[i].strip('\n').split()
	
			frequency.append(float(temp[0]))
			dT.append(float(temp[1]))
			dT_err.append(float(temp[2]))
			Eup.append(float(temp[3]))
			g.append(float(temp[4]))
			Aij.append(float(temp[5]))
			dV.append(float(temp[6]))
			dV_err.append(float(temp[7]))
			tbg.append(float(temp[8]))
		
	frequency = np.asarray(frequency)
	dT = np.asarray(dT)
	dT_err = np.asarray(dT_err)
	Eup = np.asarray(Eup)
	g = np.asarray(g)
	Aij = np.asarray(Aij)
	dV = np.asarray(dV)
	dV_err = np.asarray(dV_err)
	tbg = np.asarray(tbg)	
	
	return name,rot_cons,frequency,dT,dT_err,Eup,g,Aij,dV,dV_err,tbg
	
#calc_Q uses the Gordy & Cook approximations for linear, symmetric, or asymmetric top

def calc_Q(T,rot_cons):

	Q = 0.0
	
	#linear molecule (1 rotational constant)

	if len(rot_cons) == 1:
	
		B = rot_cons[0]*10**6
	
		Q = (k * T)/(h * B)
		
	#symmetric top (2 rotational constants: A, B=C)
	
	if len(rot_cons) == 2:
	
		A = rot_cons[0]
		B = rot_cons[1]
	
		Q = ((5.34/3) * 10**6) * (T**3/(B**2 * A))**(0.5)
		
	#asymmetric top (3 rotations constants: A, B, C)
	
	if len(rot_cons) == 3:
	
		A = rot_cons[0]
		B = rot_cons[1]
		C = rot_cons[2]
	
		Q = ((5.34) * 10**6) * (T**3/(A * B * C))**(0.5)	
		
	return Q	

#fit_func is just the linear function we're going to fit the line to.  Required to do a full curve_fit optimization with errors.

def fit_func(x, m, b):

	y = m*x + b
	
	return y
	
#do_fit does the least squares fit to the data	

def do_fit(rot_cons,frequency,dT,dT_err,Eup,g,Aij,dV,dV_err,tbg,Qs):

	#do the y-axis calculation; k in J/K, dV in km/s, frequency in MHz, Aij in s^-1; conversions done in the equation
	
	A = (8 * np.pi * k * (frequency*10**6)**2 * 0.5 * (np.pi/np.log(2))**(0.5) * dT * dV*10**5)/(g* h * ccm**3 * 10**(Aij))
	
	y = np.log(A)
	
	#calculate the error inside the logarithm
		
	A_err = A*np.sqrt((dT_err/dT)**2 + (dV_err/dV)**2)
	
	#propagate the error through the logarithm
	
	y_err = A_err/A
	
	x = np.copy(Eup)
	
	Nt_guess = np.log(1E13/g[0])
	
	T_guess = 150
	
	Q = calc_Q(T_guess,rot_cons)*Qs
	
	m_guess = -1/T_guess
	
	b_guess = np.log(Nt_guess) + np.log(Q)
	
	coeff, var_matrix = curve_fit(fit_func, x, y, [m_guess,b_guess], sigma=y_err)
	
	m = coeff[0]
	
	b = coeff[1]
	
	err_matrix = np.sqrt(np.diag(var_matrix))
	
	m_err = err_matrix[0]
	
	b_err = err_matrix[1]
	
	T_rot = -1/m
	
	T_rot_max = -1/(m - m_err)

	T_rot_err = abs(T_rot_max - T_rot)
	
	Q = calc_Q(T_rot,rot_cons)*Qs
	
	Nt = np.exp(b + np.log(Q))
	
	Nt_max = np.exp(b + b_err + np.log(Q))
	
	Nt_err = Nt_max - Nt
	
	return Nt, T_rot, y, x, m, b, y_err, T_rot_err, Nt_err
	
#make_plot plots the results

def make_plot(y, Eup, m, b, Nt, T_rot, name, y_err, T_rot_err, Nt_err):

	plt.ion()
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	plt.title(name,size='large')
	
	plt.xticks(size='large')
	plt.yticks(size='large')
	
	minorLocator = AutoMinorLocator(5)
	plt.xlabel('E$_u$ (K)',size='large')
	plt.ylabel('ln(N$_u$/g$_u$)',size='large')

	plt.locator_params(nbins=4) #Use only 4 actual numbers on the x-axis
	ax.xaxis.set_minor_locator(minorLocator) #Let the program calculate some minor ticks from that

	ax.get_xaxis().get_major_formatter().set_scientific(False) #Don't let the x-axis go into scientific notation
	ax.get_xaxis().get_major_formatter().set_useOffset(False)		
	
	plt.errorbar(Eup,y,yerr=y_err,fmt='o',color='black',capthick=1,capsize=3)
	
	fit_y = np.copy(y)
	
	fit_y = m*Eup + b
	
	ax.plot(Eup,fit_y,color='red')
	
	N_mag = int(np.floor(np.log10(Nt)))
	
	N_c = Nt * 10**(-N_mag)
	
	if round(N_c,2) >= 10.00:
	
		N_mag += 1
		N_c /= 10
	
	Nt_err_c = Nt_err / 10**N_mag
	
	N_mag_str = '{}{}{}' .format('{',N_mag,'}')
	
	plt.figtext(.15,.20,'$N_T$ = {:.1f}({:.1f}) x 10$^{}$ cm$^{}$' .format(N_c,Nt_err_c,N_mag_str,'{-2}'),size='large')
	plt.figtext(.15,.15,'$T_{}$ = {}({}) K' .format('{ex}',int(T_rot),int(T_rot_err)),size='large')
	
	plt.show()
	
#make_rot_diagram is a wrapper that takes an input file and makes a diagram

def make_rot_diagram(molecule_file,Qs=1.0):
	
	global molecule
	
	molecule = molecule_file
	
	raw_array = read_molecule(molecule_file)
	
	name,rot_cons,frequency,dT,dT_err,Eup,g,Aij,dV,dV_err,tbg = parse_molfile(raw_array)
	
	Nt, T_rot, y, x, m, b, y_err, T_rot_err, Nt_err = do_fit(rot_cons,frequency,dT,dT_err,Eup,g,Aij,dV,dV_err,tbg,Qs)
	
	make_plot(y, Eup, m, b, Nt, T_rot, name, y_err, T_rot_err, Nt_err)
	
	print('NT = {:.2e} +/- {:.2e} cm-2, T_rot = {:.2f} +/- {:.2f} K' .format(Nt,Nt_err,T_rot,T_rot_err))
	
#check_Q provides a look at the calculated Q value, to guide the user to the appropriate scaling factor

def check_Q(T):

	raw_array = read_molecule(molecule)
		
	name,rot_cons,frequency,dT,dT_err,Eup,g,Aij,dV,dV_err,tbg = parse_molfile(raw_array)
	
	return calc_Q(T,rot_cons)
	
#############################################################
#							Run Program						#
#############################################################


