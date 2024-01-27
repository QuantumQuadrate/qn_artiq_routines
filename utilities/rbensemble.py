"""
RbEnsemble
Preston Huft, Summer 2019

A class for generating Rb atom distributions to be passed to various simulations.
"""

## libraries
# from matplotlib import pyplot as plt
import numpy as np
from numpy import linspace,sin,cos,log,exp
import math as m
from math import sqrt,pi,e,exp
from random import random as rand

## local files
from utilities.physconsts import  *
from utilities.rbconsts import *

## classes 

class RbEnsemble:
	""" Ensemble of Rb atoms with a momentum distribution specified
		by temperature. Optional spatial distribution.
		
		'T'= temperature
		'size'= number of atoms
		'xdist': optional parameter specifying the initial
			position distribution
	"""	
	
	def __init__(self,T,size=None,xdist=None,statedist=None,m=mRb,units=1):
		"""
		T
		size
		xdist
		statedist
		m=mRb is the particle mass in kg. 87Rb mass by default.
		units=1 length units. 1 for meters, 1e-3 is mm, 1e-6 for um, and so on
		"""
		
		self.m = m
		self.temp = T
		self.units = units

		# For efficiency, pre-generate a specified number of atoms
		if size is not None:
			self.size = size
			self.v = self.sampling_maxboltzv(self.size,[0,1],self.temp) # rms
			self.p = self.m*self.v # rms
			self.x = np.empty(self.size)
			if xdist is None:
				self.x = np.zeros(self.size)
			elif xdist == 'normal':
				self.x = np.random.normal(0,size=self.size)
			if statedist is not None:
				self.amplitudes = self.psi_coeffs(self.size,statedist)
		else:
			self.size = 0
			self.v = np.array([]) # rms
			self.p = np.array([]) # rms
			self.x = np.array([])
			self.amplitudes = np.array([])
			
	def phasespace(self):
		""" Plots the ensemble in phase space. 1D x and p only for 
			now.
		"""
		
		m = self.m
		
		xmax = max(self.x) # like xmas but better
		xmin = min(self.x) # because i said so
		dx = xmax-xmin
		
		pmax = max(self.p)/m
		pmin = min(self.p)/m
		dp = pmax-pmin
		
		#fig, ax = plt.subplots()
		#ax.scatter(self.p/m,self.x)#,linestyle=None)
		#ax.set(xlabel='p [m/(s m)]', ylabel='r [arb]',
		#	   xlim=(pmin-.1*dp,pmax+.1*dp),
		#	   ylim=(xmin-.1*dx,xmax*+.1*dx))
		#plt.show()
		
	def vpt(self):
		""" Return a speed from Maxwell-Boltzmann dist. """
		return self.sampling_maxboltzv(1,[0,1],self.temp)

	def xpt(self,domain):
		""" Return a position from a flat dist by default. """
		
		x1,x2 = domain
		x = rand()*(x2-x1) # only works for x1,x2 > 0
		return x

	def maxboltzv(self,T,v,normalization=False):
		""" Maxwell-Boltzmann distribution of speeds for 3-dimensional
			gas. Returns f(v) for T. """
		m = self.m

		A = 4*pi*(m/(2*pi*kB*T))**(3/2) # normalization consts
		meanv = sqrt(2*kB*T/m)/self.units


		if normalization is True:
			return A
		else:
			return A*v**2*exp(-(v/meanv)**2)

	def sampling_maxboltzv(self,size=None,domain=[0,1],T=None,vectorial=False,showplot=False):
		""" Sample random speeds with a Maxwell-Boltzmann dist. 
			'size': sample size
			'domain': [v1,v2] the restricted domain of the pdf; e.g.
				a Guassian goes to zero well before inf so we could
				let the domain be a finite region
			'T': temperature
			'vectorial': 
				If False, only return a scalar. 
				Set to True to return velocity vectors with a 
				direction from a flat distribution. 
		"""
		m = self.m
		print(self.temp)
		if size is None:
			size = self.size
		if T is None:
			T = self.temp
		
		n = size 
		v1,v2 = domain # capped at 1 [m/s] by default
		print(kB)
		print(T)
		print(m)
		print(self.units)
		mean = sqrt(2 * kB * T / m) / self.units
		mean = sqrt(2*kB*T/m)/self.units
		fmax = self.maxboltzv(T,mean) # the maximum
		y_dist = np.empty(n) 
		f_dist = np.empty(n)
		v_dist = np.empty(n) # this is the distribution we want
		j = 0 # dist index
		while j < n:
			v = (v2-v1)*rand() # rand val on domain of f(x)
			f = self.maxboltzv(T,v)
			y = rand()*fmax # rand val on range of f(x)
			if y <= f:
				y_dist[j]=y
				f_dist[j]=f
				v_dist[j]=v # x vals with approximate gaussian pdf
				j+=1

		# plot distribution as a check:
		#if showplot is not False:
		#	plt.scatter(v_dist,y_dist,c='red',s=10)
		#	plt.scatter(v_dist,f_dist,c='blue',s=10)
		#	plt.show()

		return v_dist
		