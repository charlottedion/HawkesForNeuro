#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 11:13:52 2020
@author: charlotte
"""


from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from random import *
from numpy import abs, sin, cos, arange, pi, exp, linspace, random, sqrt, log, loadtxt

########################################
# Utilisation de tick, package martin Bompaire et Bacry   https://x-datainitiative.github.io/tick/modules/hawkes.html
# pour installation copier: pip install tick

from tick.plot import plot_hawkes_kernels,  plot_hawkes_kernel_norms
from tick.hawkes import SimuHawkesSumExpKernels,SimuHawkesExpKernels, SimuHawkesMulti, HawkesSumExpKern, HawkesADM4, SimuHawkes, HawkesKernelTimeFunc, HawkesExpKern, HawkesEM
from tick.base import TimeFunction
from tick.plot import plot_hawkes_baseline_and_kernels
import tkinter
from tempfile import TemporaryFile

import math as m
from math import exp, log
import scipy
from scipy import stats
from scipy import optimize
from scipy.optimize import fmin
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.stats import randint
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot_2samples

import random

from pandas import * 
import os 
import glob

import seaborn as sns
a = sns.color_palette("colorblind")
sns.palplot(sns.color_palette("colorblind"))

from scipy.io import loadmat

import random 
########################################################################################################################################################
# I Simulation of the data
########################################################################################################################################################



M = 8#20# number of neurons
n = 50# number of trials
Tmax = 100# horizon time observation

random.seed(12)

adj = np.diag(np.ones(M)*0.7)
adj[random.sample(range(M), k=int(M/2)), random.sample(range(M), k=int(M/2))] = 0.3
print(adj)

#np.savetxt('results/adj_true.txt', adj)


w, v = np.linalg.eig(adj)
print('Spectral radius true matrix', max(abs(w)))


decay = np.ones(shape=(M,M))*5
baseline = np.ones(shape=M)*0.5
#np.savetxt('results/baseline.txt', baseline) 

hawkes = SimuHawkesExpKernels(adjacency=adj, decays=decay, baseline=baseline, end_time=Tmax, verbose=False)

multi = SimuHawkesMulti(hawkes, n_simulations=n)

multi.end_time = [Tmax]*n
multi.simulate()
events = multi.timestamps

#### Description 

plt.close('all')


leng = []
for k in range(n):
    for i in range(M):
        leng = np.append(leng,len(events[k][i]))
        

hist_alltrials, ax = plt.subplots(1,1, figsize=(10,5)) 
ax.hist(leng, bins='auto', color=a[0], alpha=0.7, rwidth= 0.85, edgecolor='black')
ax.set_xlabel('Number of spikes')
ax.set_ylabel('Probability')
plt.savefig('images/histo_alltrials.pdf')


plt.close('all')

num = 1
plotspike, ax= plt.subplots(1, 1, figsize= (10,8))
for i in range(M):
    ax.plot(events[num][i], np.zeros(len(events[num][i]))+i, linestyle='', marker='+', color=a[0])
ax.set_xlabel('Spike trains')
ax.set_ylabel('Neuron')
plt.savefig('images/spikesSIMU_trial1.pdf')


#### Saving data

for k in range(n):
 events_csv = [None]*M
 for i in range(M):
    events_csv[i] = np.zeros(int(max(leng)))
    events_csv[i][range(len(events[k][i]))] = events[k][i]
    
 np.savetxt('data/events'+str(k)+'.csv', events_csv, delimiter=',', fmt='%s')

# one file .csv contains the spike trains of one trial


########################################################################################################################################################
# II Estimation with ADM4
########################################################################################################################################################

grid_decay = [0.5, 1., 10., 20., 50.] # grid in which we search the decay

score = [0]*len(grid_decay)

for i in range(len(grid_decay)):
    learnerADM4 = HawkesADM4(decay= grid_decay[i], lasso_nuclear_ratio = 1)
    learnerADM4.fit(events)
    score[i] = HawkesADM4.score(learnerADM4)

decayADM4 =  grid_decay[np.argmax(score)]
learnerADM4 = HawkesADM4(decay = decayADM4, lasso_nuclear_ratio = 1)
learnerADM4.fit(events)

adjADM4 = learnerADM4.adjacency
baselineADM4 = learnerADM4.baseline

#np.savetxt('results/adjADM4.txt', adjADM4) 

####### Graphe
learnerADM4 = plot_hawkes_kernel_norms(learnerADM4)
plt.savefig("images/learnerADM4.pdf")
#######

#print('% of coeff non-zero true matrix', np.sum(adj>0)*100/(M*M))
#print('% of coeff non-zero estimated matrix', np.sum(adjADM4>0.001)*100/(M*M))

print('Frobenius norm true matrix', np.linalg.norm(adj, 'fro'))
print('Frobenius norm estimated matrix', np.linalg.norm(adjADM4, 'fro'))

w, v = np.linalg.eig(adj)
print('Spectral radius true matrix', max(abs(w)))

w, v = np.linalg.eig(adjADM4)
print('Spectral radius estimated matrix', max(abs(w)))

########################################################################################################################################################
# III Goodness-of-fit test
########################################################################################################################################################

#------ one estimation on each trial

grid_decay = [0.1, 1., 3., 10., 20. ,50.]
decayADM4onetrial = [None]*n
adjADM4onetrial = [None]*n
baselineADM4onetrial = [None]*n

for k in range(n):
 score = [0]*len(grid_decay)
 for i in range(len(grid_decay)):
    learnerADM4 = HawkesADM4(decay= grid_decay[i], lasso_nuclear_ratio = 1)
    learnerADM4.fit(events)
    score[i] = HawkesADM4.score(learnerADM4);
    
 decayADM4onetrial[k] =  grid_decay[np.argmax(score)]
 
 learnerADM4k = HawkesADM4(decay = decayADM4onetrial[k], lasso_nuclear_ratio = 1)
 learnerADM4k.fit(events[k])
 adjADM4onetrial[k] = learnerADM4k.adjacency
 baselineADM4onetrial[k] = learnerADM4k.baseline;
 


#-------

def expcdf(t, beta=1.0):
    return 1 - np.exp(-t/beta);

def validation(neuron, sub):
    # neuron: numero du neuron que l'on regarde
    # sub: indice des sous echantillon de trials 
   residualproc_subech= [None]*len(sub)
   
   for indsub in range(Nsub):
        subjumpvec = np.r_[events[sub[indsub]][neuron], Tmax] 
        Njump = len(subjumpvec) 
        residualproc = np.zeros(Njump)
    
        for i in range(Njump): 
            residualproc[i] = baselineADM4onetrial[indsub][neuron] * subjumpvec[i]
            vecl = subjumpvec[np.where(subjumpvec < subjumpvec[i])]
            residualproc[i] +=  adjADM4onetrial[indsub][neuron, neuron] * (len(vecl) - sum(np.exp(- decayADM4onetrial[indsub] * (subjumpvec[i] - vecl))));
            # intercation
            autreneurons = list(range(M))
            del autreneurons[neuron]
            
            for j in range(len(autreneurons)):
            
                subjumpvecj = events[sub[indsub]][autreneurons[j]]
                veclj = subjumpvecj[np.where(subjumpvecj < subjumpvec[i])]
                residualproc[i] += adjADM4onetrial[indsub][neuron, autreneurons[j]] * (len(veclj) - sum(np.exp(- decayADM4onetrial[indsub] * (subjumpvec[i] - veclj))));
              
        residualproc_subech[indsub] = residualproc
   return residualproc_subech
    

#############################################################################################################################################

Nsub = int(np.floor(n**(2/3))) # size of the sub-sample used for the test

npas= 100
u = linspace(0,1, num = npas) #discretization of the sup
Nrep = 100 # number of repetition for the acceptation rate
coeff = 0.5
#------
Ztest = np.zeros((Nrep, M))
for rep in range(Nrep):
    sub = np.sort(random.sample(range(n), k = Nsub))# indices sélectionnés    
    
    test = [None] * M
    for neu in range(M):
        test[neu] = validation(neuron = neu, sub = sub); 
 
    for neu in range(M):
        residual = test[neu][0][range(len(test[neu][0])-1)] 
        if len(test[neu][0])==0:
            theta_max = 0;
        else: 
            theta_max =  max(test[neu][0]);
        
        for i in range(1,len(sub)):
                if len(residual)==0:
                  residual = np.r_[residual, test[neu][i][range(len(test[neu][i])-1)]+ theta_max];
                else:
                  residual = np.r_[residual, test[neu][i][range(len(test[neu][i])-1)]+ theta_max];
                if len(test[neu][i])==0:
                 theta_max += 0;
                else: 
                 theta_max +=  max(test[neu][i]);
      
        theta_max = (1/Nsub) * theta_max
        theta = theta_max * coeff
        seuil= Nsub * theta  
                
        resint = residual[residual<=seuil]
        
        num = len(resint)
        if num == 0:
            num = float(1)
                
        vectest = np.zeros(npas)
        for k in range(npas):
            vectest[k] = abs( (1/num) * sum(resint/seuil <= u[k]) - u[k]);
        
        Ztest[rep, neu] = sqrt(num) * max(vectest);

#------
alpha = 0.05  #risk level
KSquant = 0.325 # 0.325 for n= 50, Nsub = 13  #see tables,1.92/sqrt(Nsub) # when Nsub>40
                # 0.259 for n= 100, Nsub = 21 
      
quant = KSquant * (sqrt(Nsub) + 0.12 + 0.11/sqrt(Nsub)) #empirical quantile, #KSquant*(sqrt(Nsub)+0.12+0.11/sqrt(Nsub)) if $Nsub<45

res = np.mean(Ztest<quant, axis=0); # ACCEPTATION RATE 
print(res)    
#------

# coeff = 0.5
# [0.75 0.97 0.48 1.   0.86 0.68 0.7  0.42]

#----- Table
fig, ax = plt.subplots(figsize = (5,5))
datatable = np.zeros(shape=(M, 2))
datatable[:,0] = range(1,M+1)
datatable[:,1] = res
collabel=("Neuron", "Acceptation rate")
ax.axis('tight')
ax.axis('off')
the_table = ax.table(cellText = datatable, colLabels=collabel,loc='center')
plt.savefig('results/table_acceptationrate.pdf')



