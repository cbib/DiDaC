# coding=utf-8
from operator import itemgetter

__author__ = 'hayssam'
import scipy as sp
import scipy.stats as st
import numpy as np


# calcul de la p-valeur one sided, on s'attend à ce que la valeur observé soit plus grande que la valeur random

observed = 10
pop =np.append(np.random.normal(10,1,1000),observed)
observed_z=st.mstats.zscore(pop)[-1]
observed_z_pval = 1-st.norm.cdf(observed_z)
# empirical p-value, like
larger_than_observed = np.sum(pop>=observed)
emp_p_val = larger_than_observed*1.0/len(pop)


# on simule 1000 fois, en comparant a chaque fois la p-valeur issue du z-score avec la p-valeur issue d'un simple comptage


res=[]
for i in range(1000):
    print i
    observed = np.random.normal(20,15,1)[0]
    pop =np.append(np.random.normal(10,1,1000),observed)
    observed_z=st.mstats.zscore(pop)[-1]
    observed_z_pval = 1-st.norm.cdf(observed_z)
    # empirical p-value, like
    larger_than_observed = np.sum(pop>=observed)
    emp_p_val = larger_than_observed*1.0/len(pop)
    res.append((observed,observed_z,observed_z_pval,larger_than_observed,emp_p_val,abs(emp_p_val-observed_z_pval),observed_z_pval<0.05,emp_p_val<0.05))


res.sort(key=itemgetter(5))

# the largest differences
res[-10:]

# type 1 error: situations ou le z-test rejette
[x for x in res if x[6]!=x[7]]