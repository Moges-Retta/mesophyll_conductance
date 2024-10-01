# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:06:16 2021

@author: angel
"""
from enum import Enum

class Names(Enum):
    GAP_DECREASES=0
    GAP_NOCHANGE=1
    GAP_INCREASES=1.5
    INNER_CYTOSOL=1
    OUTER_CYTSOL=0
    EVERYWHERE=0.5
    PALISADE='Palisade'
    SPONGY='Spongy'
    CELLWALL='Cell wall'
    PLASMALLEMA='Plasma membrane'
    CHLENVELOPE='Chl. envelope'
    CYTOSOL='Cytosol'
    STROMA='Stroma'
    CO2=10
    LIGHT=11
    HIncana='H.Incana'
    BNigra='B.Nigra'
    BRapa='B.Rapa'
    AThaliana='A.Thaliana'
    HL='HL'
    LL='LL'
