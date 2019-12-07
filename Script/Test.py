# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 10:15:25 2019

@author: emyli
"""

import subprocess
import os

x = 0
for x in range(30):
    subprocess.call(["gcc", "main.c"])
    tmp=subprocess.call("DA_495.exe")
    print (tmp)
    os.rename(r'C:\Users\emyli\Desktop\DA_495\cmake-build-debug\DA.csv',r'C:\Users\emyli\Desktop\DA_495\cmake-build-debug\DA'+str(x)+'.csv')