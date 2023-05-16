#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math
import iapws
import numpy
import matplotlib.pyplot as plt
from iapws import IAPWS97 as gas
from math import sqrt
from enum import Enum

'''stores properties of any array including boundary props.'''
class Array_Data:
        def __init__(self, speed_angle_e_bounds, speed_angle_bounds, t_bounds, M_bounds, b, W_min, mu_coeffs, speed_coeff_coeffs):
            self.speed_angle_e_bounds = speed_angle_e_bounds
            self.speed_angle_bounds = speed_angle_bounds 
            self.t_bounds = t_bounds
            self.M_bounds = M_bounds
            self.b = b
            self.W_min = W_min
            self.mu_coeffs = mu_coeffs
            self.speed_coeff_coeffs = speed_coeff_coeffs
        def get_mu(self, l):
            return self.mu_coeffs[0] + self.mu_coeffs[1] * self.b / l
        def get_speed_coeff_ver(self, l):
            return self.speed_coeff_coeffs[0] + self.speed_coeff_coeffs[1] * self.b / l
        def params_in_bounds(self, speed_angle_e, speed_angle, M):
            res = ''
            if speed_angle_e < min(self.speed_angle_e_bounds) or speed_angle_e > max(self.speed_angle_e_bounds):
                res = res + ' speed angle e '
                print(self.speed_angle_e_bounds, speed_angle_e)
            if speed_angle < min(self.speed_angle_bounds) or speed_angle > max(self.speed_angle_bounds):
                res = res + ' speed angle '
                print(self.speed_angle_bounds, speed_angle)
            if M < min(self.M_bounds) or M > max(self.M_bounds): 
                res = res + ' mach '
                print(self.M_bounds, M)
            return res

