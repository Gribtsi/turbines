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
import sys
sys.path.append('..')
from Array_Data import Array_Data
from Vector_2 import vector2

MPa = 10 ** 6
kPa = 10 ** 3
unit = 1 / MPa
cm = 1 / 100
mm = 1 / 1000
kJ = 1000
to_kelvin = lambda x: x + 273.15 if x else None


class Turbine_Stage:

    '''
    constants
    '''
    mu1_comma = 0.97
    alpha1_e = 15
    delta = 0.003
    xi_vs = 0
    
    mu_a = 0.5
    delta_a = 0.0025
    mu_r = 0.75
    z = 8
    k_friction = 0.7/1000
    m_ = 1
    k_v = 0.065
    i_ = 4
    sigma_twist_max = 20 * MPa
    sigma_stretch_max = 200 * MPa
    
    class Array:
        
        class Standart_Arrays(Enum):
            c_90_15_a = Array_Data([13, 17], [70, 120], [0.72, 0.87], [0, 0.85], 51.5 * mm, 0.45 * cm ** 3, [0.982, -0.005], [0.98, -0.008])
            r_30_21_a = Array_Data([19, 24], [25, 40], [0.58, 0.68], [0, 0.90], 25.6 * mm, 0.234 * cm ** 3, [0.965, -0.01], [0.96, -0.014])
        
        def __init__(self, profile_name, topt_choosen, e_opt, d, l, M):
            self.bounds = Turbine_Stage.Array.Standart_Arrays[profile_name].value
            self.l = l
            self.M = M
                
            self.z = math.pi * d * e_opt / self.bounds.b / topt_choosen
            if int(self.z) % 2 != 0:
                self.z = round(self.z) + round(self.z) % 2
            else:
                self.z = round(self.z) - round(self.z) % 2
                
            self.t = math.pi * d * e_opt / self.bounds.b / self.z
            self.mu = self.bounds.get_mu(l)
            
            self.placement_angle = 0
            self.ksi = 0
            self.speed_coeff = 0
            self.speed_angle = 0
            self.speed_angle_e = 0

            
        def set_specified_params(self, placement_angle, ksi, speed_angle_e, speed_angle_in):

            self.placement_angle = placement_angle
            self.speed_angle_e = speed_angle_e
            self.ksi = ksi
            
            self.speed_angle_in = speed_angle_in
            
            self.speed_coeff = math.sqrt(1 - self.ksi)
            self.speed_angle_out = math.degrees(math.asin(math.sin(math.radians(self.speed_angle_e)) * self.mu / self.speed_coeff))

        def check_for_errors(self) -> str:
            res = ''
            if max([self.speed_coeff, self.bounds.get_speed_coeff_ver(self.l)]) / min([self.speed_coeff, self.bounds.get_speed_coeff_ver(self.l)]) - 1 > 0.01:
                res += f'\nSpeed coefficient varies from verification result: \n ver: {self.bounds.get_speed_coeff_ver(self.l)}, self: {self.speed_coeff}'
            if not self.bounds.params_in_bounds(self.speed_angle_e, self.speed_angle_in, self.M) == '':
                res += ('\nSome params are out of bounds; namely' + self.bounds.params_in_bounds(self.speed_angle_e, self.speed_angle_in, self.M))
            return res

        
    '''
    G - kg/s
    p0 - MPa
    t0 - C
    H0 - kJ/kg
    d - m
    n - s^-1
    '''
    def __init__(self, G, p0, t0, c0, H0, degree_of_reactivity, d, n, a0):
        self.G = G
        
        self.point_0 = gas(P=p0 * unit, T=to_kelvin(t0))
        self.u = math.pi * d * n
        
        self.c0 = c0
        self.H0 = H0
        self.degree_of_reactivity = degree_of_reactivity
        self.d = d
        self.n = n
        self.alpha0 = a0
        
        self.H0_c = (1 - self.degree_of_reactivity) * self.H0
        self.H0_r = self.degree_of_reactivity * self.H0
        
        self.point_1t = gas(h = self.point_0.h - self.H0_c, s = self.point_0.s)
        self.c1_t = math.sqrt(2 * self.H0_c * kJ)
        
        self.F1 = self.G * self.point_1t.v / self.mu1_comma / self.c1_t
        self.el1 = self.F1 / math.pi / self.d / math.sin(math.radians(self.alpha1_e))
        self.e_opt = 5 * math.sqrt(self.el1)
    
    '''setting nozzle array'''
    def set_nozzle_array(self, array_name, topt_choosen, l, M):
        self.nozzle_array = Turbine_Stage.Array(array_name, topt_choosen, self.e_opt, self.d, l, M)
    
    '''called after filling the nozzle array props according to array atlas'''
    def calc_point1(self):
        self.c1 = self.c1_t * self.nozzle_array.speed_coeff
        self.w1 = math.sqrt(self.c1 ** 2 + self.u ** 2 - 2 * self.c1 * self.u * math.cos(math.radians(self.nozzle_array.speed_angle_out)))
        self.betta1 = math.degrees(math.atan(math.sin(math.radians(self.nozzle_array.speed_angle_out)) / (math.cos(math.radians(self.nozzle_array.speed_angle_out)) - self.u / self.c1)))
        self.alpha1 = self.nozzle_array.speed_angle_out
        
        self.delta_H_c = self.c1_t ** 2 / 2 * (1 - self.nozzle_array.speed_coeff ** 2) / kJ
        self.point1 = gas(P = self.point_1t.P, h = self.point_1t.h + self.delta_H_c)
        self.w2t = math.sqrt(2 * self.H0_r * kJ + self.w1 ** 2)
        self.point_2t = gas(h = self.point1.h - self.H0_r, s = self.point1.s)
    
    '''swtting working array'''
    def set_working_array(self, array_name, topt_choosen, l, M):
        self.working_array = Turbine_Stage.Array(array_name, topt_choosen, self.e_opt, self.d, l, M)
        self.F2 = self.G * self.point_2t.v / self.working_array.mu / self.w2t
    
    '''called after the working array props according to array atlas'''
    def calc_point_2(self):
        self.w2 = self.w2t * self.working_array.speed_coeff
        self.c2 = math.sqrt(self.w2 ** 2 + self.u ** 2 - 2 * self.w2 * self.u * math.cos(math.radians(self.working_array.speed_angle_out)))
        self.alpha2 = math.degrees(math.atan(math.sin(math.radians(self.working_array.speed_angle_out)) / (math.cos(math.radians(self.working_array.speed_angle_out)) - self.u / self.w2)))
        self.betta2 = self.working_array.speed_angle_out
        self.delta_H_r = self.w2t ** 2 / 2 * (1 - self.working_array.speed_coeff ** 2) / kJ
    
    def calc_efficiency(self) -> str:
        
        self.delta_H_vs = self.c2 ** 2 / 2 / kJ
        self.E0 = self.H0 - self.xi_vs * self.delta_H_vs
        
        self.etta_ol = (self.E0 - self.delta_H_c - self.delta_H_r - (1 - self.xi_vs) * self.delta_H_vs) / self.E0
        self.etta_ol_ver = self.u * (self.c1 * math.cos(math.radians(self.alpha1)) + self.c2 * math.cos(math.radians(self.alpha2))) / self.E0 / kJ
        if math.fabs(self.etta_ol / self.etta_ol_ver - 1) > 0.01:
            res =  f'Efficiencys are too different: \nself:{self.etta_ol}, ver:{self.etta_ol_ver}'
        else:
            res =  ''
        
        self.cf = math.sqrt(2 * self.H0 * kJ)
        self.u_cf = self.u / self.cf
        self.u_cf_opt = self.nozzle_array.speed_coeff * math.cos(math.radians(self.alpha1)) / 2 / math.sqrt(1 - self.degree_of_reactivity)
        self.d_p = self.d + self.working_array.l
        self.delta_r = 0.001 * self.d_p
        self.delta_e = (1 / ((self.delta_a * self.mu_a) ** 2) + self.z / (self.delta_r * self.mu_r) ** 2) ** (-0.5)
        self.ksi_b_y = math.pi * self.d_p * self.delta_e * self.etta_ol / self.F1
        self.delta_H_y = self.ksi_b_y * self.E0
        self.ksi_d_tr = self.k_friction * self.d ** 2 / self.F1 * self.u_cf ** 3
        self.delta_H_tr = self.ksi_d_tr * self.E0
        self.ksi_v = self.k_v * (1 - self.e_opt) * self.u_cf ** 3 * self.m_ / self.e_opt / math.sin(math.radians(self.alpha1_e))
        self.B2 = self.working_array.bounds.b * math.sin(math.radians(self.working_array.placement_angle))
        self.ksi_segm = 0.25 * self.B2 * self.working_array.l / self.F1 * self.u_cf * self.etta_ol * self.i_
        self.ksi_partial = self.ksi_segm + self.ksi_v
        self.delta_H_partial = self.E0 * self.ksi_partial
        self.H_i = self.E0 - self.delta_H_c - self.delta_H_partial - self.delta_H_r - self.delta_H_tr - self.delta_H_vs * (1 - self.xi_vs) - self.delta_H_y
        self.etta_oi = self.H_i / self.E0
        self.N_i = self.G * self.H_i * kJ

        return res
        
    def calc_durability(self):
        self.sigma_twist = (self.G * self.H0 * self.etta_ol * self.working_array.l) / (2 * self.u * self.working_array.z * self.working_array.bounds.W_min * self.e_opt)
        self.b2_new = self.working_array.bounds.b * math.sqrt(self.sigma_twist / self.sigma_twist_max)
        self.bb = 2 * math.pi * self.n
        self.sigma_stretch = 0.5 * 7800 * self.bb ** 2 * self.d * self.working_array.l
        
    def build_triangles(self):
        w1 = vector2(-self.w1 * math.cos(math.radians(self.betta1)),
                     -self.w1 * math.sin(math.radians(self.betta1)))
        c1 = vector2(-self.c1 * math.cos(math.radians(self.alpha1)),
                     -self.c1 * math.sin(math.radians(self.alpha1)))
        w2 = vector2(self.w2 * math.cos(math.radians(self.betta2)),
                     -self.w2 * math.sin(math.radians(self.betta2)))
        c2 = vector2(self.c2 * math.cos(math.radians(self.alpha2)),
                     -self.c2 * math.sin(math.radians(self.alpha2)))

        aux_vectors = [
            vector2(c1.x1,
                    c1.y1 - 40,
                    c1.x1 + self.c1 * math.cos(math.radians(self.alpha1)) +
                    self.c2 * math.cos(math.radians(self.alpha2)),
                    c1.y1 - 40),
            vector2(c1.x1, c1.y1, c1.x1, c1.y1 - 40),
            vector2(c2.x1, c2.y1, c2.x1, c1.y1 - 40),
            vector2(w1.x1,
                    w1.y1 - 20,
                    w1.x1 + self.w1 * math.cos(math.radians(self.betta1)) +
                    self.w2 * math.cos(math.radians(self.betta2)),
                    w1.y1 - 20),
            vector2(w1.x1, w1.y1, w1.x1, w1.y1 - 20),
            vector2(w2.x1, w2.y1, w2.x1, w1.y1 - 20)
        ]

        plt.style.use('_mpl-gallery')

        fig, ax = plt.subplots(figsize=(20, 10.8), layout='constrained')

        ax.plot([w1.x0, w1.x1], [w1.y0, w1.y1], linewidth=2, color="blue")
        ax.plot([w2.x0, w2.x1], [w2.y0, w2.y1], linewidth=2, color="blue")
        ax.plot([c1.x0, c1.x1], [c1.y0, c1.y1], linewidth=2, color="red")
        ax.plot([c2.x0, c2.x1], [c2.y0, c2.y1], linewidth=2, color="red")

        ax.plot([w1.x1, c1.x1], [w1.y1, c1.y1], linewidth=2, color="black")
        ax.plot([w2.x1, c2.x1], [w2.y1, c2.y1], linewidth=2, color="black")
        for i in aux_vectors:
            ax.plot([i.x0, i.x1], [i.y0, i.y1], linewidth=1, color="grey")

        ax.text(aux_vectors[0].median_point.x1, aux_vectors[0].median_point.y1 + 5, f'c1 cos(a1) + c2 cos(a2)', fontsize=14)
        ax.text(aux_vectors[3].median_point.x1, aux_vectors[3].median_point.y1 + 5, f'w1 cos(b1) + w2 cos(b2)', fontsize=14)
        ax.text(c1.x1 + 100, c1.y1 + 5, f'u', fontsize=14)
        ax.text(c2.x1 + 100, c2.y1 + 5, f'u', fontsize=14)
        ax.text(c1.median_point.x1, c1.median_point.y1 + 5, f'c1', fontsize=14)
        ax.text(c2.median_point.x1, c2.median_point.y1 + 5, f'c2', fontsize=14)
        ax.text(w1.median_point.x1, w1.median_point.y1 + 5, f'w1', fontsize=14)
        ax.text(w2.median_point.x1, w2.median_point.y1 + 5, f'w2', fontsize=14)

        ax.set_xlabel(r"x, $\frac{м}{с}$", fontsize=14)
        ax.set_ylabel(r"y, $\frac{м}{с}$", fontsize=14)
        ax.set_title("Треугольники скоростей", fontsize=18)
        ""
        # make data

        # plot
        ax.set(xlim=(-400, 250), xticks=numpy.arange(-400, 250, 50),
               ylim=(-150, 50), yticks=numpy.arange(-150, 50, 50))

        plt.show()

