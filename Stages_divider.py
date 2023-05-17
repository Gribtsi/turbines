import math
import iapws
import numpy
import numpy as np
import matplotlib.pyplot as plt
from iapws import IAPWS97 as gas
from math import sqrt
from enum import Enum
import sys
import dataclasses
from scipy.optimize import fsolve

MPa = 10 ** 6
kPa = 10 ** 3
unit = 1 / MPa
cm = 1/100
mm = 1 / 1000
kJ = 1000
to_kelvin = lambda x: x + 273.15 if x else None

class Stages_Divider:

    def __init__(self,

speed_stage_diam = 1.1,
rotation_speed = 50,
n_stages = 10,
mass_flow = 239.275,
p0 = 18.05 * MPa,
h0 = 3276.79,
pz = 3.74 * MPa,
delta_diam = 0.2,
speed_coefficient = 0.93,
alpha_1 = 14,
root_reaction_degree = 0.05,
discharge_coefficient = 0.96,
overlapping = 0.003,
efficiency = 0.88,
veernost_1=37
):
        self.root_diameter = speed_stage_diam
        self.avg_diam_1 = speed_stage_diam - delta_diam

        self.point_0 = gas(P=p0 * unit, h=h0)

        self.avg_reaction_degree_1 = self.get_reaction_degree(root_reaction_degree, veernost_1)
        self.u_cf_1 = self.get_u_cf(self.avg_reaction_degree_1, alpha_1, speed_coefficient)
        self.heat_drop_1 = self.get_heat_drop(self.avg_diam_1, self.u_cf_1, rotation_speed)

        self.h1 = self.point_0.h - self.heat_drop_1
        self.point_2 = gas(h=self.h1, s=self.point_0.s)

        upper = mass_flow * self.point_2.v * self.u_cf_1
        lower = discharge_coefficient * np.sin(np.deg2rad(alpha_1)) * rotation_speed * (np.pi * self.avg_diam_1) ** 2 * (
                    1 - self.avg_reaction_degree_1) ** 0.5

        self.blade_length_1 = upper / lower
        self.blade_length_2 = self.blade_length_1 + overlapping

        print(self.avg_diam_1 / self.blade_length_1, veernost_1)
        assert np.isclose(self.avg_diam_1 / self.blade_length_1, veernost_1, rtol=0.01)

        self.root_diameter = self.avg_diam_1 - self.blade_length_2

        self.point_zt = gas(P=pz * unit, s=self.point_0.s)
        self.full_heat_drop = h0 -  self.point_zt.h
        self.actual_heat_drop = self.full_heat_drop * efficiency
        self.hz = self.point_0.h - self.actual_heat_drop
        self.point_z = gas(P=pz * unit, h=self.hz)

        self.blade_length_z = fsolve(self.equation_to_solve, 0.01)[0]

        self.avg_diam_2 = self.root_diameter + self.blade_length_z

        self.x = np.cumsum(np.ones(n_stages) * 1 / (n_stages - 1)) - 1 / (n_stages - 1)

        self.diameters = self.linear_distribution(self.avg_diam_1, self.avg_diam_2, self.x)

        self.blade_lengths = self.linear_distribution(self.blade_length_2, self.blade_length_z, self.x)

        self.veernosts = self.diameters / self.blade_lengths

        self.reaction_degrees = self.get_reaction_degree(root_dor=root_reaction_degree, veernost=self.veernosts)

        self.u_cf = self.get_u_cf(dor=self.reaction_degrees, alpha_1=alpha_1, speed_coefficient = speed_coefficient)

        self.heat_drops = self.get_heat_drop(self.diameters, self.u_cf, rotation_speed=rotation_speed)

        self.output_speed_coeff_loss = np.full_like(self.heat_drops, 0.95)
        self.output_speed_coeff_loss[0] = 1

        self.actual_heat_drops = self.output_speed_coeff_loss * self.heat_drops

        self.mean_heat_drop = np.mean(self.actual_heat_drops)

        self.reheat_factor = 4.8 * 10 ** (-4) * (1 - efficiency) * self.full_heat_drop * (n_stages - 1) / n_stages

        print(self.full_heat_drop * (1 + self.reheat_factor) / self.mean_heat_drop)
        assert round(self.full_heat_drop * (1 + self.reheat_factor) / self.mean_heat_drop) == n_stages

        self.bias = self.full_heat_drop * (1 + self.reheat_factor) - np.sum(self.actual_heat_drops)
        self.bias = self.bias / n_stages

        self.new_actual_heat_drop = self.actual_heat_drops + self.bias

    def get_reaction_degree(self, root_dor, veernost):
        return root_dor + (1.8 / (veernost + 1.8))

    def get_u_cf(self, dor, alpha_1, speed_coefficient):
        cos = np.cos(np.deg2rad(alpha_1))
        return speed_coefficient * cos / (2 * (1 - dor) ** 0.5)

    def get_heat_drop(self, diameter, u_cf, rotation_speed):
        first = (diameter / u_cf) ** 2
        second = (rotation_speed / 50) ** 2
        return 12.3 * first * second

    def linear_distribution(self, left, right, x):
        return (right - left) * x + left

    def equation_to_solve(self, x):
        return x ** 2 + x * self.root_diameter - self.avg_diam_1 * self.blade_length_2 * self.point_z.v / self.point_2.v

