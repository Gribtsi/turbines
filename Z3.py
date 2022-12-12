import math
import iapws
from iapws import IAPWS97 as gas
import matplotlib.pyplot as plt

MPa = 10 ** 6
kPa = 10 ** 3
unit = 1 / MPa
to_kelvin = lambda x: x + 273.15 if x else None





def get_r_list(d_avg, d_root, count):
    r_root = d_root/2
    r_tail = (2*d_avg-d_root)/2
    step = (r_tail-r_root)/(count-1)
    return [r_root+step*i for i in range(count)]

def calculate_degree_of_reaction(degree_of_reaction_avg, r,r_avg):
    return (degree_of_reaction_avg-1)*(1/((r/r_avg)**1.7)+1/(degree_of_reaction_avg-1))

def calculate_alpha1(alpha_1_avg, r,r_avg):
    return (r/r_avg)**2*alpha_1_avg

def get_triangle_pair(p_0, t_0, p_2, avg_diameter, degree_of_reaction, alpha_1_deg,delta_beta_deg, fi,psi,rotation_speed):

    triangle_pair = []
    inlet_point = gas(P=p_0 * unit, T=t_0)
    outlet_point = gas(P=p_2 * unit, s=inlet_point.s)

    theoretical_heat_drop = inlet_point.h - outlet_point.h
    stator_heat_drop = theoretical_heat_drop * (1 - degree_of_reaction)
    rotor_heat_drop = theoretical_heat_drop * degree_of_reaction

    c_1t = (2 * 1000 * stator_heat_drop) ** 0.5
    c_1 = c_1t * fi

    u = math.pi * avg_diameter * rotation_speed

    sin_alpha_1 = math.sin(math.radians(alpha_1_deg))
    cos_alpha_1 = math.cos(math.radians(alpha_1_deg))

    w_1 = (c_1 ** 2 + u ** 2 - 2 * c_1 * u * cos_alpha_1) ** 0.5

    w_2t = (w_1 ** 2 + 2 * rotor_heat_drop * 1000) ** 0.5
    w_2 = w_2t * psi

    beta_1 = math.atan(sin_alpha_1 / (cos_alpha_1 - u / c_1))
    beta_1_deg = math.degrees(beta_1)
    beta_2_deg = beta_1_deg - delta_beta_deg

    sin_beta_2 = math.sin(math.radians(beta_2_deg))
    cos_beta_2 = math.cos(math.radians(beta_2_deg))

    c_2 = (w_2 ** 2 + u ** 2 - 2 * w_2 * u * cos_beta_2) ** 0.5

    alpha_2 = math.atan(sin_beta_2 / (cos_beta_2 - u / w_2))
    alpha_2_deg = math.degrees(alpha_2)

    c1_plot = [[0, -c_1 * cos_alpha_1], [0, -c_1 * sin_alpha_1]]
    u1_plot = [[-c_1 * cos_alpha_1, -c_1 * cos_alpha_1 + u], [-c_1 * sin_alpha_1, -c_1 * sin_alpha_1]]
    w1_plot = [[0, -c_1 * cos_alpha_1 + u], [0, -c_1 * sin_alpha_1]]
    w2_plot = [[0, w_2 * cos_beta_2], [0, -w_2 * sin_beta_2]]
    u2_plot = [[w_2 * cos_beta_2, w_2 * cos_beta_2 - u], [-w_2 * sin_beta_2, -w_2 * sin_beta_2]]
    c2_plot = [[0, w_2 * cos_beta_2 - u], [0, -w_2 * sin_beta_2]]

    triangle_pair.append(c1_plot)
    triangle_pair.append(u1_plot)
    triangle_pair.append(w1_plot)
    triangle_pair.append(w2_plot)
    triangle_pair.append(u2_plot)
    triangle_pair.append(c2_plot)

    return triangle_pair


p_0 = 16.7 * MPa
t_0 = to_kelvin(520)
p_2 = 14.5 * MPa

d_avg = 0.892
d_root = 0.8
degree_of_reaction_avg = 0.2
alpha_1_avg = 13

delta_beta_deg = 5
fi = 0.97
psi = 0.935
rotation_speed = 50

k = 10

tringle_list = [get_triangle_pair(p_0, t_0, p_2, i*2,
                                  calculate_degree_of_reaction(degree_of_reaction_avg, i,d_avg/2),
                                  calculate_alpha1(alpha_1_avg, i,d_avg/2),
                                  delta_beta_deg,fi,psi,rotation_speed
                                  ) for i in get_r_list(d_avg,d_root,k)]
plot = plt
fig, ax = plot.subplots(1, 1, figsize=(15, 5))

i = tringle_list[0]

ax.plot(i[0][0], i[0][1], label='C_1', c='red')
ax.plot(i[1][0], i[1][1], label='u_1', c='blue')
ax.plot(i[2][0], i[2][1], label='W_1', c='green')

ax.plot(i[3][0], i[3][1], label='W_2', c='green')
ax.plot(i[4][0], i[4][1], label='u_2', c='blue')
ax.plot(i[5][0], i[5][1], label='C_2', c='red')

tringle_list.pop(0)

for i in tringle_list:

    ax.plot(i[0][0], i[0][1], c='red')
    ax.plot(i[1][0], i[1][1], c='blue')
    ax.plot(i[2][0], i[2][1], c='green')

    ax.plot(i[3][0], i[3][1], c='green')
    ax.plot(i[4][0], i[4][1], c='blue')
    ax.plot(i[5][0], i[5][1], c='red')

ax.set_title("Треугольник скоростей")
ax.legend()

plot.show()