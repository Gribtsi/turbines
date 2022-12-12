import iapws
from iapws import IAPWS97 as gas

MPa = 10 ** 6
kPa = 10 ** 3
unit = 1 / MPa
to_kelvin = lambda x: x + 273.15 if x else None

electrical_power = 250 * (10 ** 6)
p0 = 12.8 * MPa
t0 = 540

pk = 6.9 * kPa
t_feed_water = 263
p_feed_water = 1.4 * p0
z = 8

internal_efficiency = 0.85
mechanical_efficiency = 0.992
generator_efficiency = 0.99

delta_p0 = 0.05 * p0
real_p0 = p0 - delta_p0

_point_0 = gas(P=p0 * unit, T=to_kelvin(t0))
point_0 = gas(P=real_p0 * unit, h=_point_0.h)
point_kt = gas(P=pk * unit, s=_point_0.s)

hp_heat_drop = (_point_0.h - point_kt.h) * internal_efficiency
h_1 = point_0.h - hp_heat_drop
point_k = gas(P=pk * unit, h=h_1)

efficiency_hp = (_point_0.h - point_k.h) / (_point_0.h - point_kt.h)

print(efficiency_hp)

point_k_water = gas(P=pk * unit, x=0)
point_feed_water = gas(P=p_feed_water * unit, T=to_kelvin(t_feed_water))


numenator_without = point_k.T * (_point_0.s - point_k_water.s)
denumenator_without = (point_0.h - point_kt.h) + (point_0.h - point_k_water.h)
without_part = 1 - (numenator_without / denumenator_without)

numenator_infinity = point_k.T * (_point_0.s - point_feed_water.s)
denumenator_infinity = (point_0.h - point_kt.h) + (point_0.h - point_feed_water.h)
infinity_part = 1 - (numenator_infinity / denumenator_infinity)

ksi_infinity = 1 - (without_part / infinity_part)
print(ksi_infinity)



coeff = (point_feed_water.T - point_k.T) / (to_kelvin(374.2) - point_k.T)
print(coeff)

ksi = 0.9 * ksi_infinity



eff_num = hp_heat_drop
eff_denum = hp_heat_drop + (point_k.h - point_k_water.h)

efficiency = (eff_num / eff_denum) * (1 / (1 - ksi))
print(efficiency)

estimated_heat_drop = efficiency * ((point_0.h - point_feed_water.h))
print(estimated_heat_drop)


inlet_mass_flow = electrical_power / (estimated_heat_drop * 1000 * mechanical_efficiency * generator_efficiency)
print(inlet_mass_flow)

condenser_mass_flow = (electrical_power /((point_k.h - point_k_water.h) * 1000 * mechanical_efficiency *
                                          generator_efficiency) * ((1 / efficiency) - 1))
print(condenser_mass_flow)

print("Массовый расход в турбину на входе", inlet_mass_flow)
print("Массовый расход в конденсатор:", condenser_mass_flow)













