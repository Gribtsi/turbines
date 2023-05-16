#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math
import iapws
import numpy
import matplotlib.pyplot as plt
from iapws import IAPWS97 as gas
from iapws import IAPWS97
from math import sqrt
from enum import Enum
from typing import List, Tuple, Optional
from iapws import IAPWS97 as gas
from iapws import IAPWS97
import matplotlib.pyplot as plt
import numpy as np

MPa = 10 ** 6
kPa = 10 ** 3
unit = 1 / MPa
cm = 1 / 100
mm = 1 / 1000
kJ = 1000
to_kelvin = lambda x: x + 273.15 if x else None


class STM:

    delta_p0_coeff = 0.05
    delta_p_middle = 0.1
    delta_p_1 = 0.03
    def __init__(self, electrical_power, p0, t0, p_middle, t_middle, pk, t_feed_water, p_feed_water, z, internal_efficiency, mechanical_efficiency, generator_efficiency):
        self.electrical_power = electrical_power
        self.p0 = p0
        self.t0 = t0
        self.p_middle = p_middle
        self.t_middle = t_middle
        self.pk = pk
        self.t_feed_water = t_feed_water
        self.p_feed_water = p_feed_water
        self.z = z

        self.internal_efficiency = internal_efficiency
        self.mechanical_efficiency = mechanical_efficiency
        self.generator_efficiency = generator_efficiency
        
        self.delta_p0 = self.delta_p0_coeff * self.p0
        self.delta_p_middle = self.delta_p_middle * self.p_middle
        self.delta_p_1 = self.delta_p_1 * self.p_middle

        self.real_p0 = self.p0 - self.delta_p0
        self.real_p1t = self.p_middle + self.delta_p_middle
        self.real_p_middle = self.p_middle - self.delta_p_1
        
        self._point_0 = gas(P = self.p0 * unit, T = to_kelvin(self.t0))
        self.point_0 = gas(P = self.real_p0 * unit, h = self._point_0.h)
        self.point_1t = gas(P = self.real_p1t * unit, s = self._point_0.s)

        self.hp_heat_drop = (self._point_0.h - self.point_1t.h) * self.internal_efficiency
        self.h_1 = self.point_0.h - self.hp_heat_drop
        self.point_1 = gas(P = self.real_p1t * unit, h = self.h_1)
        
        self._point_middle = gas(P = self.p_middle * unit, T = to_kelvin(self.t_middle))
        self.point_middle = gas(P = self.real_p_middle * unit, h = self._point_middle.h)
        self.point_2t = gas(P = self.pk * unit, s = self._point_middle.s)

        self.lp_heat_drop = (self._point_middle.h - self.point_2t.h) * self.internal_efficiency
        self.h_2 = self.point_middle.h - self.lp_heat_drop
        self.point_2 = gas(P = self.pk * unit, h = self.h_2)

        self.efficiency_hp = (self._point_0.h - self.point_1.h) / (self._point_0.h - self.point_1t.h)
        self.efficiency_lp = (self._point_middle.h - self.point_2.h) / (self._point_middle.h - self.point_2t.h)
        
        self.point_k_water = gas(P = self.pk * unit, x=0)
        self.point_feed_water = gas(P = self.p_feed_water * unit, T = to_kelvin(self.t_feed_water))
        
        numenator_without = self.point_2.T * (self._point_middle.s - self.point_k_water.s)
        denumenator_without = (self.point_0.h - self.point_1t.h) + (self.point_middle.h - self.point_k_water.h)
        without_part = 1 - (numenator_without / denumenator_without)

        numenator_infinity = self.point_2.T * (self._point_middle.s - self.point_feed_water.s)
        denumenator_infinity = (self.point_0.h - self.point_1t.h) + (self.point_middle.h - self.point_feed_water.h)
        infinity_part = 1 - (numenator_infinity / denumenator_infinity)

        self.ksi_infinity = 1 - (without_part / infinity_part)
        
        self.coeff = (self.point_feed_water.T - self.point_2.T) / (to_kelvin(374.2) - self.point_2.T)
        
        self.ksi = self.get_ksi(self.z, self.coeff) * self.ksi_infinity
        
        eff_num = self.hp_heat_drop + self.lp_heat_drop
        eff_denum = self.hp_heat_drop + (self.point_middle.h - self.point_k_water.h)

        self.efficiency = (eff_num / eff_denum) * (1 / (1 - self.ksi))
        
        self.estimated_heat_drop = self.efficiency * ((self.point_0.h - self.point_feed_water.h) + (self.point_middle.h - self.point_1.h))
        self.inlet_mass_flow = self.electrical_power / (self.estimated_heat_drop * 1000 * self.mechanical_efficiency * self.generator_efficiency)
        
        self.condenser_mass_flow = (
        self.electrical_power /
        ((self.point_2.h - self.point_k_water.h) * 1000 * self.mechanical_efficiency * self.generator_efficiency) * ((1 / self.efficiency) - 1)
)
        
    def get_ksi(self, z, x):
        _x = numpy.array([0.6, 0.7, 0.8, 0.9, 1])
        ksi = { 8: numpy.array([0.84, 0.9, 0.94, 0.96, 0.93])}
        return numpy.interp(x, _x, ksi[z], left=None, right=None, period=None)
def legend_without_duplicate_labels(ax: plt.Axes) -> None:
            """
            Убирает дубликаты из легенды графика
            :param plt.Axes ax: AxesSubplot с отрисованными графиками
            :return None:
            """
            handles, labels = ax.get_legend_handles_labels()
            unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
            ax.legend(*zip(*unique))

def plot_process(ax: plt.Axes, points: List[IAPWS97], **kwargs) -> None:
            """
            Отрисовка процесса расширения по точкам
            :param plt.Axes ax: AxesSubplot с отрисованными графиками
            :param List[IAPWS97] points: Список инициализиованных точек процесса
            :param kwargs:
            :return None:
            """
            ax.plot([point.s for point in points], [point.h for point in points], **kwargs)

def get_isobar(point: IAPWS97) -> Tuple[List[float], List[float]]:
            """
            Собрать координаты изобары в hs осях
            :param IAPWS97 point: Точка для изобары
            :return Tuple[List[float], List[float]]:
            """
            s = point.s
            s_values = np.arange(s * 0.9, s * 1.1, 0.2 * s / 1000)
            h_values = [gas(P=point.P, s=_s).h for _s in s_values]
            return s_values, h_values

def _get_isoterm_steam(point: IAPWS97) -> Tuple[List[float], List[float]]:
            """
            Собрать координаты изотермы для пара в hs осях
            :param IAPWS97 point: Точка для изотермы
            :return Tuple[List[float], List[float]]:
            """
            t = point.T
            p = point.P
            s = point.s
            s_max = s * 1.2
            s_min = s * 0.8
            p_values = np.arange(p * 0.8, p * 1.2, 0.4 * p / 1000)
            h_values = np.array([gas(P=_p, T=t).h for _p in p_values])
            s_values = np.array([gas(P=_p, T=t).s for _p in p_values])
            mask = (s_values >= s_min) & (s_values <= s_max)
            return s_values[mask], h_values[mask]

def _get_isoterm_two_phases(point: IAPWS97) -> Tuple[List[float], List[float]]:
            """
            Собрать координаты изотермы для влажного пара в hs осях
            :param IAPWS97 point: Точка для изотермы
            :return Tuple[List[float], List[float]]:
            """
            x = point.x
            p = point.P
            x_values = np.arange(x * 0.9, min(x * 1.1, 1), (1 - x) / 1000)
            h_values = np.array([gas(P=p, x=_x).h for _x in x_values])
            s_values = np.array([gas(P=p, x=_x).s for _x in x_values])
            return s_values, h_values

def get_isoterm(point) -> Tuple[List[float], List[float]]:
            """
            Собрать координаты изотермы в hs осях
            :param IAPWS97 point: Точка для изотермы
            :return Tuple[List[float], List[float]]:
            """
            if point.phase == 'Two phases':
                return _get_isoterm_two_phases(point)
            return _get_isoterm_steam(point)

def plot_isolines(ax: plt.Axes, point: IAPWS97) -> None:
            """
            Отрисовка изобары и изотермы
            :param plt.Axes ax: AxesSubplot на котором изобразить линии
            :param IAPWS97 point: Точка для изображения изолиний
            :return None:
            """
            s_isobar, h_isobar = get_isobar(point)
            s_isoterm, h_isoterm = get_isoterm(point)
            ax.plot(s_isobar, h_isobar, color='green', label='Изобара')
            ax.plot(s_isoterm, h_isoterm, color='blue', label='Изотерма')

def plot_points(ax: plt.Axes, points: List[IAPWS97]) -> None:
            """
            Отрисовать точки на hs-диаграмме
            :param plt.Axes ax: AxesSubplot на котором изобразить точки
            :param List[IAPWS97] points: Точки для отображения
            return None
            """
            for point in points:
                ax.scatter(point.s, point.h, s=50, color="red")
                plot_isolines(ax, point)

def get_humidity_constant_line(
                point: IAPWS97,
                max_p: float,
                min_p: float,
                x: Optional[float] = None
        ) -> Tuple[List[float], List[float]]:
            """
            Собрать координаты линии с постоянной степенью сухости в hs осях
            :param IAPWS97 point: Точка для изолинии
            :param float max_p: Максимальное давление для линии
            :param float min_p: Минимальное давление для линии
            :param Optional[float] x: Степень сухости для отрисовки
            :return Tuple[List[float], List[float]]:
            """
            _x = x if x else point.x
            p_values = np.arange(min_p, max_p, (max_p - min_p) / 1000)
            h_values = np.array([gas(P=_p, x=_x).h for _p in p_values])
            s_values = np.array([gas(P=_p, x=_x).s for _p in p_values])
            return s_values, h_values

def plot_humidity_lines(ax: plt.Axes, points: List[IAPWS97]) -> None:
            """
            Отрисовать изолинии для степеней сухости на hs-диаграмме
            :param plt.Axes ax: AxesSubplot на котором изобразить изолинии
            :param List[IAPWS97] points: Точки для отображения
            return None
            """
            pressures = [point.P for point in points]
            min_pressure = min(pressures) if min(pressures) > 700 / 1e6 else 700 / 1e6
            max_pressure = max(pressures) if max(pressures) < 22 else 22
            for point in points:
                if point.phase == 'Two phases':
                    s_values, h_values = get_humidity_constant_line(point, max_pressure, min_pressure, x=1)
                    ax.plot(s_values, h_values, color="gray")
                    s_values, h_values = get_humidity_constant_line(point, max_pressure, min_pressure)
                    ax.plot(s_values, h_values, color="gray", label='Линия сухости')
                    ax.text(s_values[10], h_values[10], f'x={round(point.x, 2)}')

def plot_hs_diagram(ax: plt.Axes, points: List[IAPWS97]) -> None:
            """
            Построить изобары и изотермы для переданных точек. Если степень сухости у точки не равна 1, то построется
            дополнительно линия соответствующей степени сухости
            :param plt.Axes ax: AxesSubplot на котором изобразить изолинии
            :param List[IAPWS97] points: Точки для отображения
            return None
            """
            plot_points(ax, points)
            plot_humidity_lines(ax, points)
            ax.set_xlabel(r"S, $\frac{кДж}{кг * K}$", fontsize=14)
            ax.set_ylabel(r"h, $\frac{кДж}{кг}$", fontsize=14)
            ax.set_title("HS-диаграмма процесса расширения", fontsize=18)
            ax.legend()
            ax.grid()
            legend_without_duplicate_labels(ax)