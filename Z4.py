{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "54ca4c54",
   "metadata": {},
   "source": [
    "### Домашняя работа "
   ]
  },
  {
   "cell_type": "markdown",
   "id": "370b2394",
   "metadata": {},
   "source": [
    "По аналогии с решением задач из практики, построить график зависимости $\\eta_{oi} = f(H_0)$ в диапазоне $H_0$ = (50 - 150) $\\frac {kJ}{kg}$. $\\eta_{ол}$ = 78%. $u = 160 m/s$. Все остальные переменные и условия принять такие же как на практическом заняти. </font>"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "e867a884",
   "metadata": {},
   "outputs": [],
   "source": [
    "import iapws\n",
    "from iapws import IAPWS97 as gas\n",
    "import math\n",
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "\n",
    "MPa = 10 ** 6\n",
    "kPa = 10 ** 3\n",
    "unit = 1 / MPa\n",
    "to_kelvin = lambda x: x + 273.15 if x else None"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "671dcabf",
   "metadata": {},
   "outputs": [],
   "source": [
    "p = 5 * MPa\n",
    "t = to_kelvin(489)\n",
    "\n",
    "H_0 = range(50,150)\n",
    "\n",
    "u = 160\n",
    "\n",
    "F1 = 0.025\n",
    "d_r = 1.09\n",
    "s_div_r = 0.2\n",
    "\n",
    "degree__of_reaction = 0.1\n",
    "\n",
    "z_bandage = 2\n",
    "delta_r_bandage = 1.17 / 1000\n",
    "delta_a_bandage = 4 / 1000\n",
    "\n",
    "z_rotor = 5\n",
    "d_leak_rotor = 0.36\n",
    "delta_leak_rotor = 0.4 / 1000\n",
    "\n",
    "e = 0.8\n",
    "sin_alpha_1 = 0.225\n",
    "blade_width = 0.035\n",
    "blade_length = 0.035\n",
    "blade_efficiency = 0.78\n",
    "segments = 4\n",
    "F1 = 0.025\n",
    "\n",
    "# Примем, что диафрагменное уплотнение ступенчатое, т.е. Kу = 1\n",
    "K_y = 1\n",
    "\n",
    "# По графику mu = f(зазор / толщина усика) определим при delta_leak_rotor / 0.004 ~ 0.1.\n",
    "# Возьмем квадратное уплотнение\n",
    "mu_r_rotor = 0.8\n",
    "\n",
    "# Примем коэффициент расхода сопловой решетки 0.97\n",
    "mu_nozzle = 0.97\n",
    "\n",
    "# Для надбандажного уплотнения\n",
    "mu_a = 0.5\n",
    "\n",
    "# По графику mu = f(зазор / толщина усика) определим при delta_r_bandage / 0.004 ~ 0.3.\n",
    "# Возьмем квадратное уплотнение\n",
    "mu_r = 0.8"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "b3dc9c49",
   "metadata": {},
   "outputs": [],
   "source": [
    "point = gas(P=p * unit, T=t)\n",
    "kinematic_viscosity = point.nu"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "4eb64e4e",
   "metadata": {},
   "outputs": [],
   "source": [
    "  def get_Re_numbeer(u, d_r):\n",
    "        return  u * d_r * 0.5 / kinematic_viscosity"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "65ce5225",
   "metadata": {},
   "outputs": [],
   "source": [
    "    def get_k_frictions(s_div_r, re):\n",
    "        return 2.5 * 10 ** (-2) * s_div_r ** 0.1 * re **(-0.2)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "3559c3fb",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_friction_loss_pu(s_div_r, d_r, u, kinematic_viscosity, u_div_dummy_speed, F1):\n",
    "    \n",
    "    Re_number = get_Re_numbeer(u,d_r)\n",
    "    k = get_k_frictions(s_div_r, Re_number)\n",
    "    \n",
    "    friction_loss_pu = k * d_r ** 2 * u_div_dummy_speed ** 3 / F1 \n",
    "    return friction_loss_pu"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "8c40fa44",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_ventilation_loss_pu(m, k, sin, e, u_div_dummy_speed):\n",
    "    first = k / sin\n",
    "    second = (1 - e) / e\n",
    "    third = u_div_dummy_speed ** 3\n",
    "    return first * second * third * m"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "a714830a",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_segment_loss_pu(B, l, F, u_div_dummy_speed, blade_efficiency, segments):\n",
    "    first = 0.25 * B * l / F\n",
    "    second = u_div_dummy_speed * blade_efficiency * segments\n",
    "    return first * second"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "c672a557",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_partial_losses_pu(u_div_dummy_speed, blade_width, blade_length, F1, blade_efficiency, segments):\n",
    "    \n",
    "    ventilation_loss_pu = get_ventilation_loss_pu(\n",
    "    m = 1,\n",
    "    k = 0.065,\n",
    "    e = e,\n",
    "    u_div_dummy_speed = u_div_dummy_speed,\n",
    "    sin = sin_alpha_1\n",
    "    )\n",
    "    \n",
    "    segment_loss_pu = get_segment_loss_pu(\n",
    "    B = blade_width,\n",
    "    l = blade_length,\n",
    "    F = F1,\n",
    "    u_div_dummy_speed = u_div_dummy_speed,\n",
    "    blade_efficiency = blade_efficiency,\n",
    "    segments = segments\n",
    "    )\n",
    "    \n",
    "    partial_losses_pu = segment_loss_pu + ventilation_loss_pu\n",
    "    return partial_losses_pu\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "9e709a2a",
   "metadata": {},
   "outputs": [],
   "source": [
    "def compute_equal_gap(z, delta_r, mu_r_rotor, delta_a, mu_a):\n",
    "    first = 1 / (mu_a * delta_a) ** 2\n",
    "    second = z / (mu_r_rotor * delta_r) ** 2\n",
    "    return (first + second) ** (-0.5)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "5361311b",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_bandage_leak_loss_pu(d_shroud, delta_eq, F, dor, l, efficiency):\n",
    "    d_avg = d_shroud - l\n",
    "    first = math.pi * d_shroud * delta_eq / F\n",
    "    second = dor + 1.8 * (l / d_avg)\n",
    "    return first * (second) ** 0.5 * efficiency"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "7e435c93",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_disk_leak_loss_pu(K, F, mu_r_rotor, mu_nozzle, F_nozzle, z, efficiency):\n",
    "    upper = mu_r_rotor * K * F * efficiency\n",
    "    lower = mu_nozzle * F_nozzle * z ** 0.5\n",
    "    return upper / lower"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "efd0d581",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_leak_losses_pu(z_bandage, delta_r_bandage, mu_r_rotor, delta_a_bandage, z_rotor, mu_a, mu_nozzle, F1, degree__of_reaction, blade_length, delta_leak_rotor, d_leak_rotor, blade_efficiency, K_y):\n",
    "    \n",
    "    d_shroud = delta_r_bandage / 0.001 # По определению задачи\n",
    "    delta_eq_bandage = compute_equal_gap(\n",
    "        z = z_bandage,\n",
    "        delta_r = delta_r_bandage,\n",
    "        mu_r_rotor = mu_r_rotor,\n",
    "        delta_a = delta_a_bandage,\n",
    "        mu_a = mu_a\n",
    "    )\n",
    "\n",
    "    bandage_leak_loss_pu = get_bandage_leak_loss_pu(\n",
    "        d_shroud = d_shroud,\n",
    "        delta_eq = delta_eq_bandage,\n",
    "        F = F1,\n",
    "        dor = degree__of_reaction,\n",
    "        l = blade_length,\n",
    "        efficiency = blade_efficiency\n",
    "        )\n",
    "        \n",
    "    F_leak_rotor = math.pi * d_leak_rotor * delta_leak_rotor\n",
    "    \n",
    "    disk_leak_loss_pu = get_disk_leak_loss_pu(\n",
    "        K = K_y,\n",
    "        F = F_leak_rotor,\n",
    "        mu_r_rotor = mu_r_rotor,\n",
    "        mu_nozzle = mu_nozzle,\n",
    "        F_nozzle = F1,\n",
    "        z = z_rotor,\n",
    "        efficiency = blade_efficiency\n",
    "    )\n",
    "    \n",
    "    leak_losses_pu = disk_leak_loss_pu + bandage_leak_loss_pu\n",
    "    return leak_losses_pu"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "3ff64f86",
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_internal_efficiency(blade_efficiency, friction_loss_pu, partial_losses_pu, leak_losses_pu):\n",
    "    \n",
    "    internal_efficiency = blade_efficiency - friction_loss_pu - partial_losses_pu - leak_losses_pu\n",
    "    \n",
    "    return internal_efficiency"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "65cb66f2",
   "metadata": {},
   "outputs": [],
   "source": [
    "leak_losses_pu = get_leak_losses_pu(\n",
    "    z_bandage = z_bandage,\n",
    "    delta_r_bandage = delta_r_bandage,\n",
    "    mu_r_rotor = mu_r_rotor,\n",
    "    delta_a_bandage = delta_a_bandage,\n",
    "    z_rotor = z_rotor,\n",
    "    mu_a = mu_a,\n",
    "    mu_nozzle = mu_nozzle,\n",
    "    F1 = F1,\n",
    "    degree__of_reaction = degree__of_reaction,\n",
    "    blade_length = blade_length,\n",
    "    delta_leak_rotor = delta_leak_rotor,\n",
    "    d_leak_rotor  = d_leak_rotor,\n",
    "    blade_efficiency = blade_efficiency,\n",
    "    K_y = K_y\n",
    ")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "a48abaa8",
   "metadata": {},
   "outputs": [],
   "source": [
    "internal_efficiency = []\n",
    "for H_0value in H_0:\n",
    "    H_0value=H_0value*1000\n",
    "    dummy_speed = (2 * H_0value) ** 0.5\n",
    "    u_div_dummy_speed = u / dummy_speed\n",
    "    \n",
    "    friction_loss_pu = get_friction_loss_pu(\n",
    "        s_div_r = s_div_r,\n",
    "        d_r = d_r, \n",
    "        u = u,\n",
    "        kinematic_viscosity = kinematic_viscosity,\n",
    "        u_div_dummy_speed = u_div_dummy_speed,\n",
    "        F1 = F1\n",
    "    )\n",
    "    \n",
    "    partial_losses_pu = get_partial_losses_pu(\n",
    "        u_div_dummy_speed = u_div_dummy_speed,\n",
    "        blade_width = blade_width,\n",
    "        blade_length = blade_length,\n",
    "        F1 = F1,\n",
    "        blade_efficiency = blade_efficiency,\n",
    "        segments = segments\n",
    "    )\n",
    "    internal_efficiency.append(get_internal_efficiency(blade_efficiency,friction_loss_pu,partial_losses_pu,leak_losses_pu))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "e0b05e83",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAosAAAHrCAYAAACn9tfQAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjYuMCwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy89olMNAAAACXBIWXMAAA9hAAAPYQGoP6dpAABu20lEQVR4nO3deVhUZf8G8HsGGDZZBGQVFTdURFAJRNMwFVPLqNw10EzrffNXyVup9br2lq1mi0kLmpWGmWamZiBqaiIoiIgLggooyC67DMPM+f1BTE7DsCh4Brg/18V1Oc95zpnvnGfQ27M8RyIIggAiIiIionpIxS6AiIiIiPQXwyIRERER6cSwSEREREQ6MSwSERERkU4Mi0RERESkE8MiEREREenEsEhEREREOjEsEhEREZFODItEREREpBPDIhERERHpxLBIRERERDoxLBKJJCwsDOPHj4eDgwOMjIzg6OiIhx56CN9++y1UKpXY5RFRK/rmm28gkUhw+vTpepcHBARg4MCBzd6uXC7HkiVL4OzsDFNTU/j5+SEqKupey6UOzlDsAog6qi1btsDJyQnLly+HpaUliouLcfLkScydOxe//fYbfvjhB7FLJKI2Zu7cufjpp5/w8ssvo0+fPvjmm28wceJEHD58GA8++KDY5VEbxbBIJJKjR4/CyMhIo+3FF1+Era0tPvvsM6xduxY9evQQpzgianPi4uIQERGB999/H6+88goAIDg4GAMHDsRrr72GEydOiFwhtVU8DU0kkn8GxTp1AVEqrf31zMjIwL///W+4u7vD1NQUtra2mDp1KtLT0zXWW7VqFSQSifrHwsICvr6+2L17d7397lReXg5HR0dIJBIcOXJEY1lWVhbmz58PZ2dnGBsbw83NDf/6179QXV3d7O3V9bW3t4dCodBY54cfflDXXlBQoG4/c+YMJkyYAEtLS3Tq1AljxozByZMn6913DdX6z/1T38+RI0fU/e6sobm2bNkCiUSCM2fO4MUXX4SjoyNMTU3x2GOP4datW3e93X9qbN809TPrsnHjRnh5ecHKygrm5ubw8vJCeHi4Rh9d++v06dOQSCT45ptv1G1N/S4HBAQgICBAoy09PV1re0DtmD/zzDNwcHCAsbExPDw8sGnTpruusbm/H4cPH8bIkSPRuXNnjf26aNEi3G8//fQTDAwMsHDhQnWbiYkJ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KlUfMW1urU0xatQo2NnZAQAmT56MP//8EwcOHGhyWPxnUDl37pxGe01NDVavXo3w8HBUV1cjODgYH374ISQSCT7++GMkJSVp3Xlcpzn75l6ZmZmpv1NPPPEE0tPT8eabb2qFxTv3FwBYW1trbaup32UTExPs27cPly9fxvXr1yEIAnJzczFnzhx1ny5dusDCwgJKpVLjO9+QptRYp6nf51deeQWpqanYuXMnvv32W8hksnr/Y3Y/eHt74/DhwygtLdUIxbGxserlRHeD1ywSiWT//v31tn/11VeQSCTq05UGBgYap+sA4NNPP4VSqWz0PZRKJaqrq+s9NbV06VKMGDECkydPrnddqVSKoKAg/Prrr/U+ZeKfNTW2vX+ysbHB448/jqSkJDzzzDNayw0MDBAYGIhffvlFY2qV3NxcbNu2DQ8++KD6H8Tm1tpcgiBAEIQmheBz586ha9eu6Ny5s0Z7UlISJBKJ+tTwf/7zH5w/fx7nz59HamoqDh48iB07dqj7NnTNWnP2TUtSKpW4detWvd+npmjud7lv374YM2YMxo4dixEjRmht66mnnsLOnTvr/U9Cfn7+XdVYp6nf519//RVffvklvv76a0ycOLHJwbU1TJkyBUqlUuO0vFwux+bNm+Hn58c7oemu8cgikUhmzZqFfv364YknnoCDgwPy8/Px22+/4fDhw3jjjTfUpysfffRRfPfdd7CyssKAAQMQExODgwcP6ryo/vvvvwdQexp69+7dSE9Pr3cKkcjISPz5558N1vj2228jMjISDz30EBYuXIj+/fvj5s2b2LFjB44fP65xZKYp2/unb775Bhs2bNA42nOn//3vf4iKisKDDz6If//73zA0NMQXX3wBuVyO9957765rbYpDhw5pnIZOS0tr0lQsuoLeuXPn0LNnT5ibm+PGjRv49ttvkZ6eDisrKwDAhAkTEB8fj2nTpiEpKUnjKFp9mrNv7taoUaMQEBCAbt26oby8HD/99BPOnDmDDz744K6219zvcmPeeecdHD58GH5+fliwYAEGDBiAoqIiJCQk4ODBgygqKrqr7QJN+z7n5ORg/vz5ePbZZxEUFHTX79VS/Pz8MHXqVCxbtgx5eXno3bs3tmzZgvT0dJ1HqYmagmGRSCTvvPMOfv31V3zyySfIy8tDp06d4Ofnh/3792PChAnqfh9//DEMDAyw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cWlZjf2u/FN9YZFj0nQSQbjjEQ1ERERERHfgNYtEREREpBPDIhERERHpxLBIRERERDoxLBIRERGRTgyLRERERKQTwyIRERER6cSwSEREREQ6MSwSERERkU4Mi0RERESkE8MiEZGI5s6di6CgIK32I0eOQCKRoLi4uNFtJCUlYeTIkTAxMYGrqyvee++9li+UiDoshkUiojastLQUgYGB6N69O+Lj4/H+++9j1apV+PLLL8UujYjaCUOxCyAioru3detWVFdXY9OmTZDJZPDw8EBiYiLWrVuHhQsXil0eEbUDPLJIRNSGxcTEYNSoUZDJZOq28ePHIyUlBbdu3RKxMiJqL3hkkYhIZHv37kWnTp002pRKZZPWzcnJgZubm0abg4ODelnnzp1bpkgi6rAYFomIRDZ69Ghs3LhRoy02NhZz5swRqSIior8xLBIRiczc3By9e/fWaLtx40aT1nV0dERubq5GW91rR0fHlimQiDo0XrNIRNSG+fv74+jRo1AoFOq2qKgouLu78xQ0EbUIhkUiojZs1qxZkMlkmD9/Ps6fP4/t27fj448/RmhoqNilEVE7wdPQRERtmJWVFSIjI/HCCy9g6NChsLOzw4oVKzhtDhG1GIkgCILYRRARERGRfuJpaCIiIiLSiWGRiEiPTZgwAZ06dar35+233xa7PCLqAHgamohIj2VlZeH27dv1LrOxsYGNjc19roiIOhqGRSIiIiLSiaehiYiIiEgnhkUiIiIi0olhkYiIiIh0YlgkIiIiIp0YFomIiIhIJ4ZFIiIiItKJYZGIiIiIdPp/j5t5KU7T9rkAAAAASUVORK5CYII=\n",
      "text/plain": [
       "<Figure size 640x480 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "x = H_0\n",
    "y = internal_efficiency\n",
    "plt.figure(layout = 'constrained')\n",
    "plt.plot(x,y)  \n",
    "plt.xlabel('H_0')\n",
    "plt.ylabel('$\\eta_{oi}$')\n",
    "plt.title(\"Зависимость $\\eta_{oi}$ от значения H_0\")\n",
    "plt.grid()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0d846148",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.13"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
