import matplotlib.pyplot as plt
from numpy import arange
from weight import get_fuel_moment

arms_fuel = []
masses_fuel = arange(2000, 0, -20)

for m in masses_fuel:
    arms_fuel.append(get_fuel_moment(m))

x = [min(masses_fuel), max(masses_fuel)]
y = [min(arms_fuel), max(arms_fuel)]

fig = plt.figure(figsize=(9, 4), dpi=100)

plt.plot(x, y, color='k', linestyle='dashed', linewidth=1.2)

plt.scatter(masses_fuel, arms_fuel, color='orange', marker='.', s=35)

axes = plt.gca()
axes.set_xlim(min(masses_fuel), max(masses_fuel))
axes.set_ylim(min(arms_fuel), max(arms_fuel))

plt.ylabel('Moment [kg m]', fontsize=10)
plt.xlabel('Fuel Mass in Tanks [kg]', fontsize=10)
plt.minorticks_on()
plt.grid(b=True, which='major', color='grey', linestyle='-')
plt.grid(b=True, which='minor', color='lightgrey', linestyle='dotted')

plt.savefig('../plots/fuel_moment.png', dpi=600, orientation='portrait')

plt.show()
