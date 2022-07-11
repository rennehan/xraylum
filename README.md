# xraylum

Calculate X-ray luminosity given some gas properties, from AtomDB tables. This is just a small wrapper for pyatomdb since usually I just want the X-ray luminosity of some gas particles of temperature $T$ and metallicity $Z$. 

An example:

```python
import matplotlib.pyplot as plt
import numpy as np
from xraylum.emission import get_power_grid

# 0 is Brems. emission
# H, He, C, O, Si, Fe
elements = [0, 1, 2, 6, 8, 14, 26]
all_power = get_power_grid(elements, 0.001, 100.0)

temps_to_plot = np.linspace(4, 9, 100)

total_power = np.zeros(len(temps_to_plot))
for Z in elements:
    total_power += 10**all_power[Z](temps_to_plot)

plt.figure(figsize = (7, 6), dpi = 90, facecolor = 'w')
plt.xlabel('log(Temperature) [K]')
plt.ylabel('log(Power) [erg s$^{-1}$ cm$^{3}$]')
plt.plot(temps_to_plot, np.log10(total_power))
plt.show()
plt.close()
```
