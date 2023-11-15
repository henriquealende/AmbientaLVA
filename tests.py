from Modules.Propagation.absorption import Absorption

import matplotlib.pyplot as plt
import numpy as np

absorption = Absorption()

h = [10, 20, 30, 40, 50, 60]
freqVector = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

alpha = []
for n in range(len(h)):
    alpha.append(absorption.get_air_absorption(101.325, 20, h[n]))

plt.figure(1)
plt.semilogx(freqVector, alpha[0], label = 'h = 10%')
plt.semilogx(freqVector, alpha[1], label = 'h = 20%')
plt.semilogx(freqVector, alpha[2], label = 'h = 30%')
plt.semilogx(freqVector, alpha[3], label = 'h = 40%')
plt.semilogx(freqVector, alpha[4], label = 'h = 50%')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Attenuation [dB/km]')
plt.legend()
plt.show()

