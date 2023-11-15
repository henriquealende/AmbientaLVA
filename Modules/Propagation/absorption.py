import numpy as np

class Absorption():
    def __init__(self):
        super(Absorption, self).__init__()

    def get_ground_absorption(self, Gs, Gr, Gm, dp, hs, hr):
        """
        Calculate ground absorption for source-receptor.

        Parameters:
        - Gs: Ground factor for the source
        - Gr: Ground factor for the receptor
        - Gm: Ground factor for the middle region
        - dp: Distance between the source and the receptor
        - hs: Height of the source
        - hr: Height of the receptor

        Returns:
        - Asr: Ground absorption for source-receptor matrix
        - Am: Middle region absorption vector
        """

        # Frequency Vector
        freqVector = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

        # Ground Factors
        G = np.array([Gs, Gr, Gm])

        # Source-Receptor Heights
        h = np.array([hs, hr])

        # Ground Absorption for Source-Receptor Matrix
        Asr = np.empty([len(h), len(freqVector)])

        # Middle Region Absorption Vector
        Am = np.empty([len(freqVector)])

        # Check condition for q
        if dp <= 30 * (hs + hr):
            q = 0
        else:
            q = 1 - (30 * (hs + hr) / dp)

        # Coefficients
        for n in range(len(h)):
            a = 1.5 + 3 * (np.exp(-0.12 * ((h[n] - 5) ** 2))) * (1 - np.exp(-(dp / 50))) + 5.7 * np.exp(
                -0.09 * h[n] ** 2) * (1 - np.exp(-2.8 * 10 ** (-6) * dp ** 2))
            b = 1.5 + 8.6 * (np.exp(-0.09 * ((h[n]) ** 2))) * (1 - np.exp(-(dp / 50)))
            c = 1.5 + 14 * (np.exp(-0.46 * ((h[n]) ** 2))) * (1 - np.exp(-(dp / 50)))
            d = 1.5 + 5 * (np.exp(-0.9 * ((h[n]) ** 2))) * (1 - np.exp(-(dp / 50)))

            for m in range(len(freqVector)):
                if m == 0:
                    Asr[n, m] = -1.5
                    Am[m] = -3 * q
                elif m == 1:
                    Asr[n, m] = -1.5 + G[n] * a
                    Am[m] = -3 * q * (1 - G[2])
                elif m == 2:
                    Asr[n, m] = -1.5 + G[n] * b
                    Am[m] = -3 * q * (1 - G[2])
                elif m == 4:
                    Asr[n, m] = -1.5 + G[n] * c
                    Am[m] = -3 * q * (1 - G[2])
                elif m == 5:
                    Asr[n, m] = -1.5 + G[n] * d
                    Am[m] = -3 * q * (1 - G[2])
                else:
                    Asr[n, m] = -1.5 * (1 - G[-1])
                    Am[m] = -3 * q * (1 - G[2])
        
        Ag = np.sum([Asr[0], Asr[1], Am],axis=0)
        
        return Ag
        
    def get_barrier_absorption(self, source, r, d, h):
        """
        Calculate barrier absorption for an acoustic barrier.

        Parameters:
        - source: Type of source ('point' or 'line')
        - r: Distance from the source to the barrier
        - d: Distance from the barrier to the receptor
        - h: Height of the barrier

        Returns:
        - Abarr: Barrier absorption coefficient
        """

        # Frequency Vector
        freqVector = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

        # Speed of sound in air
        c0 = 344

        # Coefficients for different source types
        if source == 'point':
            C1 = 0.75
            C2 = 1
        elif source == 'line':
            C1 = 1
            C2 = 1

        # Distances and calculations
        A = np.sqrt((r**2) + (h**2))
        B = np.sqrt((d**2) + (h**2))
        C = r + d
        lambda_c = c0 / freqVector  # Wavelength
        N = (A + B - C) / (lambda_c / 2)

        # Barrier absorption calculation
        Abarr = 20 * C1 * np.log10((np.sqrt(2 * np.pi * N)) / (np.tanh(C2 * np.sqrt(2 * np.pi * N)))) + 5
        Abarr = np.minimum(Abarr, 20)
        # Print the result for testing purposes
        return Abarr
    

        
    def get_vegetation_absorption(self, n, rveg):
        """
        Calculate vegetation absorption coefficient.

        Parameters:
        - n: Vegetation density parameter
        - rveg: Half-width of the vegetation zone (in meters)

        Returns:
        - Aveg: Vegetation absorption coefficient
        """

        # Frequency Vector
        freqVector = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

        # Vegetation absorption calculation
        Aveg = n * ((freqVector / 1000) ** (1/3)) * (rveg/100)
        Aveg = np.minimum(Aveg, 10)
        return Aveg
    
    def get_air_absorption(self, pa, T, ht):
        """
        Calculate air absorption for different frequencies.

        Args:
        - pa (float): Atmospheric pressure (Pa).
        - T (float): Ambient temperature (°C).
        - ht (float): Relative humidity.

        Returns:
        - alpha (numpy.ndarray): Array containing air absorption values for each frequency.
        """
        # Frequency Vector
        freqVector = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])
        pr = 101.323
        T01 = 273.16
        T = T + T01
        T0 = 293.15
        C = -6.8346 * ((T01 / T) ** 1.261) + 4.6151
        psat = pr * (10 ** C)
        h = ht * (psat / pr) * ((pa / pr) ** -1)

        # Oxygen relaxation frequency
        f_r_o = ((pa / pr) * (24 + (4.04 * (10 ** 4)) * h * ((0.02 + h) / (0.391 + h))))

        # Nitrogen relaxation frequency
        f_r_n = (pa / pr) * ((T / T0) ** 0.5) * (9 + 280 * h * np.exp(-4.17 * (((T / T0) ** (-1 / 3) - 1))))

        # Calculate alpha
        alpha_classic = 8686 * (freqVector**2) * ((1.84 * (10**(-11)) * ((pa/pr)**(-1)) * ((T/T0)**(1/2))))
        
        alpha = 8686 * (freqVector**2) * ((1.84 * (10**(-11)) * ((pa/pr)**(-1)) * ((T/T0)**(1/2))) + 
                        ((T/T0)**(-5/2))* (0.01275 * ((np.exp(-2239.1/T))/(f_r_o + ((freqVector ** 2)/f_r_o))) +
                        0.1068 * ((np.exp(-3352.0/T)) / (f_r_n + ((freqVector ** 2)/f_r_n)))))
        
        return alpha