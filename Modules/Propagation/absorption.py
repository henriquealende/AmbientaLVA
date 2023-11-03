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
        freqVector = np.array([63, 125, 250, 500, 1000, 2000, 4000, 8000])

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
        freqVector = np.array([63, 125, 250, 500, 1000, 2000, 4000, 8000])

        # Speed of sound in air
        c0 = 344

        # Coefficients for different source types
        if source == 'point':
            C1 = 0.75
            C2 = 1
        else:
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

        # Print the result for testing purposes
        print(Abarr)

        return Abarr