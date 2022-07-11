import pyatomdb
import numpy as np
from scipy import interpolate

def _calculate_power(cie, elements, temperatures):
    res = {}
    res['power'] = {}
    res['temperature'] = []
    
    # Zero out all abundances
    cie.set_abund(np.arange(1, 31), 0.0)
    kT_list = temperatures * pyatomdb.const.KBOLTZ
    average_energies = (cie.ebins_out[1:] + cie.ebins_out[:-1]) / 2.0

    start_idx = 0
    if 0 in elements:
        start_idx = 1
    
    for i, kT in enumerate(kT_list):
        T = temperatures[i]

        res['temperature'].append(T)
        res['power'][i] = {}

        for Z in elements:                
            if Z == 0:
                # This is the electron-electron bremstrahlung component alone
                # Set all abundances to 1 (We need a full census of electrons in the plasma for e-e brems)
                cie.set_abund(elements[start_idx:], 1.0)
                # Turn on e-e bremsstrahlung
                cie.set_eebrems(True)
                spectrum = cie.return_spectrum(kT, dolines=False, docont=False, dopseudo=False)
            else:
                # This is everything else (remember it is element by element)
                # First, turn off all the elements.
                cie.set_abund(elements[start_idx:], 0.0)
                
                # Turn back on only THIS element (Z)
                cie.set_abund(Z, 1.0)
                
                # Turn off e-e bremsstrahlung (avoid double counting)
                cie.set_eebrems(False)
                spectrum = cie.return_spectrum(kT)

            # convert to keV cm3 s-1, sum
            res['power'][i][Z] = np.sum(spectrum * average_energies)

    return res

# D. Rennehan: Modified from https://atomdb.readthedocs.io/en/master/examples.html#make-cooling-curve
def get_power_grid(elements, energy_min, energy_max):
    if energy_min < 0.001:
        print('Energy below minimum 0.001keV.')
        raise ValueError
    if energy_max > 100.0:
        print('Energy above maximum 100keV.')
        raise ValueError

    # Temperatures at which to calculate curve (K)
    temperatures = np.logspace(4, 9, 51)

    # Set up the spectrum
    cie = pyatomdb.spectrum.CIESession()
    energy_bins = np.linspace(energy_min, energy_max, 10001)
    cie.set_response(energy_bins, raw=True)
    cie.set_eebrems(True)


    result = calculate_power(cie, elements, temperatures)

    # Save a 2D grid of Element, Temperature so that we can interpolate
    power_grid = np.zeros((len(elements), len(temperatures)))
    for i, Z in enumerate(elements):
        for j, T in enumerate(temperatures):
            # result['power'] is indexed by the actual element Z.
            # We reverse the indexing here, so that it is Element,Temperature indexed
            # by the ACTUAL index rather than the atomic number or Tindex.
            power_grid[i, j] = result['power'][j][Z] if result['power'][j][Z] > 0 else 1.0e-99

    # Interpolate over all temperatures. Final result will be in erg * cm**3 * s**-1.
    # We go back to the convention of having atomic number as the accessor in the
    # result dictionary.
    #
    # result['interpolated_power'][Z] contains an interp1d object that can be
    # used to interpolate between the temperatures. That should allow you to
    # compute the X-ray luminosity of a gas particle with any temperature.
    result.update({'interpolated_power': {}})
    for i, Z in enumerate(elements):
        result['interpolated_power'].update({Z: None})
        result['interpolated_power'][Z] = interpolate.interp1d(np.log10(result['temperature']), 
                                                               np.log10(power_grid[i, :] * pyatomdb.const.ERG_KEV), 
                                                               kind='linear')
        
    return result['interpolated_power']

