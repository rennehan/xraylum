import astropy.io.fits as fits
import numpy as np
import h5py
from scipy.interpolate import interp1d

def calculate(emissivity_files, gas_data):
    # gas_data is a dictionary with:
    # gas_data['volume']
    # gas_data['densities']
    # gas_data['temperatures']
    # gas_data['carbon']
    # gas_data['oxygen']
    # gas_data['silicon']
    # gas_data['iron']
    # gas_data['ne']
    # gas_data['nh']

    bands = np.loadtxt(bands_file)

    line_file = 'atomdb/apec_line.fits'
    coco_file = 'atomdb/apec_coco.fits'

    line = fits.open(line_file)
    continuum = fits.open(coco_file)

    # H, He, C, O, Si, Fe
    # These have to be ordered by atomic number!
    # 26 will hold everything else including iron
    atomic = [1, 2, 6, 8, 14, 26]

    min_temperature = 5.0e5

    # Solar values that Liang et al. (2016) employ
    solar = {'hydrogen': 0.7381,
             'helium': 0.2485,
             'carbon': 0.00307,
            'oxygen': 0.009618,
            'silicon': 0.0007109,
            'iron': 0.001267}



    temps = continuum[1].data['kT'] * 11604.505 * 1.0e3  # Result in K

    num_atomic = len(atomic)
    num_temps = len(temps)
    log_temps = np.log10(temps)

    def calculate_lambda(emissivities_file):
        eps_total = np.zeros((num_atomic, len(gas_data['temperature'])))
        eps_total_wide_band = np.copy(eps_total)

        with h5py.File(emissivities_file, 'r') as f:
            cont_emissivities = np.array(f['/continuum_emission'])
            line_emissivities = np.array(f['/line_emission'])
    
        log_cont = np.log10(cont_emissivities)
        log_line = np.log10(line_emissivities)

        num_interp_points = num_temps * 10
        new_temps = np.linspace(np.amin(log_temps), np.amax(log_temps), num_interp_points)

        new_cont = np.zeros((num_atomic, num_interp_points))
        new_line = np.copy(new_cont)

        # Loop over the atomic number axis
        for j in range(num_atomic):
            # do continuum emission
            non_zero_idx = np.where(cont_emissivities[j] != 0)
            if len(non_zero_idx[0]) >= 2:
                f = interp1d(log_temps[non_zero_idx], log_cont[j][non_zero_idx])

                new_temps_idx = np.where((new_temps >= np.amin(log_temps[non_zero_idx])) &
                                        (new_temps <= np.amax(log_temps[non_zero_idx])))
                new_cont[j][new_temps_idx] = f(new_temps[new_temps_idx])

            # now do line emission
            non_zero_idx = np.where(line_emissivities[j] != 0)
            if len(non_zero_idx[0]) >= 2:
                f = interp1d(log_temps[non_zero_idx], log_line[j][non_zero_idx])

                new_temps_idx = np.where((new_temps >= np.amin(log_temps[non_zero_idx])) &
                                        (new_temps <= np.amax(log_temps[non_zero_idx])))
                new_line[j][new_temps_idx] = f(new_temps[new_temps_idx])

        new_cont[np.where(new_cont == 0)] = np.nan
        new_line[np.where(new_line == 0)] = np.nan

        new_temps = 10**new_temps
        new_cont = 10**new_cont
        new_line = 10**new_line

        new_cont[np.isnan(new_cont)] = 0
        new_line[np.isnan(new_line)] = 0

        for i in range(len(new_temps)):
            if new_temps[i] < min_temperature:
                new_temps[i] = min_temperature

            print('Do temperature: %g K' % new_temps[i])

            if i == 0:
                index = np.where(gas_data['temperature'] <= new_temps[i])
            elif i == (len(new_temps) - 1):
                index = np.where(gas_data['temperature'] > new_temps[i])
            else:
                index = np.where((gas_data['temperature'] >= new_temps[i - 1]) & (gas_data['temperature'] < new_temps[i]))

            for j in range(num_atomic):
                print('Emission j=%d, eps=%g' % (j, float(new_cont[j][i] + new_line[j][i])))
                eps_total[j][index] += new_cont[j][i] + new_line[j][i]

        lambda_total = eps_total[0] + eps_total[1]
        lambda_total += (gas_data['carbon'] / solar['carbon']) * eps_total[2]
        lambda_total += (gas_data['oxygen'] / solar['oxygen']) * eps_total[3]
        lambda_total += (gas_data['silicon'] / solar['silicon']) * eps_total[4]
        lambda_total += (gas_data['iron'] / solar['iron']) * eps_total[5]

        return lambda_total

    lambda_factor = gas_data['electron_number_density'] * gas_data['hydrogen_number_density']
    lambda_factor *= gas_data['volume']

    luminosities = []
    for i in range(len(bands)):
        luminosities.append(calculate_lambda(emissivity_files[i]) * lambda_factor)

