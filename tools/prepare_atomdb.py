import pyatomdb  # Requires Python-2.7 NOT Python-3.X
import numpy as np
import argparse as ap
import h5py
import astropy.io.fits as fits

# Because pyatomdb requires Python-2.7, we can't use gallus.
# Therefore the reduction directory has to be put in manually.
parser = ap.ArgumentParser()
parser.add_argument('reduction_dir')
parser.add_argument('bands_file')
args = parser.parse_args()

# Column 0 is lower band
# Column 1 is upper band
bands = np.loadtxt(args.bands_file)

emissivity_files = []
for i in range(len(bands)):
    emissivity_files.append('%s/emissivities_%s_%s.hdf5' % (args.reduction_dir, bands[i, 0], bands[i, 1]))

line_file = 'atomdb/apec_line.fits'
coco_file = 'atomdb/apec_coco.fits'

line = fits.open(line_file)
continuum = fits.open(coco_file)

# ergperkev
erg_per_kev = 1.60205062e-9

# H, He, C, O, Si, Fe(26)
# These have to be ordered by atomic number!
atomic = [1, 2, 6, 8, 14]
atomic_in_atomdb = continuum[2].data['Z']

keV_per_kboltz = 11604.505 * 1.0e3
line_cut = line[1].data['kT'] * keV_per_kboltz
cut = continuum[1].data['kT'] * keV_per_kboltz

def calculate_continuum_emission(energy_bins):
    bin_width = energy_bins[1:] - energy_bins[:-1]
    avg_energy_in_bins = (energy_bins[1:] + energy_bins[:-1]) / 2.0

    cont_emission = np.zeros((len(atomic) + 1, len(cut)))

    for i in range(0, len(cut)):
        print('Cont. Temperature: %g K' % cut[i])
        k = 0
        for a in atomic_in_atomdb:
            spec = pyatomdb.spectrum.make_spectrum(energy_bins, i + 2, dolines = False, docont = True, dopseudo = True,
                                                   elements = [a],
                                                   linefile = line_file,
                                                   cocofile = coco_file)

            if a in atomic:
                real_idx = k
                k += 1
            else:
                real_idx = -1

            cont_emission[real_idx][i] += np.sum(spec * avg_energy_in_bins * bin_width)

    return cont_emission * erg_per_kev

def calculate_line_emission(energy_bins):
    bin_width = energy_bins[1:] - energy_bins[:-1]
    avg_energy_in_bins = (energy_bins[1:] + energy_bins[:-1]) / 2.0

    line_emission = np.zeros((len(atomic) + 1, len(line_cut)))

    for i in range(0, len(line_cut)):
        print('Line Temperature: %g K' % line_cut[i])
        k = 0
        for a in atomic_in_atomdb:
            spec = pyatomdb.spectrum.make_spectrum(energy_bins, i + 2, dolines = True, docont = False, dopseudo = True,
                                                   elements = [a],
                                                   linefile = line_file,
                                                   cocofile = coco_file)

            if a in atomic:
                real_idx = k
                k += 1
            else:
                real_idx = -1

            line_emission[real_idx][i] += np.sum(spec * avg_energy_in_bins * bin_width)

    return line_emission * erg_per_kev

for i in range(len(bands)):
    energy_bins = np.linspace(bands[i, 0], bands[i, 1], 1001)

    continuum_emission = calculate_continuum_emission(energy_bins)
    line_emission = calculate_line_emission(energy_bins)

    with h5py.File(emissivity_files[i], 'w') as f:
        h = f.create_group('Units')
        h.attrs['Emissivity'] = 'erg*cm**3/s'
        h.attrs['Band'] = '%skeV-%skeV' % (bands[i, 0], bands[i, 1])
        f.create_dataset('line_emission', data = small_band_line)
        f.create_dataset('continuum_emission', data = small_band_cont)

