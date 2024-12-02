import mlspectra

spec_path = 'ML/Spectra'

max_wavelength = 800
min_wavelength = 200

N_bin = 24
read_P = False            # If true, binned spectra will be loaded from 'file_spec', if false will be calculated and stored in this file
resolution = int((max_wavelength-min_wavelength)/N_bin)
file_spec = str(resolution)+'nmResolution_spectra.dat'
file_geom = 'fingerprint.dat'
N_fing = 40
#Int_lam, lambda_min, dlambda = mlspectra.bin_spectra_uniform(spec_path, read_P, file_spec, file_geom, e_min, e_max, N_bin, N_fing, 'qc_Opt.dat', 'dihedral_b3lyp_Opt.dat')

Int_lam, lambda_min, dlambda = mlspectra.bin_spectra_waveuniform(spec_path, read_P, file_spec, file_geom, 200, 800, N_bin, N_fing, 'qc_Opt.dat', 'dihedral_b3lyp_Opt.dat')
