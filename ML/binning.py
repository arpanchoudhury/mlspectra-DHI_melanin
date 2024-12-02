import mlspectra

spec_path = 'ML/Spectra'

max_wavelength = 800
min_wavelength = 200

resolution = 25
N_bin = int((max_wavelength-min_wavelength)/N_bin)
read_P = False            # If true, binned spectra will be loaded from 'file_spec', if false will be calculated and stored in this file
file_spec = str(resolution)+'nmResolution_spectra.dat'
file_geom = 'fingerprint_descrptor.dat'
N_fing = 40

Int_lam, lambda_min, dlambda = mlspectra.bin_spectra_waveuniform(spec_path, read_P, file_spec, file_geom, min_wavelength, max_wavelength, N_bin, N_fing, 'qc_Opt.dat', 'dihedral_b3lyp_Opt.dat')
