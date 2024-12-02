import mlspectra

spec_path = '/home/debashree/arpan/melanin/dhi-melanin/tetramer_opt/prop'

#e_min = 1.5498
#e_max = 6.1992

N_bin = 60
read_P = False            # If true, binned spectra will be loaded from 'file_spec', if false will be calculated and stored in this file
file_spec = 'spectra_60bins_waveuniform.dat'
file_geom = 'fingerprint_40.dat'
N_fing = 40
#Int_lam, lambda_min, dlambda = mlspectra.bin_spectra_uniform(spec_path, read_P, file_spec, file_geom, e_min, e_max, N_bin, N_fing, 'qc_Opt.dat', 'dihedral_b3lyp_Opt.dat')

Int_lam, lambda_min, dlambda = mlspectra.bin_spectra_waveuniform(spec_path, read_P, file_spec, file_geom, 200, 800, N_bin, N_fing, 'qc_Opt.dat', 'dihedral_b3lyp_Opt.dat')
