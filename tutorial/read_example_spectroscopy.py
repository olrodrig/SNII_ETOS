import numpy as np

def read_example_spectroscopy(sn):
  
  spectra_file = 'example_data/'+sn+'.spectra'
  JD_spectra   = np.atleast_1d(np.genfromtxt(spectra_file, usecols=0))
  spectra      = np.atleast_1d(np.genfromtxt(spectra_file, usecols=1, dtype=str))
  
  wl_spectra, f_spectra = [], []
  for spectrum in spectra:
      spectrum_file = 'example_data/'+spectrum
      wl_spectrum = np.genfromtxt(spectrum_file, usecols=0)
      f_spectrum  = np.genfromtxt(spectrum_file, usecols=1)
      wl_spectra.append(wl_spectrum)
      f_spectra.append(f_spectrum)
  
  print('Number of input spectra:', len(JD_spectra))
  
  return JD_spectra, wl_spectra, f_spectra

