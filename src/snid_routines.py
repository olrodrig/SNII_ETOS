import os
import numpy as np

#select a region of the spectrum within [wmin, wmax]
def spectrum_in_range(l_original, f_original, wmin, wmax):
    
  l_trim, f_trim = [], []
  for l, f in zip(l_original, f_original):
      if l >= wmin and l <= wmax: 
          l_trim.append(l)
          f_trim.append(f)
  l_trim, f_trim = np.array(l_trim), np.array(f_trim)

  return l_trim, f_trim

#run SNID
def run_snid(l, f, snid_options, snid_templates):
    
  #save the input spectrum in a file
  input_spectrum = 'input.txt'
  file_spectrum  = open(input_spectrum, 'w')
  for ll, ff in zip(l, f):  file_spectrum.write('%8.3f  %.7e\n' % (ll, ff))
  file_spectrum.close()
    
  #run SNID
  os.system('snid '+snid_options+' '+snid_templates+' '+input_spectrum)
  os.system('rm '+input_spectrum)
  os.system('rm snid.param')
  
  snid_file = 'input_snid.output'
  if os.path.exists(snid_file) == False:
      print('ERROR in run_snid: there is no SNID output')
      sne_bm, rlaps, zs, snid_phases = np.array([]), np.array([]), np.array([]), np.array([])
  else:
      #read important parameters
      sne_bm      = np.atleast_1d(np.genfromtxt(snid_file, skip_header=67, usecols=1, dtype=str))
      rlaps       = np.atleast_1d(np.genfromtxt(snid_file, skip_header=67, usecols=4))
      zs          = np.atleast_1d(np.genfromtxt(snid_file, skip_header=67, usecols=5))
      snid_phases = np.atleast_1d(np.genfromtxt(snid_file, skip_header=67, usecols=7))
      os.system('rm '+snid_file)
  return sne_bm, rlaps, zs, snid_phases 

#average snid results
def average_snid_results(phases, rms_phases):
    
  phase, err_phase = float('nan'), float('nan')
  if len(phases) != 0:
      
      if len(phases) == 1:
          phase, err_phase = phases[0], rms_phases[0]
      else:
          #compute the phase and its error minimizing the log-likelihood
          residuals = phases - np.mean(phases)
          rms       = np.sqrt(np.sum(residuals**2)/float(len(residuals)-1))
          err_0s  = np.linspace(0.0, rms, 100)
          m2ln_L_min = 1.0e90
          for err_0 in err_0s:
              var    = rms_phases**2 + err_0**2
              phase  = np.sum(phases/var)/np.sum(1.0/var)
              m2ln_L = np.sum(np.log(var)+(phases-phase)**2/var)
              if m2ln_L < m2ln_L_min:
                  m2ln_L_min = m2ln_L
                  phase_min  = phase
                  err_0_min  = err_0
          err_0 = err_0_min
          phase = phase_min
          
          var       = rms_phases**2 + err_0**2
          err_phase = np.sqrt(np.sum(rms_phases**2/var**2))/np.sum(1.0/var)
          err_phase = np.sqrt(err_phase**2)
              
  return phase, err_phase

#select results based on rlap and zfilter values
def filter_snid_results(t_specs, sne_match, rlaps, zs, t0s, rlap_min, zfilter, N_bests):

  r, p, s, t, zz, i = [], [], [], [], [], 0
  for rlap, z, t0, t_spec, sn_match in zip(rlaps, zs, t0s, t_specs, sne_match): 
      if rlap >= rlap_min and abs(z)<zfilter:
          if i<N_bests:
              r.append(rlap)
              zz.append(z)
              p.append(t0)
              t.append(t_spec)
              s.append(sn_match)
              i = i + 1             

  rlaps, zs, t0s, t_specs, sne_match = np.array(r), np.array(zz), np.array(p), np.array(t), s
  
  return t_specs, sne_match, rlaps, zs, t0s

#compute snid phases
def compute_snid_time(t_specs, sne_match, rlaps, zs, t0s, N_bests, JD_spectrum=-9999):

  zfilter      = 0.01
  rlap_min     = 5.0
  rlap2err_par = np.array([3.0, 72.0, -75.0]) #rlap to rms(t0) conversion

  #select good results
  t_specs, sne_match, rlaps, zs, t0s = filter_snid_results(t_specs, sne_match, rlaps, zs, t0s, rlap_min, zfilter, N_bests)
  
  if len(t0s) != 0:
      #convert rlaps to rms(t0)
      rms_t0s = 0.0
      for i in range(0, len(rlap2err_par)):  rms_t0s = rms_t0s + rlap2err_par[i]/rlaps**i
  
      #average results
      t0, err_t0 = average_snid_results(t0s, rms_t0s)
  else:
      if JD_spectrum == -9999:
          print('ERROR in compute_snid_time: no match meeting our SNID requirements')
      else:
          print("spectrum at JD="+str(JD_spectrum)+" does not have best-matches meeting our SNID requirements")
      t0, err_t0 = float('nan'), float('nan')
      t_specs, sne_match, rlaps, zs, t0s, rms_t0s = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
  
  return t_specs, sne_match, rlaps, zs, t0s, rms_t0s, t0, err_t0

#compute t0 and rms(t0) with SNID and the input spectra
def snid_results(templates_path, sn, z, JD_fd, JD_spectra, wl_spectra, f_spectra, snid_verbose, display_snid_plots):
  
  wmin, wmax = 4100.0, 7000.0
  N_bests = 10
  
  if snid_verbose:
      snid_verbose = '1'
  else:
      snid_verbose = '0'
      
  if display_snid_plots:
      display_snid_plots = '1'
  else:
      display_snid_plots = '0'
  
  #aplicaremos SNID a todos los espectros input
  t_specs_all, sne_match_all, rlaps_all, zs_all, t0s_all, N_good_spectra = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), 0
  JD_list, t_specs_list, sne_match_list, rlaps_list, zs_list, t0s_list = [], {}, {}, {}, {}, {}
  for JD_spectrum, wl_spectrum, f_spectrum in zip(JD_spectra, wl_spectra, f_spectra):

      #select data within a wl range
      wl, f = spectrum_in_range(wl_spectrum, f_spectrum, wmin, wmax)

      #define the SNID options                
      snid_options = 'verbose='+snid_verbose+' plot='+display_snid_plots+' forcez=0.0 rlapmin=1.0 k1=1 k2=4 k3=85 k4=102'
      
      #if the input SN is among the templates, then we avoid their templates
      sn_templates = np.genfromtxt(templates_path+'/templist', usecols=0, dtype=str)
      if sn+'.lnw' in sn_templates:  snid_options = snid_options + ' avoid='+sn
      
      snid_templates = 'tempdir='+templates_path+'/'
      
      #run snid
      sne_match, rlaps, zs, snid_phases = run_snid(wl, f, snid_options, snid_templates)

      if len(sne_match) != 0:
          #convert snid phases to t0:  JD_spec-JD_t0=phase -> t0_JD = JD_spec-phase -> t0 = JD_spec-phase - JD_fd = (JD_spec-JD_fd)/(1.0+z) - phase
          t0s     = (JD_spectrum-JD_fd)/(1.0+z)-snid_phases
          t_specs = snid_phases-snid_phases+JD_spectrum
          
          #compute preliminar result
          t_spec_best, sne_best, rlap_best, zs_best, t0_best, err_t0_best, t0, err_t0 = compute_snid_time(t_specs, sne_match, rlaps, zs, t0s, N_bests, JD_spectrum=JD_spectrum)
          
          if np.isnan(t0) == False:            
              JD_list.append(JD_spectrum)
              t_specs_list[str(JD_spectrum)]   = t_spec_best
              sne_match_list[str(JD_spectrum)] = sne_best
              rlaps_list[str(JD_spectrum)]     = rlap_best
              zs_list[str(JD_spectrum)]        = zs_best
              t0s_list[str(JD_spectrum)]       = t0_best
              
              N_good_spectra = N_good_spectra+1
              t_specs_all, sne_match_all = np.append(t_specs_all, t_spec_best), np.append(sne_match_all, sne_best)
              rlaps_all, zs_all, t0s_all = np.append(rlaps_all, rlap_best), np.append(zs_all, zs_best), np.append(t0s_all, t0_best)

  #sort by rlap
  inds       = rlaps_all.argsort()[::-1]
  t_specs_all, sne_match_all = t_specs_all[inds], sne_match_all[inds]
  rlaps_all, zs_all, t0s_all = rlaps_all[inds], zs_all[inds], t0s_all[inds]
     
  t_spec_best, sne_best, rlap_best, zs_best, t0_best, err_t0_best, t0, err_t0 = compute_snid_time(t_specs_all, sne_match_all, rlaps_all, zs_all, t0s_all, N_bests*N_good_spectra)
  
  t_specs_all, sne_match_all, rlaps_all, zs_all, t0s_all, N_good_spectra = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), 0
  for JD in JD_list:
      if (JD - JD_fd)/(1.0 + z) - t0 > 40.0:
          print('Discarding spectrum at JD='+str(JD)+': it is at >40 days since explosion.')
      else:
          N_good_spectra = N_good_spectra+1
          t_specs_all, sne_match_all = np.append(t_specs_all, t_specs_list[str(JD)]), np.append(sne_match_all, sne_match_list[str(JD)])
          rlaps_all, zs_all, t0s_all = np.append(rlaps_all, rlaps_list[str(JD)]), np.append(zs_all, zs_list[str(JD)]), np.append(t0s_all, t0s_list[str(JD)])
  
  #sort by rlap
  inds       = rlaps_all.argsort()[::-1]
  t_specs_all, sne_match_all = t_specs_all[inds], sne_match_all[inds]
  rlaps_all, zs_all, t0s_all = rlaps_all[inds], zs_all[inds], t0s_all[inds]
  
  t_spec_best, sne_best, rlap_best, zs_best, t0_best, err_t0_best, t0, err_t0 = compute_snid_time(t_specs_all, sne_match_all, rlaps_all, zs_all, t0s_all, N_bests*N_good_spectra)
          
  if np.isnan(t0) == True:
      return np.array([]), np.array([]), np.array([]), np.array([])
  else:
      return t_spec_best, sne_best, t0_best, err_t0_best
