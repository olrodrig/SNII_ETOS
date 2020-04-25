import numpy as np
from sklearn import mixture

def snii_templates_epochs():
    
  JD_ln_template, JD_fd_template = {}, {}
  JD_ln_template['1986L']  , JD_fd_template['1986L']   = 2446705.5  ,  2446711.1
  JD_ln_template['1990E']  , JD_fd_template['1990E']   = 2447932.5  ,  2447937.62
  JD_ln_template['1999br'] , JD_fd_template['1999br']  = 2451272.9  ,  2451280.9
  JD_ln_template['1999em'] , JD_fd_template['1999em']  = 2451471.95 ,  2451479.51
  JD_ln_template['1999gi'] , JD_fd_template['1999gi']  = 2451515.68 ,  2451522.32
  JD_ln_template['1999go'] , JD_fd_template['1999go']  = 2451527.7  ,  2451535.7
  JD_ln_template['2000dc'] , JD_fd_template['2000dc']  = 2451758.8  ,  2451765.8
  JD_ln_template['2000dj'] , JD_fd_template['2000dj']  = 2451785.487,  2451795.9
  JD_ln_template['2000el'] , JD_fd_template['2000el']  = 2451835.7  ,  2451840.6
  JD_ln_template['2001X']  , JD_fd_template['2001X']   = 2451958.0  ,  2451968.3
  JD_ln_template['2001do'] , JD_fd_template['2001do']  = 2452131.7  ,  2452135.7
  JD_ln_template['2001fa'] , JD_fd_template['2001fa']  = 2452195.9  ,  2452200.9
  JD_ln_template['2002an'] , JD_fd_template['2002an']  = 2452292.04 ,  2452297.02
  JD_ln_template['2002ce'] , JD_fd_template['2002ce']  = 2452369.7  ,  2452375.378
  JD_ln_template['2002gd'] , JD_fd_template['2002gd']  = 2452549.28 ,  2452550.53
  JD_ln_template['2003Z']  , JD_fd_template['2003Z']   = 2452660.2  ,  2452669.2
  JD_ln_template['2003bn'] , JD_fd_template['2003bn']  = 2452691.5  ,  2452692.83
  JD_ln_template['2003ej'] , JD_fd_template['2003ej']  = 2452770.8  ,  2452779.8
  JD_ln_template['2003hg'] , JD_fd_template['2003hg']  = 2452860.9  ,  2452869.9
  JD_ln_template['2003hl'] , JD_fd_template['2003hl']  = 2452863.0  ,  2452872.0
  JD_ln_template['2003iq'] , JD_fd_template['2003iq']  = 2452918.47 ,  2452921.458
  JD_ln_template['2004ci'] , JD_fd_template['2004ci']  = 2453168.9  ,  2453171.8
  JD_ln_template['2004er'] , JD_fd_template['2004er']  = 2453269.88 ,  2453273.9
  JD_ln_template['2004et'] , JD_fd_template['2004et']  = 2453270.517,  2453271.483
  JD_ln_template['2004fc'] , JD_fd_template['2004fc']  = 2453292.89 ,  2453295.124
  JD_ln_template['2004fx'] , JD_fd_template['2004fx']  = 2453300.92 ,  2453306.93
  JD_ln_template['2005ay'] , JD_fd_template['2005ay']  = 2453449.121,  2453456.58
  JD_ln_template['2005cs'] , JD_fd_template['2005cs']  = 2453548.43 ,  2453549.41
  JD_ln_template['2005dz'] , JD_fd_template['2005dz']  = 2453615.8  ,  2453623.71
  JD_ln_template['2006Y']  , JD_fd_template['2006Y']   = 2453763.09 ,  2453770.08
  JD_ln_template['2006bc'] , JD_fd_template['2006bc']  = 2453811.087,  2453819.15
  JD_ln_template['2006bp'] , JD_fd_template['2006bp']  = 2453833.677,  2453834.647
  JD_ln_template['2006it'] , JD_fd_template['2006it']  = 2454004.69 ,  2454009.67
  JD_ln_template['2006iw'] , JD_fd_template['2006iw']  = 2454009.737,  2454011.798
  JD_ln_template['2007hv'] , JD_fd_template['2007hv']  = 2454342.5  ,  2454352.87
  JD_ln_template['2007il'] , JD_fd_template['2007il']  = 2454345.94 ,  2454353.95
  JD_ln_template['2007pk'] , JD_fd_template['2007pk']  = 2454409.83 ,  2454414.81
  JD_ln_template['2008bh'] , JD_fd_template['2008bh']  = 2454538.57 ,  2454548.66
  JD_ln_template['2008br'] , JD_fd_template['2008br']  = 2454559.323,  2454564.265
  JD_ln_template['2008ho'] , JD_fd_template['2008ho']  = 2454787.77 ,  2454796.61
  JD_ln_template['2008if'] , JD_fd_template['2008if']  = 2454802.73 ,  2454812.71
  JD_ln_template['2008il'] , JD_fd_template['2008il']  = 2454822.69 ,  2454827.64
  JD_ln_template['2008in'] , JD_fd_template['2008in']  = 2454824.45 ,  2454824.95
  JD_ln_template['2009ao'] , JD_fd_template['2009ao']  = 2454886.62 ,  2454894.62
  JD_ln_template['2009bz'] , JD_fd_template['2009bz']  = 2454912.03 ,  2454919.98
  JD_ln_template['2010id'] , JD_fd_template['2010id']  = 2455450.82 ,  2455454.743
  JD_ln_template['2012aw'] , JD_fd_template['2012aw']  = 2456001.769,  2456003.349
  JD_ln_template['2013am'] , JD_fd_template['2013am']  = 2456371.698,  2456373.138
  JD_ln_template['2013by'] , JD_fd_template['2013by']  = 2456402.872,  2456403.752
  JD_ln_template['2013ej'] , JD_fd_template['2013ej']  = 2456497.04 ,  2456497.625
  JD_ln_template['2013fs'] , JD_fd_template['2013fs']  = 2456570.82 ,  2456571.737
  JD_ln_template['2013hj'] , JD_fd_template['2013hj']  = 2456635.7  ,  2456638.8
  JD_ln_template['2014G']  , JD_fd_template['2014G']   = 2456668.35 ,  2456671.111
  JD_ln_template['LSQ14gv'], JD_fd_template['LSQ14gv'] = 2456670.7  ,  2456674.8
  JD_ln_template['2014cx'] , JD_fd_template['2014cx']  = 2456901.89 ,  2456902.90
  JD_ln_template['2014cy'] , JD_fd_template['2014cy']  = 2456898.8  ,  2456900.5
  JD_ln_template['2015bs'] , JD_fd_template['2015bs']  = 2456915.5  ,  2456925.5
  JD_ln_template['2016esw'], JD_fd_template['2016esw'] = 2457607.802,  2457608.814

  return JD_ln_template, JD_fd_template

#compute weighted average
def weighted_average(x, sigma_x, with_intrinsic_error=True):
    
  if len(x) > 1:
      if with_intrinsic_error:
          residuals = x - np.mean(x)
          rms       = np.sqrt(np.sum(residuals**2)/float(len(residuals)-1))
          sigma_0s  = np.linspace(0.0, rms, 100)
      else:
          sigma_0s  = np.array([0.0])
          
      m2lnL_min = 1.e90
      for sigma_0 in sigma_0s:
          
          Var   = sigma_x**2 + sigma_0**2
          w_ave = np.sum(x/Var)/np.sum(1.0/Var)
          m2lnL = np.sum(np.log(Var)+(x-w_ave)**2/Var)
          
          if m2lnL < m2lnL_min:
              m2lnL_min = m2lnL
              best_x    = w_ave
              best_error = np.sqrt(1.0/np.sum(1.0/Var))
  else:
      best_x, best_error = x[0], sigma_x[0]
      
  return best_x, best_error

#pick rangom values given a pdf
def values_from_distribution(x, pdf, N):

  x_sample = np.random.choice(x, N, p=pdf/np.sum(pdf)) #sum of probabilities must to be 1

  return x_sample

#Simpson's rule
def simpson(x,f):
  
  integral = (f[0] + f[-1]) / 3.0  #extremes
  n = len(x)
  
  four = "o"
  for i in range(1, n - 1):
      if four == "o":
          integral += f[i] * 4.0 / 3.0
          four = "x"
      else:
          integral += f[i] * 2.0 / 3.0
          four = "o"
              
  integral = (x[1] - x[0]) * integral
      
  return integral

#discard possible outliers through the Tukey's rule
def tukey_rule(x, k=1.5):

  Q1, Q3 = np.quantile(x, [0.25, 0.75])
  IQR    = Q3 - Q1
  
  x = x[x>=Q1-k*IQR]
  x = x[x<=Q3+k*IQR]
  
  return x

#return a normalized gaussian pdf
def gaussian_pdf(mu, sigma, x_sampled):
  g_sampled = np.exp(-0.5*(mu-x_sampled)**2/sigma**2)
  g_sampled = g_sampled / simpson(x_sampled, g_sampled)
  return g_sampled

#return a uniform pdf
def uniform_pdf(x_min, x_max, x):
  h = 1.0/(x_max-x_min)
  pdf = np.linspace(h, h, len(x))
  for i in range(0, len(x)):
      if x[i] < x_min or x[i] > x_max:
          pdf[i] = 0.0
  return pdf

#return a pdf computed as a mixture of Gaussians
def get_pdf(y, y_sampled, max_components=2):
  
  x, x_sampled = y.reshape(-1,1), y_sampled.reshape(-1,1)
  
  BIC_min = 1.e90
  for n_components in range(1, max_components+1):
      gmm   = mixture.GaussianMixture(n_components=n_components)
      model = gmm.fit(x)
      BIC   = model.bic(x)
      if BIC < BIC_min:
          BIC_min   = BIC
          model_min = model         
  ln_pdf = model_min.score_samples(x_sampled)
  pdf    = np.exp(ln_pdf)
  
  return pdf

#return different pdf's
def final_pdfs(z, JD_ln, JD_fd, pdfs_per_sn, x_sampled, N_sample, rms_t0):
    
  #define the uniform pdf's given by the JD_fd and JD_fd+JD_ln
  ln, fd    = (JD_ln - JD_fd)/(1.0+z), 0.0
  pdf_fd    = uniform_pdf(-9999.0, fd, x_sampled) #fd as prior
  pdf_fd_ln = uniform_pdf(ln, fd, x_sampled) #fd and ln as prior
    
  #combine the pdf of different sne
  pdf_snid = np.linspace(1.0, 1.0, len(x_sampled))  
  for pdf_per_sn in pdfs_per_sn:
      pdf_snid = pdf_snid*pdf_per_sn
  
  #add typical rms(t0) error
  err_0    = np.random.normal(0.0, rms_t0, N_sample)
  err_0    = np.random.choice(tukey_rule(err_0), N_sample)
  err_0    = err_0 - np.median(err_0)
  t0s_snid = values_from_distribution(x_sampled, pdf_snid, N_sample)
  t0s_snid = t0s_snid + err_0
  t0s_snid = np.random.choice(tukey_rule(t0s_snid),N_sample)
  
  #compute pdf's
  pdf_snid       = get_pdf(t0s_snid, x_sampled, max_components=1)
  pdf_snid_fd    = pdf_snid*pdf_fd
  pdf_snid_fd_ln = pdf_snid*pdf_fd_ln
  
  #normalize pdf's
  pdf_snid       = pdf_snid / simpson(x_sampled, pdf_snid)
  pdf_snid_fd    = pdf_snid_fd / simpson(x_sampled, pdf_snid_fd)
  pdf_snid_fd_ln = pdf_snid_fd_ln / simpson(x_sampled, pdf_snid_fd_ln)
  
  return pdf_fd_ln, pdf_snid, pdf_snid_fd, pdf_snid_fd_ln


def average_pdf_per_sn_bm_with_t0_error(sne_bm, t0s_bm, rms_t0s_bm, x_pdf, N_sample):
  
  JD_ln_template, JD_fd_template = snii_templates_epochs()
  
  pdfs_per_sn = []
  for sn_bm, spec_phase, err_spec_phase in zip(sne_bm, t0s_bm, rms_t0s_bm):
            
      delta = round(JD_fd_template[sn_bm]-JD_ln_template[sn_bm],3)
      rms_uniform = delta/np.sqrt(12.0)
      
      if rms_uniform < 0.3*err_spec_phase:
          pdf_per_sn = gaussian_pdf(spec_phase, err_spec_phase, x_pdf)
      else:
          err_t0_template = np.random.uniform(-0.5*delta,0.5*delta, N_sample)
          err_t0_template = err_t0_template - np.median(err_t0_template) #center the distibution to zero
          
          x1 = np.random.normal(spec_phase, err_spec_phase, N_sample)
          x1 = np.random.choice(tukey_rule(x1), N_sample)
          x1 = x1 - np.median(x1) + spec_phase #center the distibution to the phase
          #include values from the uniform distrution
          x= x1 + err_t0_template
          pdf_per_sn = get_pdf(x, x_pdf)
          pdf_per_sn = pdf_per_sn / simpson(x_pdf, pdf_per_sn)
      pdfs_per_sn.append(pdf_per_sn)

  return pdfs_per_sn

def average_pdf_per_sn_bm(t0s_best, rms_t0s_best, sne_best):
    
  #best matching SNe
  sne_bm = list(set(sne_best))

  t0s_bm, rms_t0s_bm = [], []
  for sn_bm in sne_bm:
      phases, err_phases = np.array([]), np.array([])
      for sn_i, spec_phase, err_spec_phase in zip(sne_best, t0s_best, rms_t0s_best):
          if sn_i == sn_bm:
              phases     = np.append(phases, spec_phase)
              err_phases = np.append(err_phases, err_spec_phase)
      t0_best, rms_t0_best =  weighted_average(phases, err_phases)
      t0s_bm.append(t0_best)
      rms_t0s_bm.append(rms_t0_best)
    
  return sne_bm, t0s_bm, rms_t0s_bm

def typical_pdf_per_sn_bm_per_spectrum(t0s_best, rms_t0s_best, sne_best, t_spec_best):

  epochs = list(set(t_spec_best))
  
  new_sne_best, new_t0_best, new_rms_t0_best = [], [], []
  for epoch in epochs:
      
      #number of templates at the epoch
      templates = []
      for t, sn_best in zip(t_spec_best, sne_best):
          if t == epoch:  templates.append(sn_best)
      templates = list(set(templates))
      
      for template in templates:
          
          phases, err_phases = np.array([]), np.array([])
          for t, sn_best, spec_phase, err_spec_phase in zip(t_spec_best, sne_best, t0s_best, rms_t0s_best):
              if t == epoch and sn_best == template:
                  phases     = np.append(phases, spec_phase)
                  err_phases = np.append(err_phases, err_spec_phase)

          t0_best, rms_t0_best =  weighted_average(phases, err_phases)
          new_sne_best.append(template)
          new_t0_best.append(t0_best)
          new_rms_t0_best.append(rms_t0_best*np.sqrt(float(len(phases))))
  
  sne_best, t0_best, rms_t0_best = np.array(new_sne_best), np.array(new_t0_best), np.array(new_rms_t0_best)
  
  return sne_best, t0_best, rms_t0_best