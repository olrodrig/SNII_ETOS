import numpy as np
from snid_routines import snid_results
from pdf_routines import typical_pdf_per_sn_bm_per_spectrum, average_pdf_per_sn_bm, average_pdf_per_sn_bm_with_t0_error, final_pdfs, values_from_distribution
from plot_ETOS import plot_ETOS

class Explosion_Time_from_Optical_Spectra:
  def __init__(self, templates_path, sn, z, JD_spectra, wl_spectra, f_spectra, JD_fd, JD_ln=-9999, snid_verbose=False, display_snid_plots=False):
    
    #transform to the rest frame
    JD, wl, f = [], [], []
    for JD_spectrum, wl_spectrum, f_spectrum in zip(JD_spectra, wl_spectra, f_spectra):
        
        wl_spectrum = wl_spectrum/(1.0+z)
        f_spectrum  = f_spectrum*(1.0+z)
        
        #discard spectra at >40 days since the first detection
        if (JD_spectrum - JD_fd) / (1.0+z) >= 40.0:
            print("discarding spectrum at JD="+str(JD_spectrum)+": it is at >40 days since the first detection")
        else:
            JD.append(JD_spectrum)
            wl.append(wl_spectrum)
            f.append(f_spectrum)
            if min(wl_spectrum) > 4100.0:  print("WARNING: spectrum at JD="+str(JD_spectrum)+" has a minimum wavelength higher than the requested (<=410 nm)")
            if max(wl_spectrum) < 7000.0:  print("WARNING: spectrum at JD="+str(JD_spectrum)+" has a maximum wavelength lower than the requested (>=700 nm)")
    JD_spectra, wl_spectra, f_spectra = JD, wl, f
    

    sne_bm, t0s_bm, rms_t0s_bm = {}, {}, {}
    
    #1) run snid for each of the input spectra
    t_spec_bm, sne_bm['1'], t0s_bm['1'], rms_t0s_bm['1'] = snid_results(templates_path, sn, z, JD_fd, JD_spectra, wl_spectra, f_spectra, snid_verbose, display_snid_plots)

    if len(t0s_bm['1']) == 0:
        print("WARNING: Explosion_Time_from_Optical_Spectra without results!")
    else:
        n_spec = len(list(set(t_spec_bm)))
        if n_spec == 1: rms_t0 = 5.0
        if n_spec >  1: rms_t0 = 4.1
        xmin, xmax, dx = -55.0, 55.0, 0.01
        x_pdf = np.linspace(xmin, xmax, int(round((xmax-xmin)/dx,0))+1)
        N_sample = 200000
    
        np.random.seed(1)
        
        #2) select the typical pdf per SN per input spectrum
        sne_bm['2'], t0s_bm['2'], rms_t0s_bm['2'] = typical_pdf_per_sn_bm_per_spectrum(t0s_bm['1'], rms_t0s_bm['1'], sne_bm['1'], t_spec_bm)

        #3) compute the average pdf per SN
        sne_bm['3'], t0s_bm['3'], rms_t0s_bm['3'] = average_pdf_per_sn_bm(t0s_bm['2'], rms_t0s_bm['2'], sne_bm['2'])

        #4) include uncertainty due to the t0
        pdfs_snid_sn_bm = average_pdf_per_sn_bm_with_t0_error(sne_bm['3'], t0s_bm['3'], rms_t0s_bm['3'], x_pdf, N_sample)

        #5) compute the final pdf's
        pdf = {}
        pdf['fd+ln'], pdf["snid"], pdf["snid+fd"], pdf["snid+fd+ln"] = final_pdfs(z, JD_ln, JD_fd, pdfs_snid_sn_bm, x_pdf, N_sample, rms_t0)  

        
        JD_t0, rms_JD_t0 = {"fd+ln":{}}, {"fd+ln":{}}
        t0, rms_t0 = {"snid":{}, "snid+fd":{}, "snid+fd+ln":{}}, {"snid":{}, "snid+fd":{}, "snid+fd+ln":{}}
        for pdf_name in ["fd+ln", "snid", "snid+fd", "snid+fd+ln"]:
            
            if pdf_name == "fd+ln":
            
                for method in ['mean', 'median', 'mode']:
                    if JD_ln != -9999:
                        JD_t0[pdf_name][method] = round((JD_ln + JD_fd)*0.5, 2)
                        rms_JD_t0[pdf_name][method] = round((JD_fd - JD_ln)/np.sqrt(12.0), 2)
                    else:
                        JD_t0[pdf_name][method] = -9999
                        rms_JD_t0[pdf_name][method] = -9999
            else:
            
                t0s = values_from_distribution(x_pdf, pdf[pdf_name], N_sample*10)
                
                JD_t0[pdf_name], rms_JD_t0[pdf_name] = {}, {}
                for method in ['mean', 'median', 'mode']:
                    
                    if method == 'mean':
                        t0[pdf_name][method] = np.mean(t0s)
                    if method == 'median':
                        pc_func = np.array([0.0])
                        for i in range(1, len(x_pdf)):
                            integral = (x_pdf[i]-x_pdf[i-1])*(pdf[pdf_name][i]+pdf[pdf_name][i-1])*0.5
                            pc_func = np.append(pc_func, pc_func[i-1]+integral)
                        pc_func = pc_func / max(pc_func)
                        t0[pdf_name][method] = np.interp(0.5, pc_func, x_pdf)
                    if method == 'mode':
                        i_mode          = pdf[pdf_name].tolist().index(max(pdf[pdf_name]))
                        xs_for_mode     = x_pdf[i_mode-2:i_mode+3]
                        t0s_for_mode    = pdf[pdf_name][i_mode-2:i_mode+3]
                        pars            = np.polyfit(xs_for_mode, t0s_for_mode, 2)[::-1]
                        t0[pdf_name][method] = -pars[1]/(2.0*pars[2])
                        
                        if t0["snid"][method] < 0.0:
                            t0["snid+fd"][method] = t0["snid"][method]
                            if t0["snid"][method] > (JD_ln - JD_fd)/(1.0+z):
                                t0["snid+fd+ln"][method] = t0["snid"][method]
                            else:
                                t0["snid+fd+ln"][method] = (JD_ln - JD_fd)/(1.0+z)
                        else:
                            t0["snid+fd"][method] = 0.0
                            t0["snid+fd+ln"][method] = 0.0
                
                    #rms
                    residuals   = t0s - t0[pdf_name][method]
                    rms_t0[pdf_name][method] = np.sqrt(np.sum(residuals**2)/float(len(residuals)-1))
                    #JD_t0 and rms
                    JD_t0[pdf_name][method]     = round(t0[pdf_name][method]*(1.0+z) + JD_fd, 2)
                    rms_JD_t0[pdf_name][method] = round(rms_t0[pdf_name][method]*(1.0+z), 2)
            
        #input values
        self.sn     = sn
        self.z      = z
        self.JD_fd  = JD_fd
        self.JD_ln  = JD_ln
        self.n_spec = n_spec
        #final pdf's
        self.x_pdf           = x_pdf
        self.pdf             = pdf
        #explosion times and rms
        self.JD_t0 = JD_t0
        self.rms_JD_t0 = rms_JD_t0
        #values to construct figures
        self.for_plots = [sne_bm, t0s_bm, rms_t0s_bm, pdfs_snid_sn_bm]
  
  #plot
  def plot(self, figure_name=''):
      
      plot_ETOS(self, figure_name)
  
  #pick N random values from a specific pdf
  def random_values(self, pdf_name, N):
           
      
      if pdf_name not in ["fd+ln", "snid", "snid+fd", "snid+fd+ln"]:
          print('ERROR in random_values: '+pdf_name+' is not valid! Use "fd+ln", "snid", "snid+fd", of "snid+fd+ln".')
      else:
          JD_pdf = self.x_pdf*(1.0+self.z) + self.JD_fd
          JD_t0s = values_from_distribution(JD_pdf, self.pdf[pdf_name], N)
          
          return JD_t0s
  
  #compute the boundaries enclosing a certain per cent of a pdf
  def pc_boundaries(self, pdf_name, method, percent):
      
      if pdf_name not in ["fd+ln", "snid", "snid+fd", "snid+fd+ln"]:
          print('ERROR in pc_boundaries: '+pdf_name+' is not valid! Use "snid", "snid+fd", of "snid+fd+ln".')
      elif method not in ["mean", "median", "mode"]:
          print('ERROR in pc_boundaries: '+method+' is not valid! Use "mean", "median", of "mode".')
      elif percent < 0.0 or percent > 100.0:
          print('ERROR in pc_boundaries: '+str(percent)+' is not valid! Use a value between 0 and 100.')
      else:
          if pdf_name == "fd+ln":
              if self.JD_ln != -9999:
                  JD_l = 0.5*(self.JD_fd*(1.0-percent/100.0)+self.JD_ln*(1.0+percent/100.0))
                  JD_u = 0.5*(self.JD_fd*(1.0+percent/100.0)+self.JD_ln*(1.0-percent/100.0))
              else:
                  JD_l, JD_u = -9999, -9999
          else:
              JD_t0 = self.JD_t0[pdf_name][method]
              pdf   = self.pdf[pdf_name]
              JD_pdf = self.x_pdf*(1.0+self.z) + self.JD_fd
                       
              pc_func = np.array([0.0])
              for i in range(1, len(JD_pdf)):
                  integral = (JD_pdf[i]-JD_pdf[i-1])*(pdf[i]+pdf[i-1])*0.5
                  pc_func = np.append(pc_func, pc_func[i-1]+integral)
              pc_func = pc_func / max(pc_func)
              
              if JD_t0 == self.JD_ln:
                  JD_l = JD_t0
                  JD_u = np.interp(percent/100.0, pc_func, JD_pdf)
              elif JD_t0 == self.JD_fd:
                  JD_u = JD_t0
                  JD_l = np.interp(1.0-percent/100.0, pc_func, JD_pdf)
              else:
                  alpha_c = round(np.interp(JD_t0, JD_pdf, pc_func),3)
                  alpha_l = alpha_c - 0.5*percent/100.0
                  alpha_u = alpha_c + 0.5*percent/100.0
                  if alpha_l <= 0.0:
                      if "ln" in pdf_name:
                          JD_l = self.JD_ln
                      else:
                          JD_l = min(JD_pdf)
                      JD_u = np.interp(percent/100.0, pc_func, JD_pdf)
                  elif alpha_u >= 1.0:
                      JD_u = self.JD_fd
                      JD_l = np.interp(1.0-percent/100.0, pc_func, JD_pdf)
                  else: 
                      JD_l = np.interp(alpha_l, pc_func, JD_pdf)
                      JD_u = np.interp(alpha_u, pc_func, JD_pdf)
              
              if percent == 100.0 and pdf_name in ["snid+fd", "snid+fd+ln"]:
                  JD_u = self.JD_fd
                  if pdf_name == "snid+fd+ln":
                      JD_l = self.JD_ln
              
          boundaries = np.array([round(JD_l,2), round(JD_u,2)])
          
          return boundaries
