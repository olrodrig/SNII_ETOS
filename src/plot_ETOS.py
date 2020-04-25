import os
import matplotlib
import matplotlib.pyplot as plt
from pdf_routines import gaussian_pdf
import matplotlib.ticker as plticker

matplotlib.rc('axes', labelsize=15.0)
textsize = 11

def draw_ticks(ax, i, x_major, x_minor):
    
  ax.yaxis.set_major_formatter(plticker.NullFormatter())
  ax.yaxis.set_ticks_position('none') 
  
  ax.xaxis.set_major_locator(plticker.MultipleLocator(base=x_major))
  ax.xaxis.set_minor_locator(plticker.MultipleLocator(base=x_minor))
  
  ax.tick_params('x', length=6, width=1, which='major', direction='in')
  ax.tick_params('x', length=3, width=1, which='minor', direction='in')
  ax.xaxis.set_ticks_position('both')
  
  ax.set_ylabel('pdf($t_0$)', labelpad=-1)
  
def color_line_style(templates):
  colors = ['b','lime','r','c','m','y','gold','maroon','darkgreen','darkgray','violet','darkkhaki','teal']
  lines  = ['-','--',':','-.']
  template_color, template_line, i_color, i_line = {}, {}, 0, 0
  for i in range(0, len(templates)):
      template_color[templates[i]] = colors[i_color]
      template_line[templates[i]]  = lines[i_line]
      i_color = i_color+1
      if i_color == len(colors):
          i_color = 0
          i_line  = i_line+1
  return template_color, template_line

def plot_gaussian_pdfs(sne_bm, t0s_bm, rms_t0s_bm, x_pdf, ax, xpos, ypc, text, template_line, template_color, textsize, panel, with_labels=False, label_pos=0.0):
  ymax = []
  for sn_bm, spec_phase, err_spec_phase in zip(sne_bm, t0s_bm, rms_t0s_bm):
      g_per_spec = gaussian_pdf(spec_phase, err_spec_phase, x_pdf)
      ax.plot(x_pdf, g_per_spec, linestyle=template_line[sn_bm], color=template_color[sn_bm], label=sn_bm)
      ymax.append(max(g_per_spec))
  ymax = max(ymax)
  ax.set_ylim(0, ymax*1.1)
  ax.text(xpos, ymax*ypc, text, fontsize=textsize)
  ax.text(-54, ymax*ypc, panel, fontsize=textsize)
  xmin, xmax = min(x_pdf), max(x_pdf)
  ax.set_xlim(xmin, xmax)
  if with_labels:  ax.legend(bbox_to_anchor=(0.0, label_pos, 1.0, 1.0), loc=3, ncol=8, borderaxespad=0.0, labelspacing=0.4, handletextpad=0.2, mode='expand', numpoints=1)
  
def plot_non_gaussian_pdfs(sne_bm, pdfs_per_sn, x_pdf, ax, xpos, ypc, text, template_line, template_color, textsize, panel):
  ymax = []
  for sn_bm, pdf_per_sn in zip(sne_bm, pdfs_per_sn):
      ax.plot(x_pdf, pdf_per_sn, linestyle=template_line[sn_bm], color=template_color[sn_bm], label=sn_bm)
      ymax.append(max(pdf_per_sn))
  ymax = max(ymax)
  ax.set_ylim(0, ymax*1.1)
  ax.text(xpos, ymax*ypc, text, fontsize=textsize)
  ax.text(-54, ymax*ypc, panel, fontsize=textsize)
  xmin, xmax = min(x_pdf), max(x_pdf)
  ax.set_xlim(xmin, xmax)

def plot_final_pdfs(pdf, x_pdf, JD_ln, n_spec, ax):
  if JD_ln != -9999:  ax.plot(x_pdf, pdf['fd+ln'], '-', color='dimgray', lw=2)
  ax.plot(x_pdf, pdf['snid'], '-k', lw=3)
  ax.plot(x_pdf, pdf['snid+fd'], '-', color='lime', lw=1)
  if JD_ln != -9999:  ax.plot(x_pdf, pdf['snid+fd+ln'], '--', color='r', lw=1)
  
  ymax = max([max(pdf['fd+ln']),max(pdf['snid']),max(pdf['snid+fd']),max(pdf['snid+fd+ln'])])
  ax.set_ylim(0, ymax*1.1)
  xpos = 14
  if JD_ln != -9999:  ax.text(xpos, ymax*0.82, "Uniform pdf between JD$_\mathrm{fd}$ and JD$_\mathrm{ln}$", color='dimgray', fontsize=textsize)
  if n_spec == 1:
      ax.text(xpos, ymax*0.65, "SNID pdf (+ a rms error of 5.0 days)", color='k', fontsize=textsize)
  else:
      ax.text(xpos, ymax*0.65, "SNID pdf (+ a rms error of 4.1 days)", color='k', fontsize=textsize)
  ax.text(xpos, ymax*0.48, "SNID pdf with JD$_\mathrm{fd}$ as prior", color='lime', fontsize=textsize)
  if JD_ln != -9999:  ax.text(xpos, ymax*0.31, "SNID pdf with JD$_\mathrm{fd}$ and JD$_\mathrm{ln}$ as priors", color='r', fontsize=textsize)
  ax.set_xlim(min(x_pdf), max(x_pdf))
  ax.set_xlabel('$(t-\mathrm{JD}_\mathrm{fd}$)/($1+z$) [day]')
  
  ax.text(-54, ymax*0.82, '(e)', fontsize=textsize)

def plot_ETOS(ETOS, figure_name):
  
  sn     = ETOS.sn
  n_spec = ETOS.n_spec
  sne_bm          = ETOS.for_plots[0]
  t0s_bm          = ETOS.for_plots[1]
  rms_t0s_bm      = ETOS.for_plots[2]
  pdfs_snid_sn_bm = ETOS.for_plots[3]
  x_pdf = ETOS.x_pdf
  pdf   = ETOS.pdf
  
  #define colors and line styles
  template_color, template_line = color_line_style(sne_bm['3'])
  
  fig = plt.figure(figsize=(8, 6))
  if len(sne_bm['3']) <= 8:
      fig.subplots_adjust(wspace=0.0, hspace=0.0, left=0.04, bottom=0.075, right=0.99, top=0.92)
  elif len(sne_bm['3']) <= 16:
      fig.subplots_adjust(wspace=0.0, hspace=0.0, left=0.04, bottom=0.075, right=0.99, top=0.889)
  else:
      fig.subplots_adjust(wspace=0.0, hspace=0.0, left=0.04, bottom=0.075, right=0.99, top=0.858)
      
  ax = {}
  for i in range(1, 6):  ax[str(i)] = fig.add_subplot(5,1,i)

  plot_gaussian_pdfs(sne_bm['1'], t0s_bm['1'], rms_t0s_bm['1'], x_pdf, ax['1'], 19, 0.82, "all pdf's per SN$_\mathrm{bm}$ per spectrum", template_line, template_color, textsize, '(a)')
  plot_gaussian_pdfs(sne_bm['2'], t0s_bm['2'], rms_t0s_bm['2'], x_pdf, ax['2'], 16, 0.82, "typical pdf per SN$_\mathrm{bm}$ per spectrum", template_line, template_color, textsize, '(b)')
  plot_gaussian_pdfs(sne_bm['3'], t0s_bm['3'], rms_t0s_bm['3'], x_pdf, ax['3'], 35, 0.82, "<pdf> per SN$_\mathrm{bm}$",template_line,template_color, textsize, '(c)', with_labels=True,label_pos=3.02)
  plot_non_gaussian_pdfs(sne_bm['3'], pdfs_snid_sn_bm, x_pdf, ax['4'], 20, 0.82, "<pdf> per SN$_\mathrm{bm}$, with $t_0$ error", template_line, template_color, textsize, '(d)')
  plot_final_pdfs(pdf, x_pdf, ETOS.JD_ln, n_spec, ax['5'])
   
  for i in range(1, 6):  draw_ticks(ax[str(i)], i, 10.0, 2.0)
  
  if n_spec == 1:
      title = 'SN '+sn+'; 1 input spectrum'
  else:
      title = 'SN '+sn+'; '+str(n_spec)+' input spectra'
  if len(sne_bm['3']) <= 8:
      ax['1'].set_title(title, position=(0.5,1.25))
  elif len(sne_bm['3']) <= 16:
      ax['1'].set_title(title, position=(0.5,1.45))
  else:
      ax['1'].set_title(title, position=(0.5,1.65))
  
  if figure_name != '':  fig.savefig(figure_name+'.pdf', dpi=300)
