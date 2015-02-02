import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

def fit_line(x_data,y_data,errors):
	'''Fits a line using weighted least squares
	   returns a data fit statsmodels object and the residuals'''
	X = np.column_stack((x_data,np.ones(len(x_data))))
	
	errors = np.array([0.1 for i in xrange(len(x_data))])
	
	data_fit = sm.WLS(y_data,X,weights=1/errors**2).fit()
	bse = data_fit.bse
	
	##(1/n)*(|O - E|/O)
	resids2 = abs(np.exp(data_fit.fittedvalues) - np.exp(y_data))
	jstat = np.sum(resids2/np.exp(y_data))/len(y_data)
	
	##ssr is the sum of the residuals over the errors squared ie chi_squared
	##divide by N - 2 as fitting two paramaters to get chi_reduced
	chi_red = data_fit.ssr/(len(y_data)-2)
	if str(chi_red)=='inf': chi_red = 0
	
	return data_fit,jstat,bse,chi_red

def extrap(freq,SI,intercept):
	return np.exp((np.log(freq)*SI)+intercept)

def do_weights(sources):
	flux_s = [source[6] for source in sources]
	weights = [flux/sum(flux_s) for flux in flux_s]
	return weights

def weight_flux(weights,f_to_wait):
	f_weighted = np.array(weights)*f_to_wait
	print f_weighted
	
def comb_pos(sources):
	weights = do_weights(sources)
	ras_to_comb = [source[1] for source in sources]
	decs_to_comb = [source[3] for source in sources]
	ra_w = np.dot(ras_to_comb,weights)
	dec_w = np.dot(decs_to_comb,weights)
	print ra_w,dec_w
	ax2.plot(ra_w,dec_w,marker='+',color='k',linestyle='none',label='Combined',linewidth=2.0,markeredgewidth=3,markersize=14)

fig = plt.figure(figsize=(15,8))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

v3 = ['J223133-082422', 337.89135, 0.00239, -8.40619, 0.00215, 74.00000, 1.37858, 0.19718]
n4 = ['223133-082431', 337.88808, 0.00017, -8.40886, 0.00017, 1400.00000, 0.23530, 0.00770]
m1 = ['2229-086', 337.94410, 0.00125, -8.40920, 0.00139, 408.00000, 0.89000, 0.04000]
v1 = ['J223144-082445', 337.93437, 0.00171, -8.41250, 0.00162, 74.00000, 3.64113, 0.30827]
v2 = ['J223150-082444', 337.96222, 0.00171, -8.41246, 0.00162, 74.00000, 0.81108, 0.11785]
n1 = ['223143-082433', 337.93233, 0.00182, -8.40931, 0.00182, 1400.00000, 0.24400, 0.00760]
n2 = ['223147-082539', 337.94637, 0.00182, -8.42769, 0.00182, 1400.00000, 0.00980, 0.00090]
n3 = ['223150-082442', 337.96187, 0.00017, -8.41167, 0.00017, 1400.00000, 0.26650, 0.00870]

g1 = [3.64113,0.24400] ##v1,n1
f1 = [74.0,1400.0]

g2 = [0.81108,0.26650] ##v2,n3
f2 = [74.0,1400.0]

g3 = [1.37858,0.23530] ##v3,n4
f3 = [74.0,1400.0]

vlssr = [v1,v2,v3]
nvss = [n1,n2,n3,n4]
sumss = []
A154 = []
A182 = []
mrc = [m1]
mwacs = []

#nvss_source1 = ['meh','','','','','',0.0899,'']
#nvss_source2 = ['meh','','','','','',0.0954,'']
#sumss_source1 = ['meh','','','','','',0.1564,'']

#weights1 = do_weights([s1,s2])
#weight_flux(weights1,0.51246)

#weights1 = do_weights([v1,v2])
#weights2 = do_weights([n2,n3])
###weights3 = do_weights([n3,n2])
#weights = 0.5*(np.array(weights1)+np.array(weights2))
###weights = (1/3.0)*(np.array(weights1)+np.array(weights2)+np.array(weights3))
#weight_flux(weights,1.11000)

#comb_pos([n1,n2])
#comb_pos([n3,n4])

comb_flux = [g1,g2,g3]
comb_freqs = [f1,f2,f3]
for flux,freq in zip(comb_flux,comb_freqs):
	ax1.plot(np.log(freq),np.log(flux),label='g%d' %(comb_flux.index(flux)+1) )
	data_fit,jstat,bse,chi_red = fit_line(np.log(freq),np.log(flux),'toot')
	print data_fit.params[0]
	#ext_freqs = [120.0,180.0]
	#ext_fluxes = extrap(np.array(ext_freqs),data_fit.params[0],data_fit.params[1])
	#ax1.plot(np.log(freq),np.log(extrap(np.array(freq),data_fit.params[0],data_fit.params[1])),'c--',label='fit')
	#print ext_fluxes
	#print 'FREQ 120.0e+6 %.5f 0 0 0' %ext_fluxes[0]
	#print 'FREQ 180.0e+6 %.5f 0 0 0' %ext_fluxes[1]

for i in xrange(len(comb_flux)):
	print '-----------'
	for j in xrange(len(comb_flux[i])):
		print 'FREQ %.1fe+6 %.5f 0 0 0' %(comb_freqs[i][j],comb_flux[i][j])
		
####==========================================================================

	
#flux_s = [0.2525,0.05520]
#weights = [flux/sum(flux_s) for flux in flux_s]
#print weights

##weights = 0.5*(np.array([0.6901036587568173, 0.30989634124318277]) + np.array([0.5034777956126271, 0.4965222043873729]))
#f_to_wait = 1.49414
#f_weighted = np.array(weights)*f_to_wait
#print f_weighted

#flux_s = [0.00880,0.09050,0.07420,0.07900]
#weights = [flux/sum(flux_s) for flux in flux_s]
###print weights

#ras_to_comb = [85.67975,85.69112,85.70446,85.72279]
#decs_to_comb = [-21.46803,-21.44894,-21.45406,-21.45631]
#ra_w = np.dot(ras_to_comb,weights)
#dec_w = np.dot(decs_to_comb,weights)
#print ra_w,dec_w
#ax2.plot(ra_w,dec_w,marker='+',color='k',linestyle='none',label='Combined',linewidth=2.0,markeredgewidth=3,markersize=14)
#ax2.plot(78.4117713944,-30.4298500816,marker='+',color='k',linestyle='none',label='Combined',linewidth=2.0,markeredgewidth=3,markersize=14)

shapes = ['o','^','*','s','D','>','<','p','h','1','2','3','8','+','x']
#shapes = ['$1$','$2$','$3$','$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$12$','$13$','$14$','$15$',]
sizes = [6,8,8,6,6,8,8,8,8,8,8,8,8,8,8]

for source in mwacs:
	s_ind = mwacs.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/flux,marker=shapes[s_ind],color='m',linestyle='none',label='mw%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='m',linestyle='none',ms=sizes[s_ind])
	
for source in vlssr:
	s_ind = vlssr.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/flux,marker=shapes[s_ind],color='b',linestyle='none',label='v%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='b',linestyle='none',ms=sizes[s_ind])
	
for source in nvss:
	s_ind = nvss.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/flux,marker=shapes[s_ind],color='g',linestyle='none',label='n%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='g',linestyle='none',ms=sizes[s_ind])
	
for source in sumss:
	s_ind = sumss.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/flux,marker=shapes[s_ind],color='r',linestyle='none',label='s%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='r',linestyle='none',ms=sizes[s_ind])
	
for source in mrc:
	s_ind = mrc.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/flux,marker=shapes[s_ind],color='y',linestyle='none',label='m%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='y',linestyle='none',ms=sizes[s_ind])
	
for source in A154:
	s_ind = A154.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/flux,marker=shapes[s_ind],color='c',linestyle='none',label='A15_%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='c',linestyle='none',ms=sizes[s_ind])
	
for source in A182:
	s_ind = A182.index(source)
	name,ra,rerr,dec,derr,freq,flux,ferr = source
	ax1.errorbar(np.log(freq),np.log(flux),ferr/np.log(flux),marker=shapes[s_ind],color='#CCFF00',linestyle='none',label='A18_%d' %(s_ind+1),ms=sizes[s_ind])
	ax2.errorbar(ra,dec,rerr,derr,marker=shapes[s_ind],color='#CCFF00',linestyle='none',ms=sizes[s_ind])
	
ax1.set_ylabel('log flux')
ax1.set_xlabel('log freq')

ax2.set_xlabel('RA')
ax2.set_ylabel('Dec')
	
ax2.invert_xaxis()

fig.subplots_adjust(top=0.8)
ax1.legend(bbox_to_anchor=(1.0,1.2),loc='upper center',prop={'size':11},ncol=7) 

plt.show()
