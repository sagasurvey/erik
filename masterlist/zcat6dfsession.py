# coding: utf-8
get_ipython().system(u'ls -F -G ')
get_ipython().magic(u'run masterlist.py')
d=ascii.read('6dF.csv')
d=ascii.read('6dF.csv',guess=False,delimiter=',')
d.dtype
d['z_helio']]
d['z_helio']
d['z_helio']*cnst.c.to('km/s').value
d['z_helio']*cnst.c.to('km/s').value
v=d['z_helio']*cnst.c.to('km/s').value
v<3000
sum(_)
len(v)
d['ra']
d.dtype
import masterlist
reload(masterlist)
res=masterlist.x_match_tests(masterlist)
res=masterlist.x_match_tests(mastercat)
res=masterlist.x_match_tests(mastercat,vcuts=4000*u.km/u.s)
4000*u.km/u.s/cnst.c
_.decompose()
reload(masterlist)
res=masterlist.x_match_tests(mastercat,vcuts=4000*u.km/u.s)
res=masterlist.x_match_tests(mastercat,vcuts=3000*u.km/u.s)
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s)
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=1*u.arcmin)
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=10*u.arcmin)
res=masterlist.x_match_tests(mastercat1,vcuts=4000*u.km/u.s,tol=10*u.arcmin)
res.keys()
res['sixdf_nomatch']
res['sixdf_nomatch']['targetname']
res.keys()
res['sixdf_nomatch']
res['sixdf_catidx']
res.keys()
res=masterlist.x_match_tests(mastercat1,vcuts=4000*u.km/u.s,tol=10*u.arcmin)
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=10*u.arcmin)
res['sixdf_catidx']
len(_)
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=1*u.arcmin)
res['sixdf_matches']
res.keys()
res['sixdf_match']
hist(res['sixdf_match']['obsdec'])
histstep(res['sixdf_match']['obsdec'])
clf()
histstep(res['sixdf_match']['obsdec'])
histstep(res['sixdf_nomatch']['obsdec'])
histstep(res['sixdf_nomatch']['obsdec'],norm=True)
clf()
histstep(res['sixdf_nomatch']['obsdec'],normed=True)
histstep(res['sixdf_match']['obsdec'],normed=True)
ctwomass = ICRS(u.deg*twomassxsc['ra'].view(np.ndarray), u.deg*twomassxsc['dec'].view(np.ndarray))
ctwomass = coo.ICRS(u.deg*twomassxsc['ra'].view(np.ndarray), u.deg*twomassxsc['dec'].view(np.ndarray))
ctwomass.ra.shape
ctwomass.ra.shape
noma=res['sixdf_nomatch']
noma6df=res['sixdf_nomatch']
noma6dfc=coo.ICRS(noma6df['obsra']*u.deg,u.deg*noma6df['obsdec'])
idx,dd,d3d=noma6dfc.match_to_catalog_sky(ctwomass)
dd
clf()
histstep(dd.arcmin)
histstep(dd.arcmin,bins=100,range=(0,5))
clf()
histstep(dd.arcmin,bins=100,range=(0,5))
histstep(dd.arcmin,bins=100,range=(0,5))
clf()
histstep(dd.arcmin,bins=100,range=(0,5))
histstep(dd.arcmin,bins=100,range=(0,1))
clf()
histstep(dd.arcmin,bins=100,range=(0,1))
histstep(dd.arcmin,bins=100,range=(0,1),log=True)
clf()
histstep(dd.arcmin,bins=100,range=(0,1),log=True)
histstep(dd.arcmin,bins=1000,range=(0,1),log=True)
clf()
histstep(dd.arcmin,bins=1000,range=(0,1),log=True)
histstep(dd.arcmin,bins=1000,range=(0,10),log=True)
clf()
histstep(dd.arcmin,bins=1000,range=(0,10),log=True)
histstep(dd.arcsec,bins=1000,range=(0,10),log=True)
clf()
histstep(dd.arcsec,bins=1000,range=(0,10),log=True)
histstep(dd.arcsec,bins=1000,range=(0,3),log=True)
clf()
histstep(dd.arcsec,bins=1000,range=(0,3),log=True)
histstep(dd.arcsec,bins=100,range=(0,3),log=True)
clf()
histstep(dd.arcsec,bins=100,range=(0,3),log=True)
histstep(dd.arcsec,bins=100,range=(0,5),log=True)
clf()
histstep(dd.arcsec,bins=100,range=(0,5),log=True)
msk=dd.arcsec<2
twomass[idx[msk]]
twomass[idx[msk]]
#ctwomass = coo.ICRS(u.deg*twomassxsc['ra'].view(np.ndarray), u.deg*twomassxsc['dec'].view(np.ndarray))
twomassxsc[idx[msk]]
twomassxsc[idx[msk]]['dec']
histstep(twomassxsc[idx[msk]]['dec'])
clf()
histstep(twomassxsc[idx[msk]]['dec'])
histstep(twomassxsc[idx[msk]]['fasd'])
twomassxsc.dtype.names
histstep(twomassxsc[idx[msk]]['k_m_ext'])
clf()
histstep(twomassxsc[idx[msk]]['k_m_ext'])
#idx,dd,d3d=noma6dfc.match_to_catalog_sky(ctwomass)
noma6dfc['z_helio']
noma6df['z_helio']
WMAP9.distmod(noma6df['z_helio'])
#WMAP9.distmod(noma6df['z_helio'])
twomassxsc[idx[msk]]['k_m_ext'] - WMAP9.distmod(noma6df['z_helio'])
twomassxsc[idx[msk]]['k_m_ext'] - WMAP9.distmod(noma6df['z_helio'][msk])
clf()
twomassxsc[idx[msk]]['k_m_ext'] - WMAP9.distmod(noma6df['z_helio'][msk])
hist(_)
clf()
Kabs=twomassxsc[idx[msk]]['k_m_ext'] - WMAP9.distmod(noma6df['z_helio'][msk])
clf()
Kabs
histstep(Kabs[~np.isnan(Kabs)])
clf()
histstep(Kabs[~np.isnan(Kabs)],bins=50)
figure()
histstep(twomassxsc[idx[msk]]['k_m_ext'])
#idx,dd,d3d=noma6dfc.match_to_catalog_sky(ctwomass)
#noma6dfc=coo.ICRS(noma6df['obsra']*u.deg,u.deg*noma6df['obsdec'])
ma6df=res['sixdf_match']
ma6dfc=coo.ICRS(ma6df['obsra']*u.deg,u.deg*ma6df['obsdec'])
idx2,dd2,d3d2=ma6dfc.match_to_catalog_sky(ctwomass)
figure(1)
histstep(dd2.arcsec,bins=100,range=(0,5),log=True)
clf()
histstep(dd2.arcsec,bins=100,range=(0,5),log=True)
clf()
msk2=dd2.arcsec<2
histstep(twomassxsc[idx2[msk2]]['k_m_ext'])
clf()
histstep(twomassxsc[idx2[msk2]]['k_m_ext'])
histstep(twomassxsc[idx2[msk2]]['k_m_ext'],bins=linspace(16,5,50))
clf()
histstep(twomassxsc[idx2[msk2]]['k_m_ext'],bins=linspace(5,16,50))
histstep(twomassxsc[idx2[msk2]]['k_m_ext'],bins=linspace(5,16,50),label='match')
clf()
histstep(twomassxsc[idx2[msk2]]['k_m_ext'],bins=linspace(5,16,50),label='match')
histstep(twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
figure()
histstep(twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
#WMAP9.distmod(noma6df['z_helio'])
histstep(WMAP9.distmod(noma6df['z_helio'])-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
histstep(WMAP9.distmod(noma6df['z_helio'][msk])-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
clf()
histstep(WMAP9.distmod(noma6df['z_helio'][msk])-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
histstep(WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
clf()
histstep(WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
histstep(WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
WMAP9.distmod(noma6df['z_helio'][msk])
histstep(WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext'],bins=linspace(5,16,50),label='no match')
WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext']
WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext']
histstep(twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value,bins=linspace(5,16,50),label='no match')
histstep(twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value,bins=linspace(5,16,50),label='no match')
histstep(twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value,bins=linspace(5,16,50),label='no match')
WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext']
twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value-twomassxsc[idx[msk]]['k_m_ext']
twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value
twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value
histstep(twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value,label='no match')
#v=histstep(twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value,label='no match')
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match')
clf()
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match')
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match')
clf()
histstep(twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value,label='no match')
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match')
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match',bins=linspace(-22,-4,5))
clf()
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match',bins=linspace(-22,-4,50))
v=twomassxsc[idx2[msk2]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk2]).value;histstep(v[~np.isnan(v)], label='match',bins=linspace(-22,-4,50))
v=twomassxsc[idx2[msk2]]['k_m_ext']-WMAP9.distmod(ma6df['z_helio'][msk2]).value;histstep(v[~np.isnan(v)], label='match',bins=linspace(-22,-4,50))
v=twomassxsc[idx2[msk2]]['k_m_ext']-WMAP9.distmod(ma6df['z_helio'][msk2]).value;histstep(v[~np.isnan(v)], label='match',bins=linspace(-22,-4,50))
clf()
v=twomassxsc[idx2[msk2]]['k_m_ext']-WMAP9.distmod(ma6df['z_helio'][msk2]).value;histstep(v[~np.isnan(v)], label='match',bins=linspace(-22,-4,50))
v=twomassxsc[idx[msk]]['k_m_ext']-WMAP9.distmod(noma6df['z_helio'][msk]).value;histstep(v[~np.isnan(v)], label='no match',bins=linspace(-22,-4,50))
legend(loc=0)
figure(2)
legend(loc=0)
figure(2)
xlabel('K')
figure(3)
xlabel('$M_K$')
get_ipython().system(u'ls -F -G ')
zc=fits.open('zcat.fits')
len(zc)
zc
zc[1].data
zc=fits.getdata('zcat.fits')
zc
len(zc)
zc['_RAJ2000']
zc=fits.getdata('zcat.fits')
len(zc)
zc['_RAJ2000']
get_ipython().system(u'ls RC3')
get_ipython().system(u'ls RC3/rc3')
get_ipython().system(u'ls RC3/rc3.py')
get_ipython().system(u'head RC3/rc3.py')
get_ipython().system(u'ls -F -G ')
get_ipython().system(u'ls -F -G ')
get_ipython().system(u'rm zcat.dat')
zc=fits.getdata('zcat.fits')
zc
zc['Vh']
zcat.dtype
zc.dtype
#res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=1*u.arcmin)
reload(masterlist)
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=1*u.arcmin)
res['zcat_nomatch]
res['zcat_nomatch']
res=masterlist.x_match_tests(mastercat1,vcuts=3000*u.km/u.s,tol=10*u.arcmin)
res['zcat_nomatch']
res['zcat_nomatch'].dtype
res['zcat_nomatch']['Bmag']
figure()
histstep(res['zcat_nomatch']['Bmag'])
v=res['zcat_nomatch']['Bmag'];;histstep(res['zcat_nomatch']['Bmag'])
v=res['zcat_nomatch']['Bmag'];histstep(res['zcat_nomatch']['Bmag'])
sum(msk)
len(msk)
#msk=dd.arcsec<2
get_ipython().system(u'ls -F -G ')
v=res['zcat_nomatch']['Bmag'];histstep(res['zcat_nomatch']['Bmag'])
v=res['zcat_nomatch']['Bmag'];histstep(res['zcat_nomatch']['Bmag'])
figure(1)
v=res['zcat_nomatch']['Bmag'];histstep(res['zcat_nomatch']['Bmag'])
v=res['zcat_nomatch']['Bmag'];histstep(v[~np.isnan(v)])
clf()
v=res['zcat_nomatch']['Bmag'];histstep(v[~np.isnan(v)])
v=res['zcat_nomatch']['Bmag'];histstep(v[~np.isnan(v)],linspace(4,20,50))
clf()
v=res['zcat_nomatch']['Bmag'];histstep(v[~np.isnan(v)],linspace(4,20,50))
v=res['zcat_match']['Bmag'];histstep(v[~np.isnan(v)],linspace(4,20,50))
clf()
v=res['zcat_match']['Bmag'];histstep(v[~np.isnan(v)],linspace(4,20,50),label='match')
v=res['zcat_nomatch']['Bmag'];histstep(v[~np.isnan(v)],linspace(4,20,50),label='no match')
legend(loc=0)
xlabel('B')
get_ipython().system(u'ls -F -G ')
figure(2)
title('masterlist vs. 6dF')
figure(3)
title('masterlist vs. 6dF')
title('masterlist vs. 6dF (abs)')
figure(1)
title('masterlist vs. zcat')
#res['zcat_nomatch']['Bmag']
savefig('zcatcomp.pdf')
figure(2)
savefig('6dfcomp.pdf')
figure(3)
savefig('6dfcompabs.pdf')
