# coding: utf-8
get_ipython().system(u'ls -F -G ')
f=open('EDDtable25Jan2014002749.csv.bk')
f.readline()
len(','.split(_))
f.readline()
len(_.split(','))
f=open('EDDtable25Jan2014002749.csv.bk')
l=f.readline();print l;len(','.split(l))
l=f.readline();print l;len(;.split(','))
l=f.readline();print l;len(l.split(','))
f.close()
f=open('EDDtable25Jan2014002749.csv.bk')
l=f.readline();print l;len(l.split(','))
l=f.readline();print l;len(l.split(','))
l=f.readline();print l;len(l.split(','))
l=f.readline();print l;len(l.split(','))
l=f.readline();print l;len(l.split(','))
l=f.readline();print l;len(l.split(','))
d=ascii.read('EDDtable27Jan2014053708.txt.csv',delimiter=',',Reader=ascii.Basic)
f=open('EDDtable27Jan2014053708.txt.csv')
[len(f.readline().split(',')) for i in range(5)]
[len(f.readline().split(',')) for i in range(5)]
f=open('EDDtable27Jan2014053708.txt.csv')
[len(f.readline().split(',')) for i in range(5)]
f=open('EDDtable27Jan2014054332.txt.csv')
len(_.split(','))
[len(f.readline().split(',')) for i in range(5)]
[len(f.readline().split(',')) for i in range(5)]
f=open('EDDtable27Jan2014054332.txt.csv')
f.readline()
f.readline()
[len(f.readline().split('|')) for i in range(5)]
f=open('EDDtable27Jan2014054332.txt.csv')
[len(f.readline().split('|')) for i in range(5)]
f=open('EDDtable27Jan2014054836.txt.csv')
f=open('EDDtable27Jan2014054836.csv')
[len(f.readline().split('|')) for i in range(5)]
[len(f.readline().split('|')) for i in range(5)]
d=ascii.read('EDDtable27Jan2014054836.csv',delimiter=',',Reader=ascii.Basic)
get_ipython().magic(u'pinfo votable.parse')
votable.
from astropy.io import votable
get_ipython().magic(u'pinfo votable.parse')
get_ipython().system(u'ls -F -G ')
d=votable.parse('EDDtable27Jan2014055034.votable')
t=d.get_first_table()
t
get_ipython().system(u'ls -F -G ')
d=ascii.read('EDDtable27Jan2014054836.csv',delimiter=',',Reader=ascii.Basic)
d=ascii.read('EDDtable27Jan2014055230.csv',delimiter=',',Reader=ascii.Basic)
d=ascii.read('EDDtable27Jan2014055230.csv',delimiter=',',Reader=ascii.Basic)
d=ascii.read('EDDtable27Jan2014055230.csv',delimiter=',',Reader=ascii.Basic)
d=ascii.read('EDDtable27Jan2014055230.csv',delimiter=',',Reader=ascii.Basic,guess=False)
d
d.dtype
get_ipython().system(u'ls -F -G ')
t
t
t.array
t.array.dtype
d.dtype
d['RA2000']
get_ipython().system(u'ls -F -G ')
d.dtype
get_ipython().system(u'ls -F -G ..')
f=fits.open('../nsa_v0_1_2.fits')
f
f.close()
n=fits.getdata('../nsa_v0_1_2.fits',1)
n
n.dtype
len(n)
b['RA']
n['RA"]
n['RA']
nc=coo.ICRS(u.deg*n['RA'],n['DEC']*u.deg)
d.dtype
edc=coo.ICRS(u.deg*d['RA2000'],d['DEC2000']*u.deg)
d['RA2000']
d.mask
d['RA2000'].mask
raa=d['RA2000'];deca=d['DEC2000']edc=coo.ICRS(u.deg*raa[~raa.mask],deca[~deca.mask]*u.deg)
raa=d['RA2000'];deca=d['DEC2000'];edc=coo.ICRS(u.deg*raa[~raa.mask],deca[~deca.mask]*u.deg)
u.deg*raa[~raa.mask]
u.deg*raa[~raa.mask]
deca[~deca.mask]*u.deg
deca[~deca.mask]
raa[~raa.mask]
u.deg*raa[~raa.mask]
u.deg*deca[~deca.mask]
raa=d['RA2000'];deca=d['DEC2000'];edc=coo.ICRS(u.deg*raa[~raa.mask],u.deg*deca[~deca.mask])
raa=d['RA2000'];deca=d['DEC2000'];edc=coo.ICRS(u.deg*raa[~raa.mask],u.deg*deca[~deca.mask])
raa=d['RA2000'];deca=d['DEC2000'];edc=coo.ICRS(u.deg*raa[~raa.mask],u.deg*deca[~deca.mask])
deca
raa
d=ascii.read('EDDtable27Jan2014060337.csv',delimiter=',',Reader=ascii.Basic,guess=False)
raa=d['RA2000'];deca=d['DEC2000'];edc=coo.ICRS(u.deg*raa[~raa.mask],u.deg*deca[~deca.mask])
d.dtype
raa=d['al2000'];deca=d['de2000'];edc=coo.ICRS(u.hour*raa[~raa.mask],u.deg*deca[~deca.mask])
edc
len(edc)
len(edc.ra)
len(nc.ra)
nc.ra
nc.dtype
nc
n.dtype
n['ZDIST']
n['ZDIST']
n['ZDIST']
n['ZDIST']*3e5
msk = n['ZDIST']*3e5 < 3000
sum(msk)
nc[msk]
nc[msk].ra.shape
ncm=nc[msk]
edc.dtype
get_ipython().magic(u'ed ')
d.dtype
d['v']
d['v'] < 3000
len(nc.ra)
len(edc)
nc[msk].ra.shape
len(edc.ra)
idx,dd,d3d=ncm.match_to_catalog_sky(edc)
dd
dd.arcmin
hist(dd.arcmin,bins=100)
hist(dd.arcmin,bins=100,histtype='step')
hist(dd.arcmin,bins=100,histtype='step',range(0,10))
clf()
hist(dd.arcmin,bins=100,histtype='step',range=(0,10))
hist(dd.arcmin,bins=100,histtype='step',range=(0,2))
clf()
hist(dd.arcmin,bins=100,histtype='step',range=(0,2))
hist(dd.arcmin,bins=100,histtype='step',range=(0,1))
clf()
hist(dd.arcmin,bins=100,histtype='step',range=(0,1),log=True)
sum(dd<.35)
sum(dd.arcmin<.35)
sum(dd.arcmin>.35)
sum(dd.arcmin>1)
sum(dd.arcmin>5)
sum(dd.arcmin>5)
nomtch=dd.arcmin>5
idx,dd,d3d=ncm.match_to_catalog_sky(edc)
ncm[nomtch]
#msk=dd10.arcmin<np.percentile(dd10.arcmin,10000000./len(d));dm=d[msk]
msk = n['ZDIST']*3e5 < 3000
n[msk][momtcj]
n[msk][momtch]
non=matchesn[msk][nomtch]
notm=n[msk][nomtch]
notm
notm['RA']
nc[momtcj]
nc[nomtch]
nc[nomtch][0]
lambda c:'http://www.nsatlas.org/getAtlas.html?search=radec&ra={0}&dec={1}&radius=10.0&submit_form=Submit'.format(c.ra.degree,c.dec.degree)
url=lambda c:'http://www.nsatlas.org/getAtlas.html?search=radec&ra={0}&dec={1}&radius=10.0&submit_form=Submit'.format(c.ra.degree,c.dec.degree)
url(nc[nomtch][0])
url(nc[nomtch][0])
import webbrowser
webbrowser.open(url(nc[nomtch][0]))
webbrowser.open(url(nc[nomtch][0]))
62*.7
62*.73
62*.7
webbrowser.open(url(nc[nomtch][5]))
webbrowser.open(url(nc[nomtch][15]))
webbrowser.open(url(nc[nomtch][55]))
#idx,dd,d3d=ncm.match_to_catalog_sky(edc)
webbrowser.open(url(ncm[nomtch][55]))
webbrowser.open(url(ncm[nomtch][0]))
webbrowser.open(url(ncm[nomtch][0]))
webbrowser.open(url(ncm[nomtch][1]))
webbrowser.open(url(ncm[nomtch][2]))
webbrowser.open(url(ncm[nomtch][8]))
ncm
n.dtype
n[msk][nomatch]
nnm=n[msk][nomtch]
nnm.dtype
len(nnm)
nnm['ABSMAG']
nnm['ABSMAG'].shape
nnm['ABSMAG'][-1]
hist(nnm['ABSMAG'][-1])
clf()
hist(nnm['ABSMAG'][-1],bins=100)
nmm=n[msk][!nomtch]
nmm=n[msk][~nomtch]
hist(nnm['ABSMAG'][-1],bins=100)
clf()
hist(nnm['ABSMAG'][-1],bins=100,histtype='step')
hist(nnm['ABSMAG'][-1],bins=100,histtype='step')
nnm['ABSMAG'][:,-1]
hist(nnm['ABSMAG'][:,-1],bins=100,histtype='step')
clf()
hist(nnm['ABSMAG'][:,-1],bins=100,histtype='step')
hist(nmm['ABSMAG'][:,-1],bins=100,histtype='step')
figure()
hist(nmm['ABSMAG'][:,-3],bins=100,histtype='step',label='In EDD')
hist(nnm['ABSMAG'][:,-3],bins=100,histtype='step',label='Not in EDD')
xlabel('$M_r$')
tight_layout()
xlim(-10,-23)
bins=linspace(-10,-23,101)
clf()
hist(nmm['ABSMAG'][:,-3],bins=bins,histtype='step',label='In EDD')
bins=linspace(-23,-10,101)
hist(nmm['ABSMAG'][:,-3],bins=bins,histtype='step',label='In EDD')
hist(nnm['ABSMAG'][:,-3],bins=bins,histtype='step',label='Not in EDD')
bins=linspace(-23,-10,51)
clf()
hist(nmm['ABSMAG'][:,-3],bins=bins,histtype='step',label='In EDD')
hist(nnm['ABSMAG'][:,-3],bins=bins,histtype='step',label='Not in EDD')
xlim(-10,-23)
xlabel('$M_r$')
tight_layout()
legend()
xlabel('NSA $M_r$')
xlabel('NSA $M_r$',textsize=36)
xlabel('NSA $M_r$',fontsize=36)
tight_layout()
xlabel('NSA $M_r$',fontsize=24)
tight_layout()
subplots_adjust(bottom=.1)
subplots_adjust(bottom=.3)
3000/70
subplots_adjust(bottom=.1)
xlim(-12,-23)
savefig('eddornot.pdf')
title('Things in NSA < 3000 km/s')
xlabel('$M_r$',fontsize=24)
title('Things in NSA < 3000 km/s',fontsize=24)
tight_layout()
savefig('eddornot.pdf')
3000/70
3000/73
3000/73.
3000/70.
get_ipython().system(u'ls -F -G ')
d2=ascii.read('ledaonly.csv',delimiter=',',Reader=ascii.Basic,guess=False)
d2
idx,dd,d3d=ncm.match_to_catalog_sky(edc)
lc=coo.ICRS(u.deg*d2['RA2000'],d2['DEC2000']*u.deg)
lc=coo.ICRS(u.deg*d2['al2000'],d2['dl2000']*u.deg)
lc=coo.ICRS(u.deg*d2['al2000'],d2['dl2000']*u.deg)
d2.dtype.names
lc=coo.ICRS(u.hour*d2['al2000'],d2['de2000']*u.deg)
lc=coo.ICRS(u.hour*d2['al2000'],d2['de2000']*u.deg)
#lc=coo.ICRS(u.hour*d2['al2000'],d2['de2000']*u.deg)
raa=d2['al2000'];deca=d2['de2000'];lc=coo.ICRS(u.hour*raa[~raa.mask],u.deg*deca[~deca.mask])
len(lc.ra)
coo.ICRS(46.183750,54.113333).match_to_catalog_sky(lc)
coo.ICRS(46.183750*u.deg,54.113333*u.deg).match_to_catalog_sky(lc)
coo.ICRS(50.720000*u.deg,58.653333*u.deg).match_to_catalog_sky(lc)
coo.ICRS(46.183750*u.deg,54.113333*u.deg).match_to_catalog_sky(lc)[1].arcmin
coo.ICRS(50.720000*u.deg,58.653333*u.deg).match_to_catalog_sky(lc)[1].arcmin
coo.ICRS(50.720000*u.deg,58.653333*u.deg).match_to_catalog_sky(edc)[1].arcmin
len(edc)
len(edc.ra)
len(lc.ra)
lc.dtype
d2.dtype
sum(d2['v']<3000)
nomtch=dd.arcmin>5
idx,dd,d3d=ncm.match_to_catalog_sky(lc)
nomtch=dd.arcmin>5
figure()
hist(nmm['ABSMAG'][:,-3],bins=bins,histtype='step',label='In LEDA')
nmm=n[msk][~nomtch]
nnm=n[msk][nomtch]
clf()
hist(nmm['ABSMAG'][:,-3],bins=bins,histtype='step',label='In LEDA')
hist(nnm['ABSMAG'][:,-3],bins=bins,histtype='step',label='Not in LEDA')
len(lc.ra)
clf()
hist(nmm['ABSMAG'][:,-3],bins=bins,histtype='step',label='In LEDA')
hist(nnm['ABSMAG'][:,-3],bins=bins,histtype='step',label='Not in LEDA')
xlabel('$M_r$',fontsize=24)
title('Things in NSA < 3000 km/s',fontsize=24)
legend()
tight_layout()
savefig('ledaornot.pdf')
figure(2)
legend()
legend(['a','b'])
legend(['In 2MRS','Not in 2MRS'])
legend(['In EDD's 2MRS','Not in EDD's 2MRS'])
legend(['In EDD\'s 2MRS','Not in EDD\'s 2MRS'])
savefig('eddornot.pdf')
savefig('mrsornot.pdf')
get_ipython().system(u'rm eddornot.pdf')
