# coding: utf-8
import sys
get_ipython().system(u'ls -F ')
get_ipython().magic(u'cd ~')
get_ipython().system(u'ls -F ')
get_ipython().magic(u'cd ')
get_ipython().system(u'ls -F ')
get_ipython().magic(u'cd ~/projects/distantlg')
get_ipython().system(u'ls -F ')
get_ipython().system(u'ls -F saga_catalogs/')
import hosts,targeting
hs=hosts.load_all_hosts(keyonname=True)
cnts={}
ds={}
for h in hs.values():
    print 'on',h.name
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    #tcat=targeting.select_targets(h,outercutrad=300)
    #cnts[h.name] = len(tcat)
cnts={}
ds={}
for h in hs.values():
    print 'on',h.name
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    #tcat=targeting.select_targets(h,outercutrad=300)
    #cnts[h.name] = len(tcat)
h
cnts={}
ds={}
for h in hs.values():
    print 'on',h.name
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300)
    cnts[h.name] = len(tcat)
h0=hs.values()[0]
h0._cached_sdss
h0._cached_sdss;
cnts={}
ds={}
for h in hs.values():
    print 'on',h.name
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300)
    cnts[h.name] = len(tcat)
    h._cached_sdss = None
cnts={}
ds={}
for i,h in enumerate(hs.values()):
    print 'on',h.name,'#',i+1,'of',len(hs)
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300)
    cnts[h.name] = len(tcat)
    h._cached_sdss = None
get_ipython().magic(u'pinfo targeting.select_targets')
cnts={}
ds={}
for i,h in enumerate(hs.values()):
    print 'on',h.name,'#',i+1,'of',len(hs)
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300,removegama=False)
    cnts[h.name] = len(tcat)
    h._cached_sdss = None
cnts
d
ds
x=np.array(ds.values())
y=np.array(cnts.values())
y=np.array([cnts[k] for k in ds])
x=np.array(ds.values())
y=np.array([cnts[k] for k in ds])
scatter(x,y)
semilogy()
#y=np.array([h. for h in hs.values()])
h
get_ipython().magic(u'pinfo h.physical_to_projected')
get_ipython().magic(u'pinfo h.projected_to_physical')
#y=np.array([h. for h in hs.values()])
get_ipython().magic(u'pinfo targeting.select_targets')
y=np.array([h.projected_to_physical(90) for h in hs.values()])
y
y=np.array([h.projected_to_physical(1.5) for h in hs.values()])
y
y=np.array([h.physical_to_projected(300) for h in hs.values()])
y
scatter(x,y)
clf()
scatter(x,y)
y=np.array([h.physical_to_projected(300) for h in hs.values()])
x=np.array([h.distmpc for h in hs.values()])
clf()
scatter(x,y)
y=np.array([cnts[hk] for hk in hs])
scatter(x,y)
clf()
scatter(x,y)
cat=h.get_sdss_catalog()
len(targeting.select_targets(h,removegama=False,outercutrad=300))
len(targeting.select_targets(h,removegama=False,outercutrad=-90))
len(targeting.select_targets(h,removegama=False,outercutrad=-300))
len(targeting.select_targets(h,removegama=False,outercutrad=-30))
len(targeting.select_targets(h,removegama=False,outercutrad=10))
len(targeting.select_targets(h,removegama=False,outercutrad=10))
reload(targeting)
len(targeting.select_targets(h,removegama=False,outercutrad=10))
len(targeting.select_targets(h,removegama=False,outercutrad=100))
len(targeting.select_targets(h,removegama=False,outercutrad=-100))
reload(targeting)
len(targeting.select_targets(h,removegama=False,outercutrad=-100))
len(targeting.select_targets(h,removegama=False,outercutrad=100))
len(targeting.select_targets(h,removegama=False,outercutrad=300))
len(targeting.select_targets(h,removegama=False,outercutrad=-90))
len(targeting.select_targets(h,removegama=False,outercutrad=300))
h.distmpc
degrees(.3/h.distmpc)
degrees(arcsin(.3/h.distmpc))
degrees(arcsin(1/h.distmpc))
degrees(arcsin(.3/h.distmpc))
cat
degrees(arcsin(.3/h.distmpc))
len(targeting.select_targets(h,removegama=False,outercutrad=300))
len(targeting.select_targets(h,removegama=False,outercutrad=300))
reload(targeting)
len(targeting.select_targets(h,removegama=False,outercutrad=300))
len(targeting.select_targets(h,removegama=False,outercutrad=350))
len(targeting.select_targets(h,removegama=False,outercutrad=300))
clf()
cnts={}
ds={}
for i,h in enumerate(hs.values()):
    print 'on',h.name,'#',i+1,'of',len(hs)
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300,removegama=False)
    cnts[h.name] = len(tcat)
    h._cached_sdss = None
cnts={}
ds={}
for i,h in enumerate(hs.values()):
    print 'on',h.name,'#',i+1,'of',len(hs)
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300,removegama=False)
    cnts[h.name] = len(tcat)
    h._cached_sdss = None
reload(targeting)
cnts={}
ds={}
for i,h in enumerate(hs.values()):
    print 'on',h.name,'#',i+1,'of',len(hs)
    sys.stdout.flush()
    ds[h.name] = h.distmpc
    h.fnsdss = 'saga_catalogs/'+h.name+'.dat'
    tcat=targeting.select_targets(h,outercutrad=300,removegama=False)
    cnts[h.name] = len(tcat)
    h._cached_sdss = None
scatter(x,y)
x=np.array([h.distmpc for h in hs.values()])
y=np.array([cnts[hk] for hk in hs])
scatter(x,y)
xlabel(r'$d_{\rm host}/{\rm Mpc}$')
ylabel(r'$N_{\rm sats}(<300 {\rm kpc})$')
tight_layout()
ylabel(r'$N_{\rm sats}(<300 {\rm kpc}) \times 10^3$')
ylabel(r'$N_{\rm sats}(<300 {\rm kpc}) / 10^3$')
ylabel(r'$N_{\rm sats}(<300 {\rm kpc}) / 10^3$')
clf()
scatter(x,y*10^-3)
scatter(x,y*10**-3)
clf()
scatter(x,y*10**-3)
xlabel(r'$d_{\rm host}/{\rm Mpc}$')
ylabel(r'$N_{\rm sats}(<300 {\rm kpc}) / 10^3$')
tight_layout()
xlim(10,62)
xlim(20,62)
xlim(18,62)
#savefig('/Users/erik/Desktop/Nvsd.pdf')
#savefig('/Users/erik/Desktop/Nvsd.png')
figure()
figure(1)
ylabel(r'$N_{\rm targets}(<300 {\rm kpc}) / 10^3$')
savefig('/Users/erik/Desktop/Nvsd.png')
savefig('/Users/erik/Desktop/Nvsd.pdf')
figure(2)
scatter(x,y*10**-3)
x2=np.array([h.physical_to_projected(300) for h in hs.values()])
scatter(x2,y*10**-3)
clf()
scatter(x2,y*10**-3)
xlabel(r'$R_{\rm vir,proj}/{\rm arcmin}$')
xlabel(r'$R_{\rm vir,proj}/{\rm arcmin}$')
ylabel(r'$N_{\rm targets}(<300 {\rm kpc}) / 10^3$')
tight_layout()
savefig('/Users/erik/Desktop/NvsR.pdf')
savefig('/Users/erik/Desktop/NvsR.png')
