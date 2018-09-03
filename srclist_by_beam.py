#!/usr/bin/env python
'''Script to create a single mega-patch calibrator based on the primary beam
convolved with the source fluxes'''
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
# from mwapy import ephem_utils

try:
    ##Andew's new mwa primary beam repo https://github.com/MWATelescope/mwa_pb.git
    from mwa_pb import primary_beam,beam_full_EE
except ImportError:
    ##Old MWA_tools version
    from mwapy.pb import primary_beam,beam_full_EE

from numpy import *
import sys
from optparse import OptionParser,OptionGroup
import matplotlib.pyplot as plt
import subprocess
##TO DO:
## - Add precess options?

from astropy.coordinates import Angle, SkyCoord, EarthLocation

parser = OptionParser()
parser.add_option('-x','--no_patch', action='store_true',
    help="Switch on to not put the sources in to a patch")
parser.add_option('-m','--metafits',
    help="metafits file to get obs info from")
parser.add_option('-s','--srclist',
    help="Base srclist to get source info from")
parser.add_option('-n','--num_sources',
    help="Number of sources to put in the mega patch")
parser.add_option('-p', '--plot',action='store_true',
    help='Plot the sources included - NEEDS THE MODULE wcsaxes')
parser.add_option('-c', '--cutoff', default=20.0,
    help='Distance from the pointing centre within which to accept source (cutoff in deg). Default is 20deg')
parser.add_option('-o', '--order', default='experimental',
    help='Criteria with which to order the output sources - "flux" for brightest first, "distance" for closest to pointing centre first, "experimental" for a combination, "name=*sourcename*" to force a calibrator. Default = "experimental". {To use default experimental, "experimental". Other, enter "experimental=flux,distance" with flux cutoff in Jy and distance cutoff in deg.} ')
parser.add_option('-a', '--outside',action='store_true', default=False,
    help='Switch on to only consider sources OUTSIDE of the cutoff, rather than inside')
parser.add_option('--aocalibrate',action='store_true', default=False,
    help='Output sourcelist in Andre Offringa format for use with CALIBRATE. Only works with the -x --no_patch option and does not work for shapelets (Gaussians OK)')

(options, args) = parser.parse_args()

cutoff = float(options.cutoff)

def hour_to_deg(time):      #converts hh:mm:ss.ss in to degrees, must input as a string
    negtest=time[0]
    time=time.split(':')
    degr=float(time[0])*15.0
    mins=float(time[1])*(15.0/60.0)
    secs=float(time[2])*(15.0/(3600.0))
    if negtest=='-':
        deg=degr-mins-secs
    if negtest!='-':
        deg=degr+mins+secs
    return deg

##Try opening the metafits file and complain if certain data is missing
try:
    f=pyfits.open(options.metafits)
except Exception,e:
    print 'Unable to open metafits file %s: %s' % (options.metafits,e)
    sys.exit(1)
if not 'DELAYS' in f[0].header.keys():
    print 'Cannot find DELAYS in %s' % options.metafits
    sys.exit(1)
if not 'LST' in f[0].header.keys():
    print 'Cannot find LST in %s' % options.metafits
    sys.exit(1)
if not 'GPSTIME' in f[0].header.keys():
    print 'Cannot find GPSTIME in %s' % options.metafits
    sys.exit(1)
if not 'FREQCENT' in f[0].header.keys():
    print 'Cannot find FREQCENT in %s' % options.metafits
    sys.exit(1)
if not 'RA' in f[0].header.keys():
    print 'Cannot find RA in %s' % options.metafits
    sys.exit(1)
if not 'DEC' in f[0].header.keys():
    print 'Cannot find DEC in %s' % options.metafits
    sys.exit(1)

##Gather the useful info
delays=array(map(int,f[0].header['DELAYS'].split(',')))
LST = float(f[0].header['LST'])
obsID = f[0].header['GPSTIME']
freqcent = f[0].header['FREQCENT']*1e+6
ra_point = f[0].header['RA']
dec_point = f[0].header['DEC']

##From the old MWA_tools ephem_utils============================================
##==============================================================================
def putrange(x,r=24):
    """ puts a value in the range [0,r)
    """

    if (not isinstance(x,ndarray)):
        while (x<0):
            x+=r
        while (x>=r):
            x-=r
        return x
    else:
        # numpy version
        while (any(x<0)):
            x[x<0]+=r
        while (any(x>=r)):
            x[x>=r]-=r
        return x

def eq2horz(HA, Dec, lat):
    """
    [Az,Alt]=eq2horz(HA,Dec,lat)
    equatorial to horizon coords
    all decimal degrees
    The sign convention for azimuth is north zero, east +pi/2.

    from slalib sla_e2h
    https://starlink.jach.hawaii.edu/viewvc/trunk/libraries/sla/e2h.f?revision=11739&view=markup
    https://starlink.jach.hawaii.edu/viewvc/trunk/libraries/sla/

    azimuth here is defined with N=0
    """

    if (isinstance(HA,ndarray)):
        sh=sin(HA*pi/180)
        ch=cos(HA*pi/180)
        sd=sin(Dec*pi/180)
        cd=cos(Dec*pi/180)
        sl=sin(lat*pi/180)
        cl=cos(lat*pi/180)

        # (Az,El) as (x,y,z)
        x=-ch*cd*sl+sd*cl
        y=-sh*cd
        z=ch*cd*cl+sd*sl

        # to spherical
        r=sqrt(x*x+y*y)
        a=arctan2(y,x)
        a[where(r==0)]=0
        a[where(a<0)]+=pi*2
        el=arctan2(z,r)
    else:
        sh=sin(HA*pi/180)
        ch=cos(HA*pi/180)
        sd=sin(Dec*pi/180)
        cd=cos(Dec*pi/180)
        sl=sin(lat*pi/180)
        cl=cos(lat*pi/180)

        # (Az,El) as (x,y,z)
        x=-ch*cd*sl+sd*cl
        y=-sh*cd
        z=ch*cd*cl+sd*sl

        # to spherical
        r=sqrt(x*x+y*y)
        if (r==0):
            a=0
        else:
            a=arctan2(y,x)
        a=putrange(a,2*pi)
        el=arctan2(z,r)

    return [a*180/pi, el*180/pi]

def dec2sex(x):
    """ convert decimal to sexadecimal
    note that this fails for -1<x<0: d will be 0 when it should be -0
    """

    sign=1
    if (x<0):
        sign=-1
    x=math.fabs(x)

    d=int(x)
    m=int(60*(x-d))
    s=60*(60*(x-d)-m)
    if (sign == -1):
        d*=-1

    return (d,m,s)

def dec2sexstring(x, includesign=0,digits=2,roundseconds=0):
    """
    dec2sexstring(x, includesign=0,digits=2,roundseconds=0)
    convert a decimal to a sexadecimal string
    if includesign=1, then always use a sign
    can specify number of digits on seconds (if digits>=0) or minutes (if < 0)
    """

    (d,m,s)=dec2sex(float(x))

    if (not roundseconds):
        sint=int(s)
        if (digits>0):
            sfrac=(10**digits)*(s-sint)
            ss2='%02' + 'd' + '.%0' + ('%d' % digits) + 'd'
            ss=ss2 % (sint,sfrac)
        elif (digits == 0):
            ss='%02d' % sint
        else:
            mfrac=10**(math.fabs(digits))*(s/60.0)
            ss2='%02' + 'd' + '.%0' + ('%d' % math.fabs(digits)) + 'd'
            ss=ss2 % (m,mfrac)
    else:
        sint=int(s)
        if (digits == 0):
            ss='%02d' % (round(s))
        elif (digits > 0):
            ss2='%02.' + ('%d' % digits) + 'f'
            ss=ss2 % s
            if (s < 10):
                ss='0' + ss
        else:
            ss2='%02.' + ('%d' % math.fabs(digits)) + 'f'
            ss=ss2 % (m+s/60.0)
            if (m < 10):
                ss='0' + ss


    if (not includesign):
        if (digits>=0):
            sout="%02d:%02d:%s" % (d,m,ss)
        else:
            sout="%02d:%s" % (d,ss)
        if (float(x)<0 and not sout.startswith("-")):
            sout='-' + sout
    else:
        sign='+'
        if (float(x)<0):
            sign='-'
        if (digits>=0):
            sout="%s%02d:%02d:%s" % (sign,math.fabs(d),m,ss)
        else:
            sout="%s%02d:%s" % (sign,math.fabs(d),ss)

    return sout

##==============================================================================

##Class to store source information with - set to lists
##to store component info in the same place
class rts_source():
    def __init__(self):
        self.name = ''
        self.ras = []
        self.has = []
        self.decs = []
        self.freqs = []
        self.fluxs = []
        self.extrap_fluxs = []
        self.weighted_flux = -1
        self.offset = -1
        self.shapelets = []
        self.shapelet_indexes = []
        self.shapelet_coeffs = []
        self.gaussians = []
        self.gaussian_indexes = []
        self.beam_inds = []

def extrap_flux(freqs,fluxs,extrap_freq):
    '''f1/f2 = (nu1/n2)**alpha
       alpha = ln(f1/f2) / ln(nu1/nu2)
       f1 = f2*(nu1/nu2)**alpha'''
    alpha = log(fluxs[0]/fluxs[1]) / log(freqs[0]/freqs[1])
    extrap_flux = fluxs[0]*(extrap_freq/freqs[0])**alpha
    return extrap_flux

def arcdist(RA1,RA2,Dec1,Dec2):
    '''calculates distance between two celestial points in degrees'''
    dr = pi/180.0
    in1 = (90.0 - Dec1)*dr
    in2 = (90.0 - Dec2)*dr
    RA_d = (RA1 - RA2)*dr
    cosalpha = cos(in1)*cos(in2) + sin(in1)*sin(in2)*cos(RA_d)
    alpha = arccos(cosalpha)
    return alpha/dr

##Find the mwa latitde for converting ha,dec to Az,Alt

MWAPOS = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
mwa_lat = MWAPOS.lat.deg

delays=repeat(reshape(delays,(1,16)),2,axis=0)

##Read in the srclist information
rts_srcs = open(options.srclist,'r').read().split('ENDSOURCE')
del rts_srcs[-1]


def create_source(prim_name=None, prim_ra=None, prim_dec=None, offset=None, primary_info=None,beam_ind=None):
    '''Takes the information for a source and puts it into an rts_source class for later use'''

    source = rts_source()
    source.name = prim_name
    source.ras.append(float(prim_ra))
    source.decs.append(float(prim_dec))
    all_ras.append(float(prim_ra))
    all_decs.append(float(prim_dec))
    source.beam_inds.append(beam_ind)
    beam_ind += 1

    source.offset = offset
    ##Find the fluxes and append to the source class
    prim_freqs = []
    prim_fluxs = []
    for line in primary_info:
        if 'FREQ' in line:
            ##Ignore if 0 flux
            if float(line.split()[2]) <= 0.0:
                pass
            else:
                prim_freqs.append(float(line.split()[1]))
                prim_fluxs.append(float(line.split()[2]))
    source.freqs.append(prim_freqs)
    source.fluxs.append(prim_fluxs)

    ##Split all info into lines and get rid of blank entries
    lines = split_source.split('\n')
    lines = [line for line in lines if line!='']
    ##If there are components to the source, see where the components start and end
    comp_starts = [lines.index(line) for line in lines if 'COMPONENT' in line and 'END' not in line]
    comp_ends = [i for i in xrange(len(lines)) if lines[i]=='ENDCOMPONENT']

    ##Check to see if the primary source is a gaussian or shapelet
    for line in primary_info:
        ###Check here to see if primary source is a gaussian:
        if 'GAUSSIAN' in line:
            source.gaussians.append(line)
            source.gaussian_indexes.append(0)
        ##As shapelet line comes after the
        elif 'SHAPELET' in line:
            coeffs = []
            source.shapelets.append(line)
            source.shapelet_indexes.append(0)
            ##If the primary source is a shapelet, search for shapelet coeffs in primary data,
            ##gather in a list and append to the source class
            for line in primary_info:
                if 'COEFF' in line: coeffs.append(line)
            source.shapelet_coeffs.append(coeffs)

    ##For each component, go through and find ra,dec,freqs and fluxs
    ##Also check here if the component is a gaussian or shapelet
    for start,end in zip(comp_starts,comp_ends):
        freqs = []
        fluxs = []
        coeffs = []
        for line in lines[start:end]:
            if 'COMPONENT' in line:
                source.ras.append(float(line.split()[1]))
                source.decs.append(float(line.split()[2]))
                all_ras.append(float(line.split()[1]))
                all_decs.append(float(line.split()[2]))
                source.beam_inds.append(beam_ind)
                beam_ind += 1

            elif 'FREQ' in line:
                freqs.append(float(line.split()[1]))
                fluxs.append(float(line.split()[2]))
            elif 'GAUSSIAN' in line:
                source.gaussians.append(line)
                gaus_ind = comp_starts.index(start) + 1
                source.gaussian_indexes.append(gaus_ind)
            elif 'SHAPELET' in line:
                source.shapelets.append(line)
                shap_ind = comp_starts.index(start) + 1
                source.shapelet_indexes.append(shap_ind)
            elif 'COEFF' in line:
                coeffs.append(line)
        source.fluxs.append(fluxs)
        source.freqs.append(freqs)
        if len(coeffs) > 0:
            source.shapelet_coeffs.append(coeffs)

    ##For each set of source infomation, calculate and extrapolated flux at the centra flux value
    for freqs,fluxs in zip(source.freqs,source.fluxs):

        fluxs = array(fluxs)[where(array(fluxs) > 0.0)]
        freqs = array(freqs)[where(array(fluxs) > 0.0)]

        ##If only one freq, extrap with an SI of -0.8:
        if len(freqs)==1:
            ##f1 = c*v1**-0.8
            ##c = f1 / (v1**-0.8)
            c = fluxs[0] / (freqs[0]**-0.8)
            ext_flux = c*freqcent**-0.8
        ##If extrapolating below known freqs, choose two lowest frequencies
        elif min(freqs)>freqcent:
            ext_flux = extrap_flux([freqs[0],freqs[1]],[fluxs[0],fluxs[1]],freqcent)
        ##If extrapolating above known freqs, choose two highest frequencies
        elif max(freqs)<freqcent:
            ext_flux = extrap_flux([freqs[-2],freqs[-1]],[fluxs[-2],fluxs[-1]],freqcent)
        ##Otherwise, choose the two frequencies above and below, and extrap between them
        else:
            for i in xrange(len(freqs)-1):
                if freqs[i]<freqcent and freqs[i+1]>freqcent:
                    ext_flux = extrap_flux([freqs[i],freqs[i+1]],[fluxs[i],fluxs[i+1]],freqcent)
        source.extrap_fluxs.append(ext_flux)

    source_weights = []

    sources.append(source)

    return beam_ind

def get_beam_weights(ras=None,decs=None):
    '''Takes ra and dec coords, and works out the overall beam power
    at that location using the 2016 spherical harmonic beam code from mwapy'''

    ##For each component, work out it's position, convolve with the beam and sum for the source
    #for ra,dec in zip(source.ras,source.decs):
    ##HA=LST-RA in def of ephem_utils.py
    has = LST - array(ras)*15.0  ##RTS stores things in hours
    ##Convert to zenith angle, azmuth in rad

    Az,Alt=eq2horz(has,array(decs),mwa_lat)
    za=(90-Alt)*pi/180
    az=Az*pi/180

    XX,YY = primary_beam.MWA_Tile_full_EE([za], [az], freq=freqcent, delays=delays, zenithnorm=True, power=True, interp=False)

    beam_weights = (XX[0]+YY[0]) / 2.0
    ##OLd way of combining XX and YY - end up with beam values greater than 1, not good!
    #beam_weights = sqrt(XX[0]**2+YY[0]**2)

    return beam_weights

sources = []
all_ras = []
all_decs = []

beam_ind = 0
##Go through all sources in the source list, gather their information, extrapolate
##the flux to the central frequency and weight by the beam at that position
for split_source in rts_srcs:

    source = rts_source()

    ##Find the primary source info - even if no comps, this will isolate
    ##primary source infomation
    primary_info = split_source.split('COMPONENT')[0].split('\n')
    primary_info = [info for info in primary_info if info!='']

    meh,prim_name,prim_ra,prim_dec = primary_info[0].split()

    ##Check if the primary RA,Dec is below the horizon - it will crash the RTS otherwise
    ##Skip if so
    ha_prim = LST - float(prim_ra)*15.0
    Az_prim,Alt_prim = eq2horz(ha_prim,float(prim_dec),mwa_lat)

    if Alt_prim < 0.0:
        pass
    else:
        offset = arcdist(float(ra_point),float(prim_ra)*15.0,float(dec_point),float(prim_dec))
        if options.outside:
            if offset > cutoff:
                beam_ind = create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset, primary_info=primary_info,beam_ind=beam_ind)
            else:
                pass
        else:
            if offset <= cutoff:
                beam_ind = create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset, primary_info=primary_info,beam_ind=beam_ind)
            else:
                pass

##Need to work out all beam weightings in one single calculation,
##as in each instance of the beam ~40s to run
beam_weights = get_beam_weights(ras=all_ras,decs=all_decs)

##Go through all the sources, and apply the beam weights to all
##components in the source. Dot the weights and fluxes to get
##a total weighted flux
for source in sources:
    source_weights = beam_weights[source.beam_inds]
    source.weighted_flux = dot(array(source_weights),array(source.extrap_fluxs))

##Make a list of all of the weighted_fluxes and then order the sources according to those
all_weighted_fluxs = [source.weighted_flux for source in sources]
weighted_sources = [source for flux,source in sorted(zip(all_weighted_fluxs,sources),key=lambda pair: pair[0],reverse=True)][:int(options.num_sources)]

##Apparently the RTS dies if the top source is a gaussian so do a check
##Move the gaussian out of the top spot
##Actually can't do it here - brightest source is found later in "experimental"
#if len(weighted_sources[0].gaussians):
    #print 'Base source is a gaussian - RTS does not like this'
    #move_ind = 1
    #while len(weighted_sources[move_ind].gaussians) > 0:
        #move_ind += 1
    #print 'Moving source from position %d to the top' %(move_ind+1)
    #new_weighted_sources = [weighted_sources[move_ind]]
    #for index,source in enumerate(weighted_sources):
        #if index == move_ind:
            #pass
        #else:
            #new_weighted_sources.append(source)
    #weighted_sources = new_weighted_sources

if options.no_patch:
    if not options.aocalibrate:
        print "++++++++++++++++++++++++++++++++++++++\nCreating weighted srclist - not mega-patching the sources"
        if options.outside:
            output_name = "%s_%s_outside-cutoff_peel%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
        else:
            output_name = "%s_%s_peel%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
        out_file = open(output_name,'w+')
        for source in weighted_sources[:1]:
            out_file.write('SOURCE %s %.10f %.10f' %(source.name,source.ras[0],source.decs[0]))
            ##If base source was a gaussian put it in:
            if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
                out_file.write('\n'+source.gaussians[0])
            ##write out the fluxes
            for flux,freq in zip(source.fluxs[0],source.freqs[0]):
                out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
            ##Shapelet coeffs come after the frequencies
            ##If base source was a shapelet, put in and it's coefficients in
            #print source.shapelets,source.shapelet_indexes,source.shapelet_coeffs
            if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
                #print 'here',source.shapelets[0]
                out_file.write('\n'+source.shapelets[0])
                for coeff in source.shapelet_coeffs[0]:
                    out_file.write('\n'+coeff)
            ##Cycle through any components in that primary calibator
            for i in range(1,len(source.ras)):
                out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
                ##Cycle through gaussian IDs and add the gaussian line if applicable
                for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
                    if gaus_ind == i: out_file.write('\n'+gaus_line)
                ##add the fluxes
                for flux,freq in zip(source.fluxs[i],source.freqs[i]):
                    out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
                ##Have a look and see if there is a shapelet in this source
                ##and write out if applicable
                try:
                    shap_ind = source.shapelet_indexes.index(i)
                    out_file.write('\n'+source.shapelets[shap_ind])
                    for line in source.shapelet_coeffs[shap_ind]:
                        out_file.write('\n'+line)
                except:
                    pass
                out_file.write('\nENDCOMPONENT')
        out_file.write('\nENDSOURCE')

        for source in weighted_sources[1:int(options.num_sources)]:
            out_file.write('\nSOURCE %s %.10f %.10f' %(source.name,source.ras[0],source.decs[0]))
            ##If base source was a gaussian put it in:
            if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
                out_file.write('\n'+source.gaussians[0])
            ##write out the fluxes
            for flux,freq in zip(source.fluxs[0],source.freqs[0]):
                out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
            ##Shapelet coeffs come after the frequencies
            ##If base source was a shapelet, put in and it's coefficients in
            if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
                out_file.write('\n'+source.shapelets[0])
                for coeff in source.shapelet_coeffs[0]:
                    out_file.write('\n'+coeff)
            ##Cycle through any components in that primary calibator
            for i in range(1,len(source.ras)):
                out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
                ##Cycle through gaussian IDs and add the gaussian line if applicable
                for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
                    if gaus_ind == i: out_file.write('\n'+gaus_line)
                ##add the fluxes
                for flux,freq in zip(source.fluxs[i],source.freqs[i]):
                    out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
                ##Have a look and see if there is a shapelet in this source
                ##and write out if applicable
                try:
                    shap_ind = source.shapelet_indexes.index(i)
                    out_file.write('\n'+source.shapelets[shap_ind])
                    for line in source.shapelet_coeffs[shap_ind]:
                        out_file.write('\n'+line)
                except:
                    pass
                out_file.write('\nENDCOMPONENT')
            out_file.write('\nENDSOURCE')
        out_file.close()
    #If you are doing aocalibrate
    else:
        if options.outside:
            print 'outside option not yet supported for aocalibrate - NO SOURCELIST WRITTEN'
            sys.exit()
        else:
            output_name = "%s_%s_aocal%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
        print "++++++++++++++++++++++++++++++++++++++\nCreating weighted srclist for AO calibrate - not mega-patching the sources"
        out_file = open(output_name,'w+')
        out_file.write('skymodel fileformat 1.1\n')
        for source in weighted_sources:
            ##If it is a shapelet then send a warning that it is being skipped and continue to next source
            if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
                print "WARNING: Source %s with weighted flux density %s Jy not included in model. Baselist contains a shapelet but --aocalibrate not compatible with shapelets." % (source.name,source.weighted_flux)
                continue
            out_file.write('source {\n')
            out_file.write('  name "%s"\n' % (source.name))
            #cycle through the components:
            for i in range(0,len(source.ras)):
                position_string_ra_hrs=str(source.ras[i]).split('.')[0]
                position_string_ra_remainder='0.'+str(source.ras[i]).split('.')[1]
                position_string_ra_remainder_sex=dec2sexstring(float(position_string_ra_remainder),includesign=0,digits=1,roundseconds=1)
                position_string_ra_dms='%sh%sm%ss' % (position_string_ra_hrs,position_string_ra_remainder_sex.split(':')[1],position_string_ra_remainder_sex.split(':')[2])
                position_string_dec=dec2sexstring(source.decs[i],includesign=0,digits=0,roundseconds=1)
                position_string_dec_dms='%sd%sm%ss' % (position_string_dec.split(':')[0],position_string_dec.split(':')[1],position_string_dec.split(':')[2])
                position_string = '%s %s' % (position_string_ra_dms,position_string_dec_dms)
                out_file.write('  component {\n')
                ##If source is a gaussian put in type gaussian:
                if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
                    out_file.write('    type gaussian\n')
                    out_file.write('    position %s\n' % (position_string))
                    gauss_string=source.gaussians[i]
                    gauss_maj=gauss_string.split()[2]
                    gauss_min=gauss_string.split()[3]
                    gauss_pa=gauss_string.split()[1]
                    out_file.write('    shape %s %s %s \n' % (gauss_maj,gauss_min,gauss_pa)) #(maj, min, pa)))
                    #out_file.write('\n'+source.gaussians[0])
                #else it is a point source
                else:
                    out_file.write('    type point\n')
                    out_file.write('    position %s\n' % (position_string))
                ##write out the fluxes
                for flux,freq in zip(source.fluxs[i],source.freqs[i]):
                    freq_MHz=int(freq/1000000.0)
                    out_file.write("    measurement {\n")
                    out_file.write("      frequency %i MHz\n" %(freq_MHz))
                    out_file.write("      fluxdensity Jy %.2f 0 0 0\n"%(flux))
                    out_file.write("    }\n")
                out_file.write("  }\n")
            out_file.write("}\n")
        out_file.close()

else:
    if options.aocalibrate:
        print 'Mega patching not supported for aocalibrate - NO SOURCELIST WRITTEN'
        sys.exit()

    elif options.order=='flux':
        ordered_sources = weighted_sources

    elif options.order=='distance':
        ordered_offsets = [source.offset for source in weighted_sources]
        ordered_sources = [source for offset,source in sorted(zip(ordered_offsets,weighted_sources),key=lambda pair: pair[0])]


    elif 'name' in options.order:
        name = options.order.split("=")[1]
        top_source_ind = [source.name for source in weighted_sources].index(name)
        top_source = weighted_sources[top_source_ind]
        ordered_sources = [top_source]
        for i in xrange(len(weighted_sources)):
            if i!=top_source_ind:
                ordered_sources.append(weighted_sources[i])
        print '++++++++++++++++++++++++++++++++++++++\nBase Source forced as %s with \nconvolved flux %.1fJy at a distance %.2fdeg\n---------------------------------' %(top_source.name,top_source.weighted_flux,top_source.offset)

    elif 'experimental' in options.order:
        #print 'here, here, here'


        if len(options.order.split('='))==1:
            flux_cut,dist_cut = 10.0,1.0
        else:
            flux_cut,dist_cut = options.order.split('=')[1].split(',')

        close_fluxs = []
        close_dists = []
        dist_cut = float(dist_cut)
        dist_cut_lower = 0.0

        ##Try to find all sources within distance cutoff above flux threshold - if none exist, extend search
        ##radii by 0.5 deg
        while len(close_fluxs)==0:
            for index in xrange(len(weighted_sources)):
                src = weighted_sources[index]
                flux = src.weighted_flux
                offset = src.offset
                gauss_len = len(src.gaussians)
                ##RTS does not like the top source being a gaussian
                ##Skip the gaussians - (13/11/2017)
                if dist_cut_lower<offset<dist_cut and gauss_len > 0:
                    print 'Potential primary calibrator is a gaussian - skipping (convolved flux %.1f)' %flux

                if flux>float(flux_cut) and dist_cut_lower<offset<dist_cut and gauss_len < 1:
                    close_fluxs.append(flux)
                    close_dists.append(offset)
            #print 'No primary calibrator between %.1fdeg and %.1fdeg of centre' %(dist_cut_lower,dist_cut)
            dist_cut+=0.5
            dist_cut_lower = dist_cut - 0.5
            if dist_cut_lower > cutoff:
                print "++++++++++++++++++++++++++++++++++++++\nNo source above %.2fJy within initial cutoff distance\nNO SOURCE LIST GENERATED\n++++++++++++++++++++++++++++++++++++++" %float(flux_cut)
                sys.exit()

        ##This is the brightest source within the base source cutoff distance
        brightest_close_flux = sorted(close_fluxs,reverse=True)[0]

        ##This is where the bright close source appears in the beam weighted source list
        weighted_fluxes = [source.weighted_flux for source in weighted_sources]
        brightest_ind = weighted_fluxes.index(brightest_close_flux)

        ##Use the positional offset as an identifier as well in case so other source has the
        ##same flux
        weighted_offsets = [source.offset for source in weighted_sources]
        brightest_close_offset = weighted_offsets[brightest_ind]

        ##Find out name of source
        weighted_names = [source.name for source in weighted_sources]
        brightest_close_name = weighted_names[brightest_ind]

        print "++++++++++++++++++++++++++++++++++++++\nBase source %s convolved flux is %.3fJy at a distance \nof %.3fdeg from point centre\n---------------------------------" %(brightest_close_name,brightest_close_flux,brightest_close_offset)

        ##Put this source at the top of the ordered list, and then append all other sources after
        ##NOTE - this means that apart from the top source, all other sources are flux ordered.
        ordered_sources = [weighted_sources[brightest_ind]]
        for source in weighted_sources:
            if source.weighted_flux!=brightest_close_flux and source.offset!=brightest_close_offset:
                ordered_sources.append(source)

    ##Make a new single patch source based on the user specified number of components
    if options.outside:
        output_name = "%s_%s_outside-cutoff_patch%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
    else:
        output_name = "%s_%s_patch%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
    out_file = open(output_name,'w+')

    ##Print out the strongest source as the primary calibator
    for source in ordered_sources[:1]:
        out_file.write('SOURCE %s %.10f %.10f' %(obsID,source.ras[0],source.decs[0]))
        ##If base source was a gaussian put it in:
        if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
            out_file.write('\n'+source.gaussians[0])
        ##write out the fluxes
        for flux,freq in zip(source.fluxs[0],source.freqs[0]):
            out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
        ##Shapelet coeffs come after the frequencies
        ##If base source was a shapelet, put in and it's coefficients in
        if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
            out_file.write('\n'+source.shapelets[0])
            for coeff in source.shapelet_coeffs[0]:
                out_file.write('\n'+coeff)
        ##Cycle through any components in that primary calibator
        for i in range(1,len(source.ras)):
            out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
            ##Cycle through gaussian IDs and add the gaussian line if applicable
            for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
                if gaus_ind == i: out_file.write('\n'+gaus_line)
            ##add the fluxes
            for flux,freq in zip(source.fluxs[i],source.freqs[i]):
                out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
            ##Have a look and see if there is a shapelet in this source
            ##and write out if applicable
            try:
                shap_ind = source.shapelet_indexes.index(i)
                out_file.write('\n'+source.shapelets[shap_ind])
                for line in source.shapelet_coeffs[shap_ind]:
                    out_file.write('\n'+line)
            except:
                pass
            out_file.write('\nENDCOMPONENT')

    ##For all other sources, add all information as COMPONENTS
    for source in ordered_sources[1:int(options.num_sources)]:
        for i in xrange(len(source.ras)):
            out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
            ##Cycle through gaussian IDs and add the gaussian line if applicable
            for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
                if gaus_ind == i: out_file.write('\n'+gaus_line)
            #Add the fluxes
            for flux,freq in zip(source.fluxs[i],source.freqs[i]):
                out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
            ##Have a look and see if there is a shapelet in this source
            ##and write out if applicable
            try:
                shap_ind = source.shapelet_indexes.index(i)
                out_file.write('\n'+source.shapelets[shap_ind])
                for line in source.shapelet_coeffs[shap_ind]:
                    out_file.write('\n'+line)
            except:
                pass
            out_file.write('\nENDCOMPONENT')
    out_file.write('\nENDSOURCE')
    out_file.close()

print "Created %s\n++++++++++++++++++++++++++++++++++++++" %output_name
##Finito!!

if options.plot:
    cmd = "./plot_srclist.py -m %s -s %s" %(options.metafits, output_name)
    subprocess.call(cmd,shell=True)
