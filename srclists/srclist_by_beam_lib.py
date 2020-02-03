'''Library containing functions to generate reduced RTS srclists, based on the
expected brightest calibrators when convovled with MWA primary beam'''
from __future__ import print_function
import sys
try:
    ##Andew's new mwa primary beam repo https://github.com/MWATelescope/mwa_pb.git
    from mwa_pb import primary_beam
except ImportError:
    # Old MWA_tools primary beam code
    try:
        from mwapy.pb import primary_beam
        print("WARNING: Using mwapy for MWA primary beam code.\nYou should use mwa_pb from https://github.com/MWATelescope/mwa_pb", file=sys.stderr)
    except ImportError:
        print("Could not import MWA primary beam code from mwa_pb or mwapy.\nYou may need to install mwa_pb from https://github.com/MWATelescope/mwa_pb", file=sys.stderr)
        exit(1)

from numpy import *
# from astropy.coordinates import EarthLocation

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

    from erfa eraHd2ae
    https://github.com/liberfa/erfa/blob/master/src/hd2ae.c

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

def create_source(prim_name=None, prim_ra=None, prim_dec=None, offset=None, primary_info=None,beam_ind=None,all_ras=False,all_decs=False, split_source=None,
  freqcent=None):
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
            if float(line.split()[2]) == 0.0:
                pass
            else:
                prim_freqs.append(float(line.split()[1]))
                prim_fluxs.append(float(line.split()[2]))
    source.freqs.append(prim_freqs)
    source.fluxs.append(prim_fluxs)

    ##Split all info into lines and get rid of blank entries
    lines = split_source.split('\n')
    lines = [line for line in lines if line!='' and '#' not in line]

    ##If there are components to the source, see where the components start and end
    comp_starts = [i for i in arange(len(lines)) if 'COMPONENT' in lines[i] and 'END' not in lines[i]]
    # comp_starts = where(array(comp_starts) == 'COMPONENT')
    comp_ends = [i for i in arange(len(lines)) if lines[i]=='ENDCOMPONENT']

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
        fluxs = array(fluxs)[where(array(fluxs) != 0.0)]
        freqs = array(freqs)[where(array(fluxs) != 0.0)]

        ##If only one freq, extrap with an SI of -0.8:
        if len(freqs)==1:
            ##f1 = c*v1**-0.8
            ##c = f1 / (v1**-0.8)
            c = fluxs[0] / (freqs[0]**-0.8)
            ext_flux = c*freqcent**-0.8
        ##If extrapolating below known freqs, choose two lowest frequencies
        elif min(freqs) >= freqcent:
            ext_flux = extrap_flux([freqs[0],freqs[1]],[fluxs[0],fluxs[1]],freqcent)
        ##If extrapolating above known freqs, choose two highest frequencies
        elif max(freqs) <= freqcent:
            ext_flux = extrap_flux([freqs[-2],freqs[-1]],[fluxs[-2],fluxs[-1]],freqcent)
        ##Otherwise, choose the two frequencies above and below, and extrap between them
        else:
            for i in arange(len(freqs)-1):
                if freqs[i]<freqcent and freqs[i+1]>freqcent:
                    ext_flux = extrap_flux([freqs[i],freqs[i+1]],[fluxs[i],fluxs[i+1]],freqcent)
        source.extrap_fluxs.append(ext_flux)

    source_weights = []

    return beam_ind, source, all_ras, all_decs

def get_beam_weights(ras=None,decs=None,LST=None, mwa_lat=None, freqcent=None,delays=None):
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

def write_rts_peel(weighted_sources=None,num_sources=None,output_name=None):
    '''Writes a simple peel srclist, where all the calibrators remain as
    separate sources'''
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

    for source in weighted_sources[1:num_sources]:
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

def write_rts_patch(ordered_sources=None,num_sources=None,output_name=None,
                    zeroJy_source=None,extra_calibrators=None,obsID=None):
    '''Takes the apparent flux ordered sources and writes them to a patch
       Adds a zeroJy_source as the base source if required, and adds extra
       dd calibrators after the patch if required.'''
    out_file = open(output_name,'w+')

    ##Add the strongest source as the base source
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
    ##If we added a zero Jy base source, need to include an extra source
    if zeroJy_source:
        all_sources = int(num_sources) + 1
    else:
        all_sources = int(num_sources)

    for source in ordered_sources[1:all_sources]:
        for i in arange(len(source.ras)):
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

    ##Add the extra calibrators here if required
    if extra_calibrators:
        print('Adding %d extra dd calibrators' %(len(extra_calibrators)))
        for source in extra_calibrators:

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
        print('---------------------------------')
    out_file.close()

def write_aocal(weighted_sources=None,output_name=None):
    out_file = open(output_name,'w+')
    out_file.write('skymodel fileformat 1.1\n')
    for source in weighted_sources:
        ##If it is a shapelet then send a warning that it is being skipped and continue to next source
        if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
            print("WARNING: Source %s with weighted flux density %s Jy not included in model. Baselist contains a shapelet but --aocalibrate not compatible with shapelets." % (source.name,source.weighted_flux))
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
