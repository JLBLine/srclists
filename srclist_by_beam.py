#!/usr/bin/env python
'''Script to generate reduced RTS srclists, based on the expected brightest
calibrators when convovled with MWA primary beam'''
from __future__ import print_function
try:
    import astropy.io.fits as pyfits
except:
    import pyfits

from numpy import *
import sys
import argparse
from srclist_by_beam_lib import *
from astropy.coordinates import EarthLocation

parser = argparse.ArgumentParser(description='Make RTS style srclists to either calibrate in a patch or peel with separate sources.')
parser.add_argument('-x','--no_patch', action='store_true',
    help="Switch on to not put the sources in to a patch")
parser.add_argument('-m','--metafits',
    help="metafits file to get obs info from")
parser.add_argument('-s','--srclist',
    help="Base srclist to get source info from")
parser.add_argument('-n','--num_sources',
    help="Number of sources to put in the mega patch")
parser.add_argument('-c', '--cutoff', default=20.0,
    help='Distance from the pointing centre within which to accept source (cutoff in deg). Default is 20deg')
parser.add_argument('-o', '--order', default='experimental',
    help='Criteria with which to order the output sources - defaults to --order=experimental=10,1.0 \
    Options are: "flux" for brightest first, "distance" for closest to pointing centre first, \
    "name=*sourcename*" to force a calibrator, "experimental=flux,distance" to search for a source \
    above a flux cutoff in Jy and within a distance cutoff in deg.} ')
parser.add_argument('-a', '--outside',action='store_true', default=False,
    help='Switch on to only consider sources OUTSIDE of the cutoff, rather than inside')
parser.add_argument('--aocalibrate',action='store_true', default=False,
    help='Output sourcelist in Andre Offringa format for use with CALIBRATE. Only works with the -x \
    --no_patch option and does not work for shapelets (Gaussians OK)')
parser.add_argument('-z', '--zeroJy_source',action='store_true', default=False,
    help='Add to include a zero Jy point source as the base source, located at the beam pointing centre')
parser.add_argument('-l', '--dd_cal_list', default=False,
    help='Add a specific list of sources to add as separate direction dependent calibrators, alongside the patch. \
    For example, if you want to calibrate 3 sources: a 1000 source patch; Fornax A (ben_ForA_lobe2_); \
    Pictor A (PICA_point), add --dd_cal_list=ben_ForA_lobe2_,PICA_point. \
    Overall you will have 1002 calibrators in your sky model.')
parser.add_argument('-b', '--num_extra_dd_cals', default=False,
    help='Include a number of separate direction dependent calibrators, alongside the patch, \
    from the n brightest calibrators after those included in the patch. For example, \
    --num_extra_dd_cals=10 and --num_sources=1000 when creating a patch will put sources \
    ranked 1-1000 in a patch, and add sources ranked 1001-1010 as separate dd calibrators. \
    Note this is additive to --dd_cal_list, so if you have two sources in --dd_cal_list, \
    --num_extra_dd_cals=10 and --num_sources=1000, you final sky model will have 1012 calibrators total')

args = parser.parse_args()

cutoff = float(args.cutoff)

##Try opening the metafits file and complain if certain data is missing
try:
    f=pyfits.open(args.metafits)
except:
    sys.exit('Unable to open metafits file %s' % (args.metafits))
if not 'DELAYS' in f[0].header.keys():
    sys.exit('Cannot find DELAYS in %s' % args.metafits)
if not 'LST' in f[0].header.keys():
    sys.exit('Cannot find LST in %s' % args.metafits)
if not 'GPSTIME' in f[0].header.keys():
    sys.exit('Cannot find GPSTIME in %s' % args.metafits)
if not 'FREQCENT' in f[0].header.keys():
    sys.exit('Cannot find FREQCENT in %s' % args.metafits)
if not 'RA' in f[0].header.keys():
    sys.exit('Cannot find RA in %s' % args.metafits)
if not 'DEC' in f[0].header.keys():
    sys.exit('Cannot find DEC in %s' % args.metafits)

##Gather the useful info
delays = array(f[0].header['DELAYS'].split(','), dtype=int)

if len(where(delays == 32)[0]) == 16:
    print('---------------------------------------------------------------------')
    print('srclist_by_beam.py: All beam delays in metafits are equal to 32.\nmwa_pb will fall over and the world will end, cannot create srclist\nExiting now')
    print('---------------------------------------------------------------------')
    exit()

LST = float(f[0].header['LST'])
obsID = f[0].header['GPSTIME']
freqcent = f[0].header['FREQCENT']*1e+6
ra_point = f[0].header['RA']
dec_point = f[0].header['DEC']

##Find the mwa latitde for converting ha,dec to Az,Alt

MWAPOS = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
mwa_lat = MWAPOS.lat.deg

delays=repeat(reshape(delays,(1,16)),2,axis=0)

##Read in the srclist information
rts_srcs = open(args.srclist,'r').read().split('ENDSOURCE')
del rts_srcs[-1]

##Here the user wants a number of separate calibrators alongside
##the patch - pull them out before ordering for the patch
if args.dd_cal_list:
    ##Get the names of the dd_cal_list
    dd_cal_list = args.dd_cal_list.split(',')
    extra_calibrators = []
    extra_ras = []
    extra_decs = []
else:
    dd_cal_list = []
    extra_calibrators = False

sources = []
all_ras = []
all_decs = []

beam_ind = 0
##Go through all sources in the source list, gather their information, extrapolate
##the flux to the central frequency
for split_source in rts_srcs:
    ##Find the primary source info - even if no comps, this will isolate
    ##primary source infomation
    primary_info = split_source.split('COMPONENT')[0].split('\n')
    primary_info = [info for info in primary_info if info!='']

    _,prim_name,prim_ra,prim_dec = primary_info[0].split()

    ##Check if the primary RA,Dec is below the horizon - it will crash the RTS otherwise
    ##Skip if so
    ha_prim = LST - float(prim_ra)*15.0
    Az_prim,Alt_prim = eq2horz(ha_prim,float(prim_dec),mwa_lat)

    if prim_name in dd_cal_list:
        _,source,extra_ras,extra_decs = create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset,
                                        primary_info=primary_info, beam_ind=0,all_ras=extra_ras,all_decs=extra_decs,
                                        split_source=split_source,freqcent=freqcent)

        extra_calibrators.append(source)

    else:
        if Alt_prim < 0.0:
            pass
        else:
            offset = arcdist(float(ra_point),float(prim_ra)*15.0,float(dec_point),float(prim_dec))
            if args.outside:
                if offset > cutoff:
                    beam_ind,source, all_ras, all_decs = create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset,
                                                         primary_info=primary_info,beam_ind=beam_ind,all_ras=all_ras,all_decs=all_decs,
                                                         split_source=split_source,freqcent=freqcent)
                    sources.append(source)
                else:
                    pass
            else:
                if offset <= cutoff:
                    beam_ind,source, all_ras, all_decs = create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset,
                                                         primary_info=primary_info,beam_ind=beam_ind,all_ras=all_ras,all_decs=all_decs,
                                                         split_source=split_source,freqcent=freqcent)
                    sources.append(source)
                else:
                    pass

##Need to work out all beam weightings in one single calculation,
##as in each instance of the beam ~40s to run
beam_weights = get_beam_weights(ras=all_ras,decs=all_decs,LST=LST,mwa_lat=mwa_lat,
                                freqcent=freqcent,delays=delays)

##Go through all the sources, and apply the beam weights to all
##components in the source. Dot the weights and fluxes to get
##a total weighted flux
for source in sources:
    source_weights = beam_weights[source.beam_inds]
    source.weighted_flux = dot(array(source_weights),array(source.extrap_fluxs))

##Make a list of all of the weighted_fluxes and then order the sources according to those
all_weighted_fluxs = [source.weighted_flux for source in sources]
all_sorted_sources = [source for flux,source in sorted(zip(all_weighted_fluxs,sources),key=lambda pair: pair[0],reverse=True)]

weighted_sources = all_sorted_sources[:int(args.num_sources)]

if args.num_extra_dd_cals:
    try:
        num_extra_dd_cals = int(args.num_extra_dd_cals)
        if not extra_calibrators:
            extra_calibrators = []

        for source in all_sorted_sources[int(args.num_sources):int(args.num_sources)+num_extra_dd_cals]:
            extra_calibrators.append(source)
    except ValueError:
        print('WARNING Failed to convert --num_extra_dd_cals=%s into an integer. Unable to add extra dd calibrators sensibly' %args.num_extra_dd_cals)

if args.no_patch:
    if not args.aocalibrate:
        print("++++++++++++++++++++++++++++++++++++++\nCreating weighted srclist - not mega-patching the sources")
        if args.outside:
            output_name = "%s_%s_outside-cutoff_peel%s.txt" %(args.srclist.split('/')[-1].split('.')[0],obsID,args.num_sources)
        else:
            output_name = "%s_%s_peel%s.txt" %(args.srclist.split('/')[-1].split('.')[0],obsID,args.num_sources)

        write_rts_peel(weighted_sources=weighted_sources,num_sources=int(args.num_sources),
                       output_name=output_name)
    #If you are doing aocalibrate
    else:
        if args.outside:
            sys.exit('outside option not yet supported for aocalibrate - NO SOURCELIST WRITTEN')
        else:
            output_name = "%s_%s_aocal%s.txt" %(args.srclist.split('/')[-1].split('.')[0],obsID,args.num_sources)
        print("++++++++++++++++++++++++++++++++++++++\nCreating weighted srclist for AO calibrate - not mega-patching the sources")

        write_aocal(weighted_sources=weighted_sources,output_name=output_name)

else:
    if args.aocalibrate:
        print('Mega patching not supported for aocalibrate - NO SOURCELIST WRITTEN')
        sys.exit()

    elif args.order=='flux':
        ordered_sources = weighted_sources

    elif args.order=='distance':
        ordered_offsets = [source.offset for source in weighted_sources]
        ordered_sources = [source for offset,source in sorted(zip(ordered_offsets,weighted_sources),key=lambda pair: pair[0])]

    elif 'name' in args.order:
        name = args.order.split("=")[1]
        top_source_ind = [source.name for source in weighted_sources].index(name)
        top_source = weighted_sources[top_source_ind]
        ordered_sources = [top_source]
        for i in arange(len(weighted_sources)):
            if i!=top_source_ind:
                ordered_sources.append(weighted_sources[i])
        print('++++++++++++++++++++++++++++++++++++++\nBase Source forced as %s with \nconvolved flux %.1fJy at a distance %.2fdeg\n---------------------------------' %(top_source.name,top_source.weighted_flux,top_source.offset))

    elif 'experimental' in args.order:
        if len(args.order.split('='))==1:
            flux_cut,dist_cut = 10.0,1.0
        else:
            flux_cut,dist_cut = args.order.split('=')[1].split(',')

        close_fluxs = []
        close_dists = []
        dist_cut = float(dist_cut)
        dist_cut_lower = 0.0

        ##Try to find all sources within distance cutoff above flux threshold - if none exist, extend search
        ##radii by 0.5 deg
        while len(close_fluxs)==0:
            for index in arange(len(weighted_sources)):
                src = weighted_sources[index]
                flux = src.weighted_flux
                offset = src.offset
                gauss_len = len(src.gaussians)

                if flux>float(flux_cut) and dist_cut_lower<offset<dist_cut:
                    close_fluxs.append(flux)
                    close_dists.append(offset)
            dist_cut+=0.5
            dist_cut_lower = dist_cut - 0.5
            if dist_cut_lower > cutoff:
                print("++++++++++++++++++++++++++++++++++++++\nNo source above %.2fJy within initial cutoff distance\nNO SOURCE LIST GENERATED\n++++++++++++++++++++++++++++++++++++++" %float(flux_cut))
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

        print("++++++++++++++++++++++++++++++++++++++\nBase source %s convolved flux is %.3fJy at a distance \nof %.3fdeg from point centre\n---------------------------------" %(brightest_close_name,brightest_close_flux,brightest_close_offset))

        ##Put this source at the top of the ordered list, and then append all other sources after
        ##NOTE - this means that apart from the top source, all other sources are flux ordered.
        ordered_sources = [weighted_sources[brightest_ind]]
        for source in weighted_sources:
            if source.weighted_flux!=brightest_close_flux and source.offset!=brightest_close_offset:
                ordered_sources.append(source)

    ##If added a zero Jy point source to give a beam centre direction to the patch,
    ##create the zero Jy source here and shove at the front of the ordered sources list
    if args.zeroJy_source:

        src = rts_source()
        src.name = 'pointing'
        src.ras.append(ra_point / 15.0)
        src.decs.append(dec_point)
        src.freqs.append([160e+6])
        src.fluxs.append([0.0])

        ordered_sources = [src] + ordered_sources

    ##Make a new single patch source based on the user specified number of components
    if args.outside:
        output_name = "%s_%s_outside-cutoff_patch%s.txt" %(args.srclist.split('/')[-1].split('.')[0],obsID,args.num_sources)
    else:
        output_name = "%s_%s_patch%s.txt" %(args.srclist.split('/')[-1].split('.')[0],obsID,args.num_sources)

    write_rts_patch(ordered_sources=ordered_sources,num_sources=int(args.num_sources),
                    output_name=output_name,zeroJy_source=args.zeroJy_source,
                    extra_calibrators=extra_calibrators,obsID=obsID)

print("Created %s\n++++++++++++++++++++++++++++++++++++++" %output_name)
##Finito!!
