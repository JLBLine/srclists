############################################################
## A repository for RTS srclists --------------------------#
############################################################
##02/02/2015 - J Line
 - Added the most up to date and complete PUMA srclist,
   srclist_puma-v2_complete.txt. This srclist covers the
   full southern sky, using either mwacs or mrc as a base.
   It is complete with point source extended models between
   22h < RA < 6h, -57deg < Dec < +3deg (2 hours about EoR0
   and EoR1). Elsewhere, only the sources accepted by PUMA
   are currently included.
   
 - In ./puma_v2, the separate srclists used to create
   srclist_puma-v2_complete.txt are stored and explained
   
 - In ./extended, information and plots of all the extended
   sources in srclist_puma-v2_complete.txt are included
   
 - srclist_by_beam.py will take a srclist, and return a
   given number of sources in a megapatch srclist. These 
   will be brightest sources once convolved with the beam, 
   as defined by the pointing given by a designated metafits
   file. The sources are also confined to a user specified
   distance of the pointing centre. To use the plotting option 
   (-p), wscaxes module is needed.
