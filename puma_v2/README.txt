############################################################
## A repository for RTS srclists --------------------------#
############################################################
##02/02/2015 - J Line
  - srclist_puma-v2_mwacs.txt contains all of the accepted
    matches by PUMA (v2) when matching mwacs to vlssr, mrc,
    sumss, and nvss.
    
  - srclist_puma-v2_mrcred.txt contains all of the accepted
    matches when matching mrc to vlssr, sumss, and nvss. The
    mrc base was a reduced mrc, essentially where ever mwacs
    didn't have coverage.
    
  - srclist_hack_v2.txt contains all the extended point
    source models, created by hand using the 'eyeball'
    output of PUMA, between 22h < RA < 6h, 
    -57deg < Dec < +3deg.
    
  - srclist_puma-v2_complete.txt is all of the above
    combined
    
In each catalogue, if a spectrum had a chi-squared of less
than 2, fitted values were used to keep the data as compact
as possible.
