# A repository for Aussie MWA EoR srclists 

Current suggested sky model to use is:

```
srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07.fits
```

## Hyperdrive / WODEN sky models

### 07/11/2023
`srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07.fits` \
`srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07.yaml`
 - LoBES (Lynch et al 2021) in EoR0 field
 - Procopio et al 2017 with a few added shapelet models in EoR1 field (including Fornax A model from Line et 2020)
 - Cook et all 2022 for Cen A and Galactic Plane SNR
 - Everything is either a power-law or curved power-law SED

## RTS sky models (deprecated)

### 28/06/2019
`srclist_pumav3_EoR0aegean_EoR1pietro+ForA_phase1+2.txt`
 - Added a much improved Fornax A model using both phase 1 and phase 2 extended data

`srclist_pumav3_EoR0aegean_EoR1pietro+ForA_phase1+2_TGSSgalactic.txt`
 - New ForA model plus attempt at filling galactic plane with TGSS sources

### 12/07/2018
`srclist_pumav3_EoR0aegean_EoR1pietro+ForA_TGSSgalactic.txt`

 - Attempt to fill in the galactic plane with TGSS sources

### 16/02/2017
`srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt`

 The main features of this srclist are:
   + Based off of the GLEAM catalogue
   + Matched to TGSS, VLSSr, SUMMS, and NVSS
   + Has PUMA based extended models up to 20 deg away from EoR0
   + Has Pietro's extended models up to 10 deg away from EoR1
   + Outside of those areas, it is the PUMA matched GLEAM catalogue
     if a GLEAM source failed the PUMA tests, it has been included
     but just as a GLEAM source - it has not been eyeballed
   + Also has Jenny's shapelet model for Fornax A

### 02/02/2015 - J Line
`srclist_puma-v2_complete.txt`

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
