###################################################################################
# Script name   : filterStudyAnimalSpeciesStrain.R
# Date Created  : 16-Jan-2020
# Programmer    : Bo Larsen
# --------------------------------------------------------------------------------
# Change log: 
# Programmer/date     Description
# -----------------   ------------------------------------------------------------
# <init/dd-Mon-yyyy>  <description>
#
# -------------------------------------------------------------------------------
# Purpose       : Extract studies and animals which fulfills a specified species
#                 value and optional strain value from a pooled SEND data store.
#
# Description   : Function FilterAnimalsSpeciesStrain:
# 
#                 Returns a data table with a set of animals extracted from the table of  
#                 animals given as input where the animals fits the species value(s) given 
#                 as input (a list of multiple species values may be given as input) 
#                 and also fits the strain value(s) if given as input (optional).
#                 Animals from the input set, which are included in studies where one of  
#                 of the input species values - and optionally also strain values - is 
#                 registered in TS (and TS have only this species or species/strain 
#                 value included) are included in the output set.
#                 Animals from the input set, belonging to studies 
#                   - with no species/strain value registered  or  
#                   - with more than one species/strain value registered
#                 in TS are included in the output set if they exists in DM with matching 
#                 DM.SPECIES or DM.SPECIES/STRAIN values or are included in TX in a trial set 
#                 with matching TX.SPECIES or TX.SPECIES/STRAIN (i.e. TXVAL where TXPARMCD is 
#                 SPECIES' or 'STRAIN')
#                 The comparisons of the species and strain values are done case 
#                 insensitive.
#                 
#                 If the input parameter inclUncertain flag is enabled, uncertain animals
#                 are included in the output set.
#                 These uncertain situations are identified and reported for SPECIES and STRAIN respectively 
#                 (in column UNCERTAIN_MSG):
#                  - TS parameter SPECIES/STRAIN is missing or invalid (not CT value - CDISC code list SPECIES/STRAIN) 
#                    and TX parameter SPECIES/STRAIN is missing or invalid (not CT value) and DM.SPECIES/STRAIN is 
#                    missing or invalid (not CT value)
#                  - Different values of SPECIES/STRAIN across TS, TX and DM for studies where no or only one 
#                    TS parameter SPECIES/STRAIN is registered
#                  - Multiple TS parameter SPECIES/STRAIN values are registered for study and TX parameter 
#                    SPECIES/STRAIN and/or DM.SPECIES/STRAIN do not match any of the TS values.
#                  - Multiple TS parameter SPECIES/STRAIN values are registered for study and TX parameter 
#                    SPECIES/STRAIN and DM.SPECIES/STRAIN are unequal.
#                  Non-empty UNCERTAIN_MSG values are merged with non-empty UNCERTAIN_MSG values 
#                  which may exist in the input set of animals (animalList).
#
# Input         : - The TS, TX and DM domains - are imported from the pooled SEND data store if 
#                   the don't exist in workspace.
#                 - A data tables specified in the input parameters animalList:
#                   It contains the list of animals to filter for specified species and strain value(s)
#                   - must contain these character variables:
#                       STUDYID
#                       USUBJID
#                     other variables may be included 
#                 - CDISC CT code lists SPECIES and STRAIN imported from a CDISC CT file.
#                   
# Output        : A data table with the character columns:
#                   STUDYID
#                   USUBJID
#                   SPECIES
#                   STRAIN (if strainFilter has been specified)
#                   UNCERTAIN_MSG - if input parameter inclUncertain flag is enabled
#                 plus any additional columns which may be included in the input data animalList
#
# Parameters    : animalList:     Mandatory, data table (see Input).
#                 speciesFilter:  Mandatory, character.
#                                   The species value(s) to use as criterion for filtering of the input data table.
#                                   It can be a single string, a vector or a list of multiple strings.
#                 strainFilter:   Optional, character.
#                                   The species value(s) to use as criterion for filtering of the input data table.
#                                   It can be a single string, a vector or a list of multiple strings.
#                                   Only allowed when a single value is specified for speciesFilter.
#                 inclUncertain:  Optional, Include uncertain rows or not
#                   
###################################################################################

library(data.table)

GetSpeciesStrain<-function(animalList=FALSE) {
  


  
  # import TS, TX and DM if they not already exists
  if (!exists("TS")) {
    importSENDDomains(c("TS"))
  }
  if (!exists("TX")) {
    importSENDDomains(c("TX"), animalStudies)
  }
  if (!exists("DM")) {
    importSENDDomains(c("DM"), animalStudies)
  }
  
  ##################################################################################################################
  
  
  if (!animalList) {
    animalStudies <- unique(TS[,.(STUDYID)])
    animalList <- unique(DM[,.(STUDYID, USUBJID)])
  }

  # Extract all TS rows for parameter SPECIES, rename TSVAL to SPECIES_TS 
  # - remove duplicates
  # - limit to the set of studies for input set of animals
  tsSPECIESall<-merge(unique(TS[TSPARMCD == 'SPECIES', .(STUDYID,SPECIES_TS=toupper(trimws(TSVAL)))]), animalStudies, by='STUDYID')
  
  
  
  # Add studies with no TS parameter SPECIES from the set of studies for input set of animals
  tsSPECIESall<-
    rbindlist(list(tsSPECIESall,
                   fsetdiff(animalStudies, 
                            tsSPECIESall[,.(STUDYID)])[!is.na(STUDYID),.(STUDYID, SPECIES_TS=as.character(NA))]),
              use.names=TRUE, fill=TRUE)
  
  # Add variables with 
  #  - count of number of distinct SPECIES per study
  #  - concatenation of all species per study (for studies with one species this is equal to SPECIES_TS)
  tsSPECIESall[, `:=` (ALL_SPECIES = unique(SPECIES_TS), NUM_SPECIES_TS = .N), by = STUDYID]
  tsSPECIESall[,`:=`(ALL_SPECIES_TS = c(.SD)), by = STUDYID, .SDcols='SPECIES_TS']
  
  # JOin the list of studies/species with DM to get all animal level SPECIES 
  #  - the list will contain STUDYID/USUBJID duplicates for studies with multiple SPECIES registered in TS               
  animalSPECIESall<-               
    merge(merge(tsSPECIESall,
                animalList[,.(STUDYID, USUBJID)],
                by='STUDYID', allow.cartesian = TRUE),
          DM[,.(STUDYID, USUBJID, SETCD, SPECIES_DM=ifelse(SPECIES=="",as.character(NA),toupper(trimws(SPECIES))))],
          by=c('STUDYID','USUBJID'), allow.cartesian = TRUE)
  
  # Join the list of studies/animals/species with TX to get all set level SPECIES
  #  - add variable SPECIES with the first non-empty species value from DM, TX or TS
  #  - add variable with count of unique USUBJID per study (there is expected to be one usubjid per studyid per TSPARMCD 'SPECIES' )
  animalSPECIESall<-
    merge(animalSPECIESall, 
          unique(TX[TXPARMCD=='SPECIES',.(STUDYID,SETCD,SPECIES_TX=toupper(trimws(TXVAL)))]), 
          by=c('STUDYID','SETCD'), all.x=TRUE )[,`:=` (SPECIES=fcoalesce(SPECIES_DM,SPECIES_TX,SPECIES_TS))]
  animalSPECIESall[, `:=` (NUM_ANIMALS = .N), by = .(STUDYID, USUBJID)]
  
  # Extract unique set of animals and related species - exclude uncertain animals
  animalSPECIESallUniq<-unique(animalSPECIESall[,.(STUDYID, USUBJID, SPECIES)])
  
  

  ##########################################################################################################################
  

  # List of studyid values included in the list of animals extracted at SPECIES level
  animalSPECIESStudies<-unique(animalSPECIESallUniq[,.(STUDYID)])
  
  # Extract all TS rows for parameter STRAIN, rename TSVAL to STRAIN_TS 
  # - remove duplicates
  # - limit to the set of studies for input set of animals
  tsSTRAINall<-merge(unique(TS[TSPARMCD == 'STRAIN', .(STUDYID,STRAIN_TS=toupper(trimws(TSVAL)))]), animalSPECIESStudies, by='STUDYID')
  # Add studies with no TS parameter STRAIN from the set of studies for input set of animals
  tsSTRAINall<-
    rbindlist(list(tsSTRAINall,
                   fsetdiff(animalSPECIESStudies, 
                            tsSTRAINall[,.(STUDYID)])[!is.na(STUDYID),.(STUDYID, STRAIN_TS=as.character(NA))]),
              use.names=TRUE, fill=TRUE)
  
  # Add variables with 
  #  - count of number of distinct STRAIN per study
  #  - concatenation of all strain per study (for studies with one strain this is equal to STRAIN_TS)
  tsSTRAINall[, `:=` (ALL_STRAIN = unique(STRAIN_TS), NUM_STRAIN_TS = .N), by = STUDYID]
  tsSTRAINall[,`:=`(ALL_STRAIN_TS = c(.SD)), by = STUDYID, .SDcols='STRAIN_TS']
  
  # JOin the set of studies/strain with DM to get all animal level STRAIN  limited to the found set of animals at species level
  #  - the set will contain STUDYID/USUBJID duplicates for studies with multiple STRAIN registered in TS               
  animalSTRAINall<-               
    merge(merge(tsSTRAINall,
                animalList[,.(STUDYID, USUBJID)],
                by='STUDYID'),
          merge(animalSPECIESallUniq[,.(STUDYID,USUBJID)], 
                DM[,.(STUDYID, USUBJID, SETCD, STRAIN_DM=ifelse(STRAIN=="",as.character(NA),toupper(trimws(STRAIN))))],
                by=c('STUDYID','USUBJID')),
          by=c('STUDYID','USUBJID'))
  
  # Join the list of studies/animals/strain with TX to get all set level STRAIN
  #  - add variable STRAIN with the first non-empty strain value from DM, TX or TS
  #  - add variable with count of unique USUBJID per study (there is expected to be one usubjid per studyid per TSPARMCD 'STRAIN' )
  animalSTRAINall<-
    merge(animalSTRAINall, 
          unique(TX[TXPARMCD=='STRAIN',.(STUDYID,SETCD,STRAIN_TX=toupper(trimws(TXVAL)))]), 
          by=c('STUDYID','SETCD'), all.x=TRUE )[,`:=` (STRAIN=fcoalesce(STRAIN_DM,STRAIN_TX,STRAIN_TS))]
  animalSTRAINall[, `:=` (NUM_ANIMALS = .N), by = .(STUDYID, USUBJID)]

  # Extract unique set of animals and related strain - exclude uncertain animals
  animalSTRAINallUniq<-unique(animalSTRAINall[,.(STUDYID, USUBJID, STRAIN)])

  

  # Merge set of animals with the set of animals extracted at SPECIES level to get the SPECIES values
  foundAnimals<-merge(animalSPECIESallUniq,
                      animalSTRAINallUniq,
                      by=c('STUDYID','USUBJID'))
  
  return(foundAnimals)
  
}
