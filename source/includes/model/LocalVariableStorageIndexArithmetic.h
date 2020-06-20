#ifndef CSMP_LOCAL_VARIABLE_STORAGE_INDEX_ARITHMETIC_H
#define CSMP_LOCAL_VARIABLE_STORAGE_INDEX_ARITHMETIC_H

// To be absolutely safe while preprocessor directives are always global, some checkes are done here.
#ifdef DATA_OFFSET_TO_IP
LocalVariableStorageIndexArithmetic: macro DATA_OFFSET_TO_IP for local variable storage has already been defined elsewhere!
#endif
#ifdef DATA_OFFSET_IP
LocalVariableStorageIndexArithmetic: macro DATA_OFFSET_IP for local variable storage has already been defined elsewhere!
#endif
#ifdef FLAG_OFFSET_TO_IP
LocalVariableStorageIndexArithmetic: macro FLAG_OFFSET_TO_IP for local variable storage has already been defined elsewhere!
#endif
#ifdef FLAG_OFFSET_IP
LocalVariableStorageIndexArithmetic: macro FLAG_OFFSET_IP for local variable storage has already been defined elsewhere!
#endif
#ifdef DATA_OFFSET_TO_FVIP
LocalVariableStorageIndexArithmetic: macro DATA_OFFSET_TO_FVIP for local variable storage has already been defined elsewhere!
#endif
#ifdef DATA_OFFSET_FVIP
LocalVariableStorageIndexArithmetic: macro DATA_OFFSET_FVIP for local variable storage has already been defined elsewhere!
#endif
#ifdef FLAG_OFFSET_TO_FVIP
LocalVariableStorageIndexArithmetic: macro FLAG_OFFSET_TO_FVIP for local variable storage has already been defined elsewhere!
#endif
#ifdef FLAG_OFFSET_FVIP
LocalVariableStorageIndexArithmetic: macro FLAG_OFFSET_FVIP for local variable storage has already been defined elsewhere!
#endif

// CRT calls for elements/faces/interfaces
#define FE_SES static_cast<const STOREE*>(this)->Sectors()                                  ///< sectors of element storee
#define FE_FAS static_cast<const STOREE*>(this)->Facets()                                   ///< facets of element storee
#define FE_FVIPS_PER_SECTOR static_cast<const STOREE*>(this)->IntegrationPointsPerSector()  ///< integration points per sector
#define FE_FVIPS_PER_FACET static_cast<const STOREE*>(this)->IntegrationPointsPerFacet()    ///< integration points per facet
#define IPS_SI static_cast<const STOREE*>(this)->IntegrationPoints()                        ///< simplex integration points
#define IPS_SE FE_SES*FE_FVIPS_PER_SECTOR                                                   ///< sector integration points
#define IPS_FA FE_FAS*FE_FVIPS_PER_FACET                                                    ///< facet integration points

// element/face/interface integration point offsets

// offset to the first element integration point value stored
#define DATA_OFFSET_TO_IP idx.localVariables.totalDataDepth 

#define DATA_OFFSET_IP DATA_OFFSET_TO_IP + ip * idx.integrationPointVariables.ipvSimplex.totalDataDepth + idx.dataOffset

#define FLAG_OFFSET_TO_IP idx.localVariables.totalFlagDepth 

#define FLAG_OFFSET_IP FLAG_OFFSET_TO_IP + ip * idx.integrationPointVariables.ipvSimplex.totalFlagDepth + idx.flagOffset


// element/face/interface finite volume (sector/facet) integration point offsets

#define DATA_OFFSET_TO_FVIP idx.localVariables.totalDataDepth + IPS_SI * idx.offsetFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalDataDepth\
                            + IPS_SE * idx.offsetFactorSector * idx.integrationPointVariables.ipvSector.totalDataDepth

#define DATA_OFFSET_FVIP DATA_OFFSET_TO_FVIP + sector_or_facet * idx.ipFactorSector *idx.integrationPointVariables.ipvSector.totalDataDepth\
                         + sector_or_facet * idx.ipFactorFacet *idx.integrationPointVariables.ipvFacet.totalDataDepth\
                         + ip * idx.dataDepth + idx.dataOffset

#define FLAG_OFFSET_TO_FVIP idx.localVariables.totalFlagDepth + IPS_SI * idx.offsetFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalFlagDepth\
                            + IPS_SE * idx.offsetFactorSector * idx.integrationPointVariables.ipvSector.totalFlagDepth

#define FLAG_OFFSET_FVIP FLAG_OFFSET_TO_FVIP +  sector_or_facet * idx.ipFactorSector *idx.integrationPointVariables.ipvSector.totalFlagDepth\
                         + sector_or_facet * idx.ipFactorFacet *idx.integrationPointVariables.ipvFacet.totalFlagDepth\
                         + ip * idx.flagDepth + idx.flagOffset

#endif /* CSMP_LOCAL_VARIABLE_STORAGE_INDEX_ARITHMETIC_H - macro definitions */
