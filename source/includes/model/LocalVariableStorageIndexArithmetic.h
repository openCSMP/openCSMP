
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
