#---------------------------------------------------------------------------
set(Physics_files 
   InSANEFunctionManager 
   InSANEDiffXSec 
   InSANECompositeDiffXSec 
   InSANEDiffXSecKinematicKey 
   InSANEGridDiffXSec 
   InSANEElectroProductionXSec 
   InSANEepElasticDiffXSec 
   InSANEXSections 
   InSANEInclusiveDiffXSec 
   InSANEInclusiveBornDISXSec 
   InSANEExclusiveDiffXSec 
   MAIDKinematicKey   
   MAIDInclusiveDiffXSec   
   MAIDNucleusInclusiveDiffXSec   
   MAIDPolarizedKinematicKey   
   MAIDPolarizedTargetDiffXSec   
   MAIDExclusivePionDiffXSec   
   MAIDExclusivePionDiffXSec2   
   MAIDInclusiveElectronDiffXSec   
   MAIDInclusivePionDiffXSec   
   InSANEVirtualComptonAsymmetries   
   MAIDVirtualComptonAsymmetries   
   F1F209eInclusiveDiffXSec 
   F1F209StructureFunctions 
   F1F209QuasiElasticFormFactors 
   CTEQ6eInclusiveDiffXSec 
   NMC95StructureFunctions 
   LHAPDFStructureFunctions 
   PolarizedDISXSec 
   QFSInclusiveDiffXSec 
   WiserXSection 
   EPCVXSection 
   OARPionDiffXSec 
   BETAG4StructureFunctions 
   InSANEStructureFunctions 
   InSANECompositeStructureFunctions 
   InSANECompositePolarizedStructureFunctions 
   InSANEPolarizedStructureFunctions 
   InSANEStructureFunctionsFromPDFs 
   InSANEPolSFsFromComptonAsymmetries 
   InSANEPolarizedStructureFunctionsFromPDFs 
   InSANEFormFactors 
   MSWFormFactors 
   AmrounFormFactors 
   BilenkayaFormFactors 
   InSANEBeamSpinAsymmetry 
   InSANEPhaseSpaceSampler 
   InSANEPhaseSpaceVariable 
   InSANEPhaseSpace 
   InSANEPartonDistributionFunctions 
   InSANEPolarizedPartonDistributionFunctions 
   StatisticalQuarkFits     
   StatisticalUnpolarizedPDFs 
   StatisticalPolarizedPDFs 
   BBSQuarkHelicityDistributions     
   BBSUnpolarizedPDFs 
   BBSPolarizedPDFs
   AvakianQuarkHelicityDistributions     
   AvakianUnpolarizedPDFs 
   AvakianPolarizedPDFs 
   DSSVPolarizedPDFs 
   AAC08PolarizedPDFs 
   BBPolarizedPDFs 
   JAMPolarizedPDFs 
   MHKPolarizedPDFs 
   GSPolarizedPDFs 
   DNS2005PolarizedPDFs 
   LSS2006PolarizedPDFs 
   LSS2010PolarizedPDFs 
   LHAPDFUnpolarizedPDFs 
   LHAPDFPolarizedPDFs 
   CTEQ6UnpolarizedPDFs 
   CJ12UnpolarizedPDFs 
   LCWFPartonDistributionFunctions 
   InSANEPolarizedCrossSectionDifference 
   InSANEAsymmetryBase 
   InSANEAsymmetriesFromStructureFunctions 
   InSANECrossSectionDifference 
   InSANEPOLRAD    
   InSANEPOLRADKinematics    
   InSANEPOLRADInternalPolarizedDiffXSec 
   InSANEPOLRADQuasiElasticTailDiffXSec 
   InSANERADCOR 
   InSANERADCORKinematics 
   InSANERADCORRadiatedDiffXSec 
   InSANERADCORRadiatedUnpolarizedDiffXSec 
   InSANERADCORInternalUnpolarizedDiffXSec 
   InSANERADCOR2 
   InSANERadiativeTail 
   InSANEElasticRadiativeTail 
   InSANEInelasticRadiativeTail 
   InSANEInelasticRadiativeTail2 
   InSANEPOLRADUltraRelativistic 
   QuasiElasticInclusiveDiffXSec 
   QEIntegral  
   LSS98QuarkHelicityDistributions
   LSS98UnpolarizedPDFs
   LSS98PolarizedPDFs
   InSANEPartonHelicityDistributions
   InSANEPartonDistributionFunctionsFromPHDs
   ABKM09UnpolarizedPDFs
   CTEQ10UnpolarizedPDFs
   MSTW08UnpolarizedPDFs
   MRST2001UnpolarizedPDFs
   MRST2002UnpolarizedPDFs
   MRST2006UnpolarizedPDFs
   InSANERadiatorBase
   InSANERadiator
   )

set(Physics_SRCS)
set(Physics_HEADERS)
foreach(infileName ${Physics_files})
   SET(Physics_SRCS ${Physics_SRCS} "${PROJECT_SOURCE_DIR}/src/${infileName}.C")
   SET(Physics_HEADERS ${Physics_HEADERS} "${PROJECT_SOURCE_DIR}/include/${infileName}.h")
endforeach(infileName)

set(PhysicsFortran_files
   F1F209_FAST
   sane_pol 
   wiser_func 
   wiser_fit 
   qfs_func 
   qfs_sigs 
   qfs_targs 
   ppdf 
   Cteq6Pdf2010 
   CJ12pdf 
   aac08 
   LSS2006pdf_g1 
   LSS2010_pdfs 
   polnlo 
   dssv 
   mrst2001 
   mrst2002 
   mrst2006 
   partondf 
   nmc_org 
   r1998 
   epcvnew 
   epc_or_v3_funcs 
   epc_or_v3_prog 
   readgrid 
   Amroun 
   nqfs
   abkm09
   cteq10
   mstwpdf
   )

set(PhysicsFortran_SRCS)
foreach(infileName ${PhysicsFortran_files})
   SET(PhysicsFortran_SRCS ${PhysicsFortran_SRCS} "${PROJECT_SOURCE_DIR}/src/${infileName}.f")
endforeach(infileName)


# set everything needed for the root dictonary and create the
# dictionary
set(Physics_LINKDEF ${PROJECT_SOURCE_DIR}/include/InSANEPhysics_LinkDef.h )
set(Physics_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/InSANEPhysicsDict.cxx) 
ROOT_GENERATE_DICTIONARY("${Physics_HEADERS}" "${Physics_LINKDEF}" "${Physics_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

# add the dictionary to the list of source files
SET(Physics_SRCS ${PhysicsFortran_SRCS} ${Physics_SRCS} ${Physics_DICTIONARY}) 

# Set the library version in the main CMakeLists.txt
SET(Physics_MAJOR_VERSION 0)
SET(Physics_MINOR_VERSION 0)
SET(Physics_PATCH_VERSION 0)
SET(Physics_VERSION "${Physics_MAJOR_VERSION}.${Physics_MINOR_VERSION}.${Physics_PATCH_VERSION}")
SET(Physics_LIBRARY_PROPERTIES ${Physics_LIBRARY_PROPERTIES}
    VERSION "${Physics_VERSION}"
    SOVERSION "${Physics_MAJOR_VERSION}"
    SUFFIX ".so"
)

add_library(          InSANEPhysics SHARED ${Physics_SRCS})
target_link_libraries(InSANEPhysics ${ROOT_LIBRARIES})
set_target_properties(InSANEPhysics PROPERTIES ${Physics_LIBRARY_PROPERTIES})

#install(TARGETS Physics DESTINATION ${CMAKE_BINARY_DIR}/lib)

