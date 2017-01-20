
# Alphaphase library
SET(ALPHAPHASE_src
            ${SRCALPHAPHASE}/CoreDefinition.f90
            ${SRCALPHAPHASE}/CoreUtils.f90
            ${SRCALPHAPHASE}/CoreSubsetDefinition.f90
            ${SRCALPHAPHASE}/HaplotypeLibraryDefinition.f90
            ${SRCALPHAPHASE}/HaplotypeLibraryPhasing.f90
            ${SRCALPHAPHASE}/InputOutput.f90
            ${SRCALPHAPHASE}/LongRangePhasing.f90
            ${SRCALPHAPHASE}/MemberManagerDefinition.f90
            ${SRCALPHAPHASE}/ParametersDefinition.f90
	    ${SRCALPHAPHASE}/ProgramParametersDefinition.f90
            ${SRCALPHAPHASE}/PedigreeDefinition.f90
            ${SRCALPHAPHASE}/Random.f90
            ${SRCALPHAPHASE}/SurrogateDefinition.f90
            ${SRCALPHAPHASE}/TestResultDefinition.f90        
	    ${SRCALPHAPHASE}/AlphaPhaseFunctions.f90
	    ${SRCALPHAPHASE}/AlphaPhaseResultsDefinition.f90
            ${SRCALPHAPHASE}/AlphaPhase.f90 
            
)

if (ALPHAPHASEEXE)
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${ALPHAPHASEEXE} ${ALPHAPHASE_src})
    TARGET_LINK_LIBRARIES(${ALPHAPHASEEXE} ${AHLIB})
    TARGET_LINK_LIBRARIES(${ALPHAPHASEEXE} ${RESULT})
endif()
if(APLIB)
    ADD_LIBRARY(${APLIB} STATIC ${ALPHAPHASE_src})
endif()
# SET(HDF5TEST_src ${SRCALPHAIMPUTE}/hdf5test.f90 ${SRCALPHAIMPUTE}/hdf5module.f90)
# ADD_EXECUTABLE(hdf5test ${HDF5TEST_src})

# add external links here
# 

INSTALL(TARGETS ${ALPHAPHASEEXE} RUNTIME DESTINATION bin)
