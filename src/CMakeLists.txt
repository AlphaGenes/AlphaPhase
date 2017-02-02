
# Alphaphase library
SET(ALPHAPHASE_src
            CoreDefinition.f90
            CoreUtils.f90
            CoreSubsetDefinition.f90
            HaplotypeLibraryDefinition.f90
            HaplotypeLibraryPhasing.f90
            InputOutput.f90
            LongRangePhasing.f90
            MemberManagerDefinition.f90
            ParametersDefinition.f90
	    ProgramParametersDefinition.f90
            PedigreeDefinition.f90
            Random.f90
            SurrogateDefinition.f90
            TestResultDefinition.f90        
	    AlphaPhaseFunctions.f90
	    AlphaPhaseResultsDefinition.f90
            AlphaPhase.f90 
            
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
