
# Alphaphase library
SET(ALPHAPHASE_src
            CoreModule.f90
            CoreUtils.f90
            CoreSubsetModule.f90
            HaplotypeLibraryModule.f90
            HaplotypeLibraryPhasing.f90
            InputOutput.f90
            LongRangePhasing.f90
            MemberManagerModule.f90
	    OutputParametersModule.f90
            AlphaPhaseParametersModule.f90
	    ProgramParametersModule.f90
            SurrogateModule.f90
            TestResultModule.f90        
	    AlphaPhaseFunctions.f90
	    AlphaPhaseResultsModule.f90
            AlphaPhase.f90 
            
)

if (ALPHAPHASEEXE)
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${ALPHAPHASEEXE} ${ALPHAPHASE_src})
    TARGET_LINK_LIBRARIES(${ALPHAPHASEEXE} ${AHLIB})
    TARGET_LINK_LIBRARIES(${ALPHAPHASEEXE} ${RESULT})
endif()
if(APLIB)

    SET(ALPHAPHASELIB_src  ${ALPHAPHASE_src}) 
    ADD_LIBRARY(${APLIB} STATIC ${ALPHAPHASELIB_src})
    TARGET_LINK_LIBRARIES(${APLIB} ${AHLIB})
endif()
# SET(HDF5TEST_src ${SRCALPHAIMPUTE}/hdf5test.f90 ${SRCALPHAIMPUTE}/hdf5module.f90)
# ADD_EXECUTABLE(hdf5test ${HDF5TEST_src})

# add external links here
# 

INSTALL(TARGETS ${ALPHAPHASEEXE} RUNTIME DESTINATION bin)
