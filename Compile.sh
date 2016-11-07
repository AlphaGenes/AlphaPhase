rm alphaphase
rm *.mod
ifort -fpp -static-intel -O3 -m64 -o alphaphase -DCOMMIT=" 26a29e5" src/Constants.f90 src/ParametersDefinition.f90 src/Random.f90 src/Sorting.f90 src/PedigreeDefinition.f90 src/NRMcode.f90 src/CoreUtils.f90 src/CoreDefinition.f90 src/CoreSubsetDefinition.f90 src/Clustering.f90 src/HaplotypeLibraryDefinition.f90 src/SurrogateDefinition.f90 src/LongRangePhasing.f90 src/MemberManagerDefinition.f90 src/TestResultDefinition.f90 src/InputOutput.f90 src/HaplotypeLibraryPhasing.f90 src/AlphaPhase.f90
