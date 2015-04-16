##
#  CMake script for the step-1 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "step-1")

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cpp
  InitialValuesU.cpp
  InverseDiffusivity.cpp
  DirichletBoundaryValues.cpp
  DirichletBoundaryValuesDot.cpp
  DirichletBoundaryValuesPoisson.cpp
  SourceBodyValues.cpp  
  NeumannBoundaryValues.cpp
  EllipticProblem.cpp
  #  SchurComplement.cpp
  PSigmaPStar.cpp
  Q_UStar_SStar_inverseMassQ_S_U_Q.cpp
  ApproximateSchurComplement.cpp
  SchurComplementPoisson.cpp
  ApproximateSchurComplementPoisson.cpp
  NeumanBoundaryValuesPoisson.cpp
  # You can specify additional files here!
  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../
  /home/kvanvels/share /h1/kvanvels/workspace/deal.II
/home/kvanvels/deal.II $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
