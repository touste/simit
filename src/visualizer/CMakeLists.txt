set(VIZ_PROJECT_NAME ${PROJECT_NAME}viz)

find_package(OPENGL)
find_package(GLUT) 

if (NOT OPENGL_FOUND OR NOT GLUT_FOUND)
  message("Could not find OpenGL or GLUT so simviz will not be built")
else()
  set(SIMIT_VISUALIZER_SOURCE_DIR ${SIMIT_SOURCE_DIR}/visualizer)
  file(GLOB SOURCES *.cpp ${SIMIT_VISUALIZER_SOURCE_DIR}/*.cpp)
  file(GLOB HEADERS *.h ${SIMIT_VISUALIZER_SOURCE_DIR}/*.h)
  add_definitions(-DGL_GLEXT_PROTOTYPES)

  add_library(${VIZ_PROJECT_NAME} ${SIMIT_LIBRARY_TYPE} ${HEADERS} ${SOURCES})

  target_link_libraries(${VIZ_PROJECT_NAME} ${PROJECT_NAME})

  # GL
  include_directories(${OpenGL_INCLUDE_DIRS})
  link_directories(${OpenGL_LIBRARY_DIRS})
  add_definitions(${OpenGL_DEFINITIONS})

  include_directories(${GLUT_INCLUDE_DIRS})
  link_directories(${GLUT_LIBRARY_DIRS})
  add_definitions(${GLUT_DEFINITIONS})

  target_link_libraries(${VIZ_PROJECT_NAME} GL)
  target_link_libraries(${VIZ_PROJECT_NAME} GLU)
  target_link_libraries(${VIZ_PROJECT_NAME} glut)

  include_directories(${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})
endif()