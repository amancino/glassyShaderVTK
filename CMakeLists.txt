
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(GlassyShader)

find_package(VTK COMPONENTS 
  vtkCommonColor
  vtkCommonCore
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersSources
  vtkIOGeometry
  vtkIOLegacy
  vtkIOPLY
  vtkIOImage
  vtkIOXML
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping BozoShader: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(GlassyShader MACOSX_BUNDLE GlassyShader.cxx )
  target_link_libraries(GlassyShader PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(GlassyShader MACOSX_BUNDLE GlassyShader.cxx )
  target_link_libraries(GlassyShader PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS GlassyShader
    MODULES ${VTK_LIBRARIES}
    )
endif () 
