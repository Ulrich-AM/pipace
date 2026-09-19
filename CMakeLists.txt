cmake_minimum_required(VERSION 3.20)
project(PIPACE LANGUAGES CXX)

add_executable(pipace WIN32
  main.cpp
  fluid/FluidEngine.cpp
  fluid/WorkerPool.cpp
  rigid/RigidBodyEngine.cpp
  ui/UiShell.cpp
  ui/UiLanguage.cpp
  ui/UiAssets.cpp
  ui/ToolWindow.cpp
  gas/GasEngine.cpp
  thermal/ThermalEngine.cpp
  render/WorldRenderer.cpp
  substance/SubstanceRegistry.cpp
  substance/PhaseTransfer.cpp
  substance/LiquidMixtureProperties.cpp
  chemistry/ReactionRegistry.cpp
  chemistry/ReactionMatterAccess.cpp
  chemistry/ReactionEngine.cpp
  world/WorldQuery.cpp
  world/WaterPhaseChange.cpp
  sim/SimulationQuality.cpp
)
target_include_directories(pipace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(pipace PRIVATE cxx_std_17)
target_compile_definitions(pipace PRIVATE UNICODE _UNICODE NOMINMAX WIN32_LEAN_AND_MEAN)
target_link_libraries(pipace PRIVATE user32 gdi32 gdiplus)

if(MSVC)
  target_compile_options(pipace PRIVATE /W4 /permissive-)
else()
  target_compile_options(pipace PRIVATE -Wall -Wextra -Wpedantic)
  if(MINGW)
    target_link_options(pipace PRIVATE -static-libgcc -static-libstdc++ -municode -mwindows)
  endif()
endif()

add_custom_command(TARGET pipace POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/assets $<TARGET_FILE_DIR:pipace>/assets)
