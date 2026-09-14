cmake_minimum_required(VERSION 3.20)
project(WaterTest LANGUAGES CXX)

add_executable(water_test WIN32
  main.cpp
  fluid/FluidEngine.cpp
  fluid/WorkerPool.cpp
  rigid/RigidBodyEngine.cpp
  ui/UiShell.cpp
  ui/UiLanguage.cpp
  ui/UiAssets.cpp
  gas/GasEngine.cpp
  thermal/ThermalEngine.cpp
)
target_include_directories(water_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(water_test PRIVATE cxx_std_17)
target_compile_definitions(water_test PRIVATE UNICODE _UNICODE NOMINMAX WIN32_LEAN_AND_MEAN)
target_link_libraries(water_test PRIVATE user32 gdi32 gdiplus)

if(MSVC)
  target_compile_options(water_test PRIVATE /W4 /permissive-)
else()
  target_compile_options(water_test PRIVATE -Wall -Wextra -Wpedantic)
  if(MINGW)
    target_link_options(water_test PRIVATE -static-libgcc -static-libstdc++ -municode -mwindows)
  endif()
endif()

add_custom_command(TARGET water_test POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/assets $<TARGET_FILE_DIR:water_test>/assets)
