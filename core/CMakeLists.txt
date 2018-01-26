FIND_PACKAGE(Ceres REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options filesystem REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})


add_library(core STATIC scene/Camera_Intrinsics.h scene/Camera_Intrinsics.cpp scene/Two_View.h scene/Two_View.cpp scene/Camera.cpp scene/Camera.h include/Core.h utils/Estimator.cpp utils/Estimator.h utils/Scene_IO.cpp utils/Scene_IO.h utils/Utilities.cpp utils/Utilities.h)

target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CERES_LIBRARIES} ${Boost_INCLUDE_DIRS})
target_link_libraries(core ${Boost_LIBRARIES} ${CERES_LIBRARIES})


