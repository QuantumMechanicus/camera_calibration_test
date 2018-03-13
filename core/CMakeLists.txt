FIND_PACKAGE(Ceres REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options filesystem REQUIRED)
FIND_PACKAGE(Glog REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})


add_library(core STATIC scene/Intrinsics.h scene/Two_View.h scene/Camera.h include/Core.h utils/Scene_IO.cpp utils/Scene_IO.h utils/Utilities.cpp utils/Utilities.h interfaces/Abstract_Estimator.h interfaces/ICamera.h interfaces/INode.h interfaces/IEdge.h interfaces/ITwo_View.h scene/Scene.h interfaces/IGraph.h interfaces/IScene.h interfaces/IIntrinsics.h utils/Local_Parametrization_SO3.h utils/Local_Parametrization_Sphere.h)

target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CERES_LIBRARIES} ${Boost_INCLUDE_DIRS})
target_link_libraries(core ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GLOG_LIBRARIES})


