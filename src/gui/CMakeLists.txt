add_executable(fem_visualizer
    main.cpp
    mainwindow.cpp
    meshwidget.cpp
)

target_link_libraries(fem_visualizer PRIVATE
    fem_core
    Qt6::Core
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
)