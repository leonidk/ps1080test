project(PS1080Test)
add_executable(PS1080Test main.cpp detect.cpp)
target_link_libraries(PS1080Test ${SDL2_LIBRARY} ${OPENNI2_LIBRARY})

add_custom_command(TARGET PS1080Test POST_BUILD  
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OPENNI2_REDIST_DIR}/OpenNI2.dll"      
        $<TARGET_FILE_DIR:PS1080Test>)             

add_custom_command(TARGET PS1080Test POST_BUILD  
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${SDL2_INCLUDE_DIR}/../bin/SDL2.dll"      
        $<TARGET_FILE_DIR:PS1080Test>)    