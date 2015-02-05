project(PS1080Test)

file(GLOB ps1080test_SRC
    "*.h"
	"*.hpp"
    "*.cpp"
)

add_executable(PS1080Test ${ps1080test_SRC})
target_link_libraries(PS1080Test ${SDL2_LIBRARY} ${OPENNI2_LIBRARY})

if(WIN32)
    add_custom_command(TARGET PS1080Test POST_BUILD  
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OPENNI2_REDIST_DIR}/OpenNI2.dll"      
        $<TARGET_FILE_DIR:PS1080Test>)             

    add_custom_command(TARGET PS1080Test POST_BUILD  
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${SDL2_INCLUDE_DIR}/../bin/SDL2.dll"      
        $<TARGET_FILE_DIR:PS1080Test>)    
endif()
