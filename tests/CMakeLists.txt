# Add target to download input file
add_test(NAME download_sample
    WORKING_DIRECTORY ${INSTALLDIR}
    COMMAND curl -OL https://github.com/KIT-CMS/CROWNTestingSamples/raw/main/nanoAOD.root)
set_tests_properties(download_sample PROPERTIES FIXTURES_SETUP download_sample LABELS ntuple_download)

# Generate a test for each generated target
foreach(TARGET_NAME ${TARGET_NAMES})
    message(STATUS "Add test for target ${TARGET_NAME}")
    add_test(NAME crown_ntuple_${TARGET_NAME}
             WORKING_DIRECTORY ${INSTALLDIR}
             COMMAND ${TARGET_NAME} output_${TARGET_NAME}.root nanoAOD.root)
    set_tests_properties(crown_ntuple_${TARGET_NAME} PROPERTIES FIXTURES_REQUIRED download_sample LABELS ntuple_run)
endforeach()


# add friend tree test
# Add target to download input file
add_test(NAME friends_download_sample
    WORKING_DIRECTORY ${INSTALLDIR}
    COMMAND curl -OL https://github.com/KIT-CMS/CROWNTestingSamples/raw/main/CrownNtuple.root)
set_tests_properties(friends_download_sample PROPERTIES FIXTURES_SETUP friends_download_sample LABELS friends_download)

# Generate a test for each generated target
foreach(TARGET_NAME ${TARGET_NAMES})
    message(STATUS "Add test for target ${TARGET_NAME}")
    add_test(NAME friends_ntuple_${TARGET_NAME}
             WORKING_DIRECTORY ${INSTALLDIR}
             COMMAND ${TARGET_NAME} output_${TARGET_NAME}.root CrownNtuple.root)
    set_tests_properties(friends_ntuple_${TARGET_NAME} PROPERTIES FIXTURES_REQUIRED friends_download_sample LABELS friends_run)
endforeach()
