message(STATUS "")
message(STATUS ">>> Setting up 'mcmc' directory.")


ADD_EXECUTABLE(hamMCpwa hamMCpwa.cc)
TARGET_LINK_LIBRARIES(hamMCpwa rootpwa)
TARGET_LINK_LIBRARIES(hamMCpwa ${PWA_LIBS})
TARGET_LINK_LIBRARIES(hamMCpwa ${ROOT_LIBS})

ADD_EXECUTABLE(hamMCtest hamMCtest.cc)
TARGET_LINK_LIBRARIES(hamMCtest rootpwa)
TARGET_LINK_LIBRARIES(hamMCtest ${PWA_LIBS})
TARGET_LINK_LIBRARIES(hamMCtest ${ROOT_LIBS})

ADD_EXECUTABLE(mcmc2fitbin mcmc2fitbin.cc)
TARGET_LINK_LIBRARIES(mcmc2fitbin rootpwa)
TARGET_LINK_LIBRARIES(mcmc2fitbin ${PWA_LIBS})
TARGET_LINK_LIBRARIES(mcmc2fitbin ${ROOT_LIBS})