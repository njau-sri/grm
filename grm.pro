TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_LFLAGS += -static

SOURCES += main.cpp \
    appgrm.cpp \
    util.cpp \
    vcfio.cpp \
    cmdline.cpp

HEADERS += \
    appgrm.h \
    main.h \
    strsplit.h \
    util.h \
    vcfio.h \
    cmdline.h
