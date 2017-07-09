INCLUDEPATH += $${_PRO_FILE_PWD_}/../

INCLUDEPATH *= /usr/include
INCLUDEPATH *= /usr/local/include

#Generate binary executables at the root of the build folder
DESTDIR = $${OUT_PWD}/../
DEFINES += OM_STATIC_BUILD

win32{
    #LIBS += "$${OUT_PWD}/../OpenMesh/libOpenMesh.a"
    LIBS += "$${OUT_PWD}/../OpenMesh/OpenMesh.lib"
    LIBS += -L$$PWD/ANN/lib/ -lANN
    #LIBS += opengl32.lib
    LIBS += -lglut32
    #LIBS += -LC:\glut
}

unix{
    LIBS += -L$${OUT_PWD}/../OpenMesh -lOpenMesh

    macx{
        LIBS += -framework OpenGL
    }
    else{
        LIBS += -lGLU
    }
}
