
find_path( OIIO_INCLUDE_DIR
    NAMES
        OpenImageIO/imageio.h
    PATHS
        ${OIIO_INCLUDE_PATH}
        ${OIIO_PATH}/include
        /usr/include
        /usr/local/include
        /sw/include
        /opt/local/include
        DOC "The directory where OpenImageIO/imageio.h resides")

if (DEFINED ENV{RP_oiio})
    set(OIIO_INCLUDE_DIR "$ENV{RP_oiio}/include")
endif()

find_library( OIIO_LIBRARIES
    NAMES
         OIIO OpenImageIO
    PATHS
        $ENV{RP_oiio}/${CMAKE_INSTALL_LIBDIR}
        ${OIIO_LIBRARY_PATH}
        ${OIIO_PATH}/lib64/
        ${OIIO_PATH}/lib/
        /usr/lib64
        /usr/lib
        /usr/local/lib64
        /usr/local/lib
        /sw/lib
        /opt/local/lib
        DOC "The OIIO library")

find_library( OIIO_Util_LIBRARIES
    NAMES
         OIIO_Util OpenImageIO_Util
    PATHS
        $ENV{RP_oiio}
        ${OIIO_LIBRARY_PATH}
        ${OIIO_PATH}/lib64/
        ${OIIO_PATH}/lib/
        /usr/lib64
        /usr/lib
        /usr/local/lib64
        /usr/local/lib
        /sw/lib
        /opt/local/lib
        DOC "The OIIO Util library")
