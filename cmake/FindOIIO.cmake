
find_path( OIIO_INCLUDE_DIR
    NAMES
        OpenImageIO/imageio.h
    PATHS
        ${OIIO_INCLUDE_PATH}
        ${OIIO_PATH}/include/
        /usr/include
        /usr/local/include
        /sw/include
        /opt/local/include
        DOC "The directory where OpenImageIO/imageio.h resides")

find_library( OIIO_LIBRARIES
    NAMES
         OIIO OpenImageIO
    PATHS
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