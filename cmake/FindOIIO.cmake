
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

# From OIIO changelog - Version: 2.3 (1 Sept 2021)
#   "All the utility classes are now in libOpenImageIO_Util only and libOpenImageIO depends on and links to libOpenImageIO_Util, 
#   rather than the utility classes being defined separately in both libraries."
find_library( OIIO_LIBRARIES_UTIL
    NAMES
         OpenImageIO_Util
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
        DOC "The OIIO library utility")
