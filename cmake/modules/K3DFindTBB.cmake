SET(K3D_TBB_FOUND 0)

FIND_PATH(K3D_TBB_INCLUDE_DIR parallel_for.h
	/usr/local/include
	/usr/include
	PATH_SUFFIXES tbb
	)

FIND_LIBRARY(K3D_TBB_LIBRARY tbb
	/usr/local/lib
	/usr/lib
	)

IF(K3D_TBB_INCLUDE_DIR AND K3D_TBB_LIBRARY)
	SET(K3D_TBB_FOUND 1)
ENDIF(K3D_TBB_INCLUDE_DIR AND K3D_TBB_LIBRARY)

MARK_AS_ADVANCED(
	K3D_TBB_INCLUDE_DIR
	K3D_TBB_LIBRARY
	)

