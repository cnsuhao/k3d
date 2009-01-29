SET(K3D_GNOME_VFS_FOUND 1)

INCLUDE(K3DFindPkgConfig)
PKG_CHECK_MODULES(GNOMEVFS gnome-vfs-2.0)

IF(GNOMEVFS_FOUND)
	SET(K3D_GNOME_VFS_INCLUDE_DIRS
		${GNOMEVFS_INCLUDE_DIRS}
		)

	SET(K3D_GNOME_VFS_LIB_DIRS
		${GNOMEVFS_LIBRARY_DIRS}
		)

	SET(K3D_GNOME_VFS_LIBS
		${GNOMEVFS_LIBRARIES}
		)

	SET(K3D_GNOME_VFS_FOUND 1)
ENDIF(GNOMEVFS_FOUND)

