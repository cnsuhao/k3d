#
# Try to find COLLADA libraries, and include path.
# 
#Unix configuration
IF(UNIX AND NOT APPLE)
	FIND_PATH (COLLADA_DAE_INCLUDE_PATH dae.h
		/usr/local/include/colladadom
		/usr/include/colladadom
	)
	FIND_PATH (COLLADA_DOM_INCLUDE_PATH /dom
		/usr/local/include/colladadom/1.4
		/usr/include/colladadom/1.4
	)
	
	FIND_LIBRARY (COLLADA_LIBRARY
		NAMES collada14dom
		PATHS
		/usr/local/lib
		/usr/lib
	)
	
ENDIF(UNIX AND NOT APPLE)
	
#Apple configuration
IF(UNIX AND APPLE)
		FIND_PATH (COLLADA_DAE_INCLUDE_PATH dae.h
		/Library/Frameworks/Collada14Dom.framework
		)
		FIND_PATH (COLLADA_DOM_INCLUDE_PATH dom/
		/Library/Frameworks/Collada14Dom.framework/1.4
		)
		FIND_LIBRARY (COLLADA_LIBRARY
		NAMES Collada14Dom
		PATHS
		/Library/Frameworks/Collada14Dom.framework
		)
ENDIF(UNIX AND APPLE)

#Windows configuration
IF(WIN32)
		FIND_PATH (COLLADA_DAE_INCLUDE_PATH dae.h
		C:\Program Files\domcollada\include

		)
		FIND_PATH (COLLADA_DOM_INCLUDE_PATH dom/
		C:\Program Files\domcollada\include\1.4
		)
	
		FIND_LIBRARY (COLLADA_LIBRARY
			NAMES Collada14Dom
			PATHS
		C:\Program Files\domcollada\build\vc8-1.4
		)
ENDIF(WIN32)

IF (COLLADA_DAE_INCLUDE_PATH AND COLLADA_DOM_INCLUDE_PATH AND COLLADA_LIBRARY)
	MARK_AS_ADVANCED(COLLADA_DAE_INCLUDE_PATH)
	MARK_AS_ADVANCED(COLLADA_DOM_INCLUDE_PATH)	
	MARK_AS_ADVANCED (COLLADA_FOUND COLLADA_LIBRARY)
	SET (COLLADA_FOUND TRUE)
ELSE (COLLADA_DAE_INCLUDE_PATH AND COLLADA_DOM_INCLUDE_PATH AND COLLADA_LIBRARY)
	SET (COLLADA_FOUND FALSE)
ENDIF (COLLADA_DAE_INCLUDE_PATH AND COLLADA_DOM_INCLUDE_PATH AND COLLADA_LIBRARY)

SET (COLLADA_DAE_INCLUDE_DIR ${COLLADA_DAE_INCLUDE_PATH})
SET (COLLADA_DOM_INCLUDE_DIR ${COLLADA_DOM_INCLUDE_PATH})
SET (COLLADA_LIBRARIES ${COLLADA_LIBRARY})



