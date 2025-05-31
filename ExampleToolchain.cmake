# the name of the target operating system
set(CMAKE_SYSTEM_NAME Darwin)
	
# Force compilation to M1 arm64
# NOTE: Requires all dependencies to be arm64!
set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)

# which compilers to use for C and C++
set(CMAKE_C_COMPILER /opt/homebrew/Cellar/llvm/19.1.7/bin/clang)
set(CMAKE_CXX_COMPILER /opt/homebrew/Cellar/llvm/19.1.7/bin/clang++)

# where is the target environment located
# NOTE: Boost is in /opt/homebrew/include/
# NOTE: oneAPI/TBB is in /usr/local
set(CMAKE_FIND_ROOT_PATH /opt/homebrew/Cellar/llvm/19.1.7/bin
	/opt/homebrew/include/
	/usr/local/include/)

# NOTE: Have not modified these from the cmake example
# adjust the default behavior of the FIND_XXX() commands:
# search programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)

# search headers and libraries in the target environment
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
