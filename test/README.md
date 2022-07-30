# Compilation

You can compile the package in Linux like this:

    {path}$mkdir build
    {path}$cd build
    {path}/build$cmake ..
    {path}/build$make
    {path}/build$./test3d_inter_1AC_plane ../data/data3d_inter_1AC_plane.bin
    {path}/build$./test3d_inter_2AC_plane ../data/data3d_inter_2AC_plane.bin
    {path}/build$./test3d_intra_2AC_plane ../data/data3d_intra_2AC_plane.bin
    {path}/build$./test4d_inter_2AC_vertical ../data/data4d_inter_2AC_vertical.bin
    {path}/build$./test4d_intra_2AC_vertical ../data/data4d_intra_2AC_vertical.bin