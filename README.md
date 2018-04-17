# tarbgz: externally-indexed random-access compressed tar archives

This tool uses bgzip, a random-access backward-compatible extension of gzip, to allow individual files in tar archives to be accessed efficiently.

The tar file is first created with `tar`, and then compressed with `bgzip`, producing a .tar.gz. This file is then indexed with `tarbgz.py`, producing a second .tar.gz.index file. The index contains a complete listing of the files in the tar archive, and can be used to search for files without the tar archive being available. Then, when a desired file is found, the file can be extracted from the tar archive using a well-defined byte range read on the compressed archive.

Although not implemented yet, this tool is intended to support sub-archive reads for files stored in Amazon Glacier and other similar offline or nearline archiving systems, in which it is desireable to coalesce small files into large single archives, but it is also desireable to retain the ability to retrieve individual files.
