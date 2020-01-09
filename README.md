# TPM_converting_tools
C++ code implementations for converting raw read count to TPM values
<br>Requiring boost for unzipping, zipping, lexical_cast, etc, also based on [<tclap/CmdLine.h>](https://github.com/xguerin/tclap).
<br>**Matrix_read.cpp** is a program that accepts output folder *filter_matrix* from 10xGenomics.
<br>**TpmC.cpp** is a program that accepts files output from htseq-count. (with multithreading implemented)
