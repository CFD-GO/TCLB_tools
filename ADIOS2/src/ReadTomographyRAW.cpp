/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write a global array from multiple processors.
 *
 * A global array is an N-dimensional array. A process can write a sub-array
 * into the global array by stating the N-dimensional offset and the size of
 * the sub-array. At reading, one can read back any portion of the array
 * regardless of how many processors wrote that data.
 *
 * Processes are NOT required
 * - to stay in the boundaries of the global dimensions. However, one will not
 * be able to read back data outside of the boundaries.
 * - to fill the whole global array, i.e. one can leave holes in it. At reading,
 * one will get the fill-value set for the array for those coordinates that
 * are not written by any process.
 *
 * The global dimensions of a global array MUST NOT change over time.
 * If they are, then the array should be handled as a local array. Of course, if
 * only a single output step is written to a file, that still shows up at
 * reading as a global array.
 *
 * The decomposition of the array across the processes, however, can change
 * between output steps.
 *
 * Created on: Jun 2, 2017
 *      Author: pnorbert
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

typedef std::int32_t DTYPE;


typedef std::uint8_t FTYPE; //TODO: this should be char_to_bool to save A LOT memory

int main(int argc, char *argv[])
{


    int rank = 0, nproc = 1;


    unsigned long NX = 734;
    unsigned long NY = 734;
    unsigned long NZ = 100;
    unsigned long NofBytes = 32;
    bool IsSigned = true;
    float threshold = 0.500002;



    std::fstream rf("/home/mdzik/Pobrane/05_CC_K2_dry_vx_6x6x6_dim_734x734x425.raw", std::ios::in | std::ios::binary);
    if(!rf) {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }

    const float dmax = static_cast<float>( std::numeric_limits<DTYPE>::max() );
    const float dmin = static_cast<float>( std::numeric_limits<DTYPE>::min() );

    std::vector<DTYPE> buffer(NX*NY);

    std::vector<FTYPE> masked(NX*NY);
    

    std::cout << "Bite size: " <<  sizeof(DTYPE) << std::endl;
    std::cout << "Expected filesize: (ls -la) " <<  NX*NY*NZ*sizeof(DTYPE) << std::endl;


    rf.read((char*)buffer.data(), NX*NY*sizeof(DTYPE));
    if (rf)
        std::cout << "all characters read successfully.";
    else
        std::cout << "error: only " << rf.gcount() << " could be read" << std::endl;

    rf.seekg (0, rf.beg);


    const DTYPE threshold_dtyped = static_cast<DTYPE>( dmin + threshold*(dmax-dmin) );
    std::cout << "Set threshold: " << threshold << "-> VALUE < " << threshold_dtyped << std::endl;

    unsigned int white = 0;
    unsigned int all = 0;
    DTYPE max = std::numeric_limits<DTYPE>::min()+1;
    DTYPE min = std::numeric_limits<DTYPE>::max()-1;
    
    std::cout << "Reading one Z slice" << std::endl;

    for (unsigned long i =0; i<NX*NY;i++) {
        masked[i] = (FTYPE)(buffer[i] < threshold_dtyped);

        all++;
        //std::cout << buffer[i] <<  " MIN: " << min << " MAX:" << max  << std::endl;
        if (buffer[i] > max) max = buffer[i];
        if (buffer[i] < min) min = buffer[i];

        if (masked[i]) white++;
        }

    std::cout << "All: " << all << " V < THRESHOLD (return 1) " << white << std::endl;

    std::cout << "MIN: " << min << " MAX:" << max <<std::endl;

#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
    const int NSTEPS = 5;

#if ADIOS2_USE_MPI
    adios2::ADIOS adios("./adios.xml", MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif


    try
    {
        // Get io settings from the config file or
        // create one with default settings here
        adios2::IO io = adios.DeclareIO("FlaggedGeometry");

        // should work but id doesnt
        // const std::string extent = "0 0 0 " + std::to_string(NX) + " " + std::to_string(NY) + " 0";
        // const std::string imageData = R"(
        //     <?xml version="1.0"?>
        //     <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
        //     <ImageData WholeExtent=")" + extent + R"(" Origin="0 0 0" Spacing="1 1 1">
        //         <Piece Extent=")" + extent + R"(">
        //         <CellData Scalars="FlaggedGeometry">
        //             <DataArray Name="FlaggedGeometry" />
        //             <DataArray Name="TIME">
        //                 step
        //             </DataArray>
        //         </CellData>
        //         </Piece>
        //     </ImageData>
        //     </VTKFile>)";

        // io.DefineAttribute<std::string>("vtk.xml", imageData);        

        adios2::Variable<size_t> var_step = io.DefineVariable<size_t>("step");

        /*
         * Define global array: type, name, global dimensions
         * The local process' part (start, count) can be defined now or later
         * before Write().
         */

        adios2::Variable<FTYPE> varGlobalArray = io.DefineVariable<FTYPE>("GlobalArray", {NX,NY});

        // Open file. "w" means we overwrite any existing file on disk,
        // but Advance() will append steps to the same file.
        adios2::Engine writer = io.Open("globalArray", adios2::Mode::Write);

        for (size_t step = 0; step < NZ; step++)
        {
            std::cout << "Reading Z-slice no:" << step <<std::endl;

            writer.BeginStep();
            
            writer.Put<size_t>(var_step, &step);

            rf.read((char*)buffer.data(), NX*NY*sizeof(DTYPE));
            if (!rf)
                std::cout << "error: only " << rf.gcount() << " could be read" << std::endl;


            for (size_t i = 0; i < NX*NY;i++) {
                masked[i] = buffer[i] > threshold_dtyped;
            }
            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            // adios2::SelectionBoundingBox sel();
            varGlobalArray.SetSelection(adios2::Box<adios2::Dims>( {0,0}, {NX,NY}));


            writer.Put<FTYPE>(varGlobalArray, (FTYPE*)masked.data());

            // Indicate we are done for this step.
            // Disk I/O will be performed during this call unless
            // time aggregation postpones all of that to some later step
            writer.EndStep();
        }

        // Called once: indicate that we are done with this output for the run
        writer.Close();

    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cout << "Invalid argument exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::ios_base::failure &e)
    {
        if (rank == 0)
        {
            std::cout << "System exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cout << "Exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif


    return 0;
}
