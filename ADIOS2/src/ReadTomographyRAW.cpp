/**
 * Part of the github.com/CFD-GO/TCLB
 * Example tomograph RAW reader intended for TCLB input - initial verison
 * @author mdzik
 * Example 
 * input data: https://www.digitalrocksportal.org/projects/202/origin_data/728/
 * call: ./ReadTomographyRAW ...../05_CC_K2_dry_vx_6x6x6_dim_734x734x425.raw ...../TCLB/example/adios/karman.adios.xml 734 734 425 32 1 0.500003
 * set 425 to 5 to reduce size for few slices, XY leave unchanged
 * see TCLB/example/adios for TCLB part example
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

typedef std::uint8_t FTYPE; //TODO: this should be char_to_bool to save A LOT memory

template<typename DTYPE>
int process(int argc, char*argv[]){
    
    unsigned long NX = atoi(argv[3]);
    unsigned long NY = atoi(argv[4]);
    unsigned long NZ = atoi(argv[5]);

    std::cout << "Will read RAW file: " << argv[1] <<  std::endl;
    std::cout << "Will read XML file: " << argv[2] <<  std::endl;
    std::cout << "Expected RAW shape: (" << NX << ", " << NY << ", "<< NZ << ")" << std::endl;
    

    int rank = 0, nproc = 1;

    float threshold = atof(argv[8]);

    std::fstream rf(argv[1], std::ios::in | std::ios::binary);
    if(!rf) {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }

    const float dmax = static_cast<float>( std::numeric_limits<DTYPE>::max() );
    const float dmin = static_cast<float>( std::numeric_limits<DTYPE>::min() );

    std::vector<DTYPE> buffer(NX*NY);

    std::vector<FTYPE> masked(NX*NY);
    

    std::cout << "Byte size: " <<  sizeof(DTYPE) << std::endl;
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
    adios2::ADIOS adios(argv[2], MPI_COMM_WORLD);
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
            if (!rf){
                std::cout << "error: only " << rf.gcount() << " could be read" << std::endl;
                writer.Close();
                return 1;
            }

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




int main(int argc, char *argv[])
{

    if (9 > argc) {
        std::cout << "Usage: ReadTomographRAW path_to_filename.raw path_to_adios.xml NX NY NZ BYTES SIGNED(0 or 1) RelativeThreshold" << std::endl;
        return 1;
    } 
    std::cout << "THIS IS EXPERIMENTAL/UNTESTED, DONT USE MPIRUN on IT!!" << std::endl;


    unsigned long NofBytes = atoi(argv[6]);
    bool IsSigned = atoi(argv[7]) == 1;
    float threshold = atof(argv[8]);

    std::cout << "Datatype: " << NofBytes << " bytes " << (IsSigned ? "signed":"unsigned") << std::endl;

    if (IsSigned && (NofBytes == 32))
        return process<std::int32_t>(argc, argv);

    if (!IsSigned && (NofBytes == 32))
        return process<std::uint32_t>(argc, argv);

    if (IsSigned && (NofBytes == 16))
        return process<std::int16_t>(argc, argv);

    if (!IsSigned && (NofBytes == 16))
        return process<std::uint16_t>(argc, argv);

    if (IsSigned && (NofBytes == 8))
        return process<std::int8_t>(argc, argv);

    if (!IsSigned && (NofBytes == 8))
        return process<std::uint8_t>(argc, argv);


}



