
# setup - simulation
MLBUps = 200            # single NVidia V100 (rysy) - 200 MLBUps on D3Q27Q27
nGPUS = 100             # number of GPUs
avg_velocity = 0.01     # [lattice units/iteration]
domain_length = 512    # lattice units
fluid_pass = 10        # fluid pass through the domain
iterations = (1./avg_velocity)*domain_length*fluid_pass
lattice_size = domain_length**3
# lattice_size = 3*domain_length * domain_length/2  # consider a slice of minimal width: 3 lattice units
wall_clock_time_s = iterations*lattice_size/(MLBUps*1E6*nGPUS)

print(f"Wall clock time: {wall_clock_time_s/(60*24):.2f} [days] on {nGPUS} GPUs; "
      f"where domain is a {domain_length}^3 [lu] cube")


# setup - storage
DF_Bytes = 2*27*8           # D3Q27Q27 - double
DF_temp_Bytes = DF_Bytes    # temp array
q_Bytes = 2*27*(1/4)*8        # IBB distance from intersection
flag_Bytes = 4              # node type

macroscopic_rock_porosity_field_Bytes = 8   # for isotropic, non-homogeneous Darcy model

memory_per_lattice_node_MB = (DF_Bytes
                              +DF_temp_Bytes
                              +q_Bytes # without IBB
                              # +macroscopic_rock_porosity_field_Bytes
                              +flag_Bytes ) /1E6

required_gpu_memory_MB = memory_per_lattice_node_MB*lattice_size
print(f"Required GPU memory during simulation: {required_gpu_memory_MB/1E3:.2f} [GB]")

macroscopic_velocity_field_Bytes = 3 * 8    # Vx,Vy,Vz
macroscopic_fluid_density_field_Bytes = 8   # rho

mem_per_node_MB = (macroscopic_velocity_field_Bytes
                   + macroscopic_fluid_density_field_Bytes
                   + macroscopic_rock_porosity_field_Bytes
                   + flag_Bytes) / 1E6

total_storage_size_MB = mem_per_node_MB*lattice_size
print(f"Total HDD storage size for single snapshot: {total_storage_size_MB/1E3:.2f} [GB]")
