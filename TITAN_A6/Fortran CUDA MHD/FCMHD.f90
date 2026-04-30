! include "Storage_module.f90"
! include "cuf_kernel.cuf"
! include "cuf_solvers.cuf"
! include "Solvers.f90"



program MIK
    use STORAGE
    use MY_CUDA

    call Set_Storage()
    !call Fill_data()

    call CUDA_info()
    call CUDA_START_MGD()

    call Save_Storage()

end program MIK