module cpl_exchange_mesh
  use o_DATA_TYPES
  implicit none
  save
  !**** general variables 
  integer                                     :: atm_nx,atm_ny
  real(kind=8), allocatable, dimension(:)     :: atm_lon, atm_lat 

  character*100                               :: atm_mesh_path='./'
  real(kind=8)                                :: CplDx
  real(kind=8)                                :: CplDy


  type(sparse_matrix)                         :: atm_2_oce
  type(sparse_matrix)                         :: oce_2_atm
  integer, dimension(:,:),        allocatable   :: oce_2_atm_mask
  real(kind=8), dimension(:),     allocatable   :: oce_fld
  real(kind=8), dimension(:,:),   allocatable   :: atm_fld
  real(kind=8), dimension(:,:),   allocatable   :: a2o_fcorr_stat  !flux correction statistics for the output
  
  real(kind=8)                                :: time_send(2), time_recv(2)
  integer                                     :: o2a_call_count=0
 
  real(kind=8), allocatable, dimension(:,:)   :: cplsnd
  
end module cpl_exchange_mesh

module cpl_config_param
   implicit none
   save

  ! Presetting of coupling parameters, should later be delivered by ECHAM xml - file

  ! *** enable / disable OASIS4 coupling
  real(kind=8)                  :: CplLonMin=-180.
  real(kind=8)                  :: CplLonMax=180.
  real(kind=8)                  :: CplLatMin=-90.
  real(kind=8)                  :: CplLatMax=90. 
  integer                       :: CplIDim=960!192   !320
  integer                       :: CplJDim=480!96    !160
  namelist /cpl_params/ CplLonMin,CplLonMax,CplLatMin,CplLatMax,CplIDim,CplJDim

end module cpl_config_param
