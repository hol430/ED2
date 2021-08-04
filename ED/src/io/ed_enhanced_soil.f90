Module ed_enhanced_soil
  implicit none

Contains
  subroutine wise_soil_database(headname, nsite, nlandsea, lat, lon, classout, pctout,  &
       orgc_out, c2n_out)
    use hdf5_utils, only: shdf5_open_f, shdf5_close_f, shdf5_irec_f, shdf5_info_f
    use soil_coms, only: find_soil_class
    implicit none
!  subroutine wise_soil_database(trim(soil_database(ifm)),maxsite,npoly,'soil_text'                &
!                        ,lat_list,lon_list,ntext_soil_list,ipcent_soil)
    character(len=*), intent(in) :: headname
    integer, intent(in) :: nsite
    integer, intent(in) :: nlandsea
    real, dimension(3,nlandsea), intent(in) :: lat, lon
    integer, dimension(nsite,nlandsea), intent(out) :: classout
    real, dimension(nsite,nlandsea), intent(out) :: pctout
    real, dimension(nsite,nlandsea), intent(out) :: orgc_out
    real, dimension(nsite,nlandsea), intent(out) :: c2n_out
    logical :: lexist
    integer :: nfile
    character(len=256) :: fname
    integer :: ndims
    integer, dimension(2) :: idims
    real, dimension(:,:), allocatable :: clay, sand, orgc, c2n
    integer :: ilandsea,ilat,ilon

    classout = 0
    pctout = 0.
    orgc_out = 0.
    c2n_out = 0.

    inquire(file=trim(headname),exist=lexist)
    if(.not.lexist)then
       write(*,*)'File does not exist in wise_soil_database'
       write(*,*)trim(headname)
       stop
    endif

    open(12,file=trim(headname),form='formatted',status='old')
    read(12,*)nfile

    read(12,'(a)')fname
    call shdf5_open_f(trim(fname),'R')
    call shdf5_info_f('CLAY',ndims,idims)
    allocate(clay(idims(1),idims(2)))
    allocate(sand(idims(1),idims(2)))
    allocate(orgc(idims(1),idims(2)))
    allocate(c2n(idims(1),idims(2)))
    call shdf5_irec_f(ndims,idims,'CLAY',rvara=clay)
    call shdf5_irec_f(ndims,idims,'SAND',rvara=sand)
    call shdf5_irec_f(ndims,idims,'ORGC',rvara=orgc)
    call shdf5_irec_f(ndims,idims,'CN',rvara=c2n)
    do ilandsea = 1, nlandsea
       ilat = int((lat(1,ilandsea)+90.0)*2.0)+1
       ilon = int((lon(1,ilandsea)+180.0)*2.0)+1
       if(sand(ilon,ilat) > -1)then
          classout(1,ilandsea) = find_soil_class(sand(ilon,ilat)*0.01,clay(ilon,ilat)*0.01)
          pctout(1,ilandsea) = 1.
          orgc_out(1,ilandsea) = orgc(ilon,ilat)
          c2n_out(1,ilandsea) = c2n(ilon,ilat)
       endif
    enddo

    deallocate(clay)
    deallocate(sand)
    deallocate(orgc)
    deallocate(c2n)
    call shdf5_close_f()

    close(12)

    return
  end subroutine wise_soil_database

end Module ed_enhanced_soil
