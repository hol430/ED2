Module ed_enhanced_soil
  implicit none

Contains

  subroutine wise_soil_database(headname, nsite, nlandsea, lat, lon, classout, pctout,  &
       orgc_out, c2n_out, apa_out, lab_out, occ_out, org_out, sec_out)

    use hdf5_utils, only: shdf5_open_f, shdf5_close_f, shdf5_irec_f, shdf5_info_f
    use soil_coms, only: find_soil_class

    implicit none

    character(len=*), intent(in) :: headname
    integer, intent(in) :: nsite
    integer, intent(in) :: nlandsea
    real, dimension(3,nlandsea), intent(in) :: lat, lon
    integer, dimension(nsite,nlandsea), intent(out) :: classout
    real, dimension(nsite,nlandsea), intent(out) :: pctout
    real, dimension(nsite,nlandsea), intent(out) :: orgc_out
    real, dimension(nsite,nlandsea), intent(out) :: c2n_out
    real, dimension(nsite,nlandsea), intent(out) :: apa_out, lab_out, occ_out, org_out, sec_out
    logical :: lexist
    integer :: nfile
    character(len=256) :: fname
    integer :: ndims
    integer, dimension(2) :: idims
    real, dimension(:,:), allocatable :: clay, sand, orgc, c2n, apa, lab, occ, org, sec
    integer :: ilandsea,ilat,ilon

    classout = 0
    pctout = 0.
    orgc_out = 0.
    c2n_out = 0.

    apa_out = 0.
    lab_out = 0.
    occ_out = 0.
    org_out = 0.
    sec_out = 0.

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

    read(12,'(a)')fname
    call shdf5_open_f(trim(fname),'R')
    call shdf5_info_f('apa',ndims,idims)
    allocate(apa(idims(1),idims(2)))
    allocate(lab(idims(1),idims(2)))
    allocate(occ(idims(1),idims(2)))
    allocate(org(idims(1),idims(2)))
    allocate(sec(idims(1),idims(2)))
    call shdf5_irec_f(ndims,idims,'apa',rvara=apa)
    call shdf5_irec_f(ndims,idims,'lab',rvara=lab)
    call shdf5_irec_f(ndims,idims,'occ',rvara=occ)
    call shdf5_irec_f(ndims,idims,'org',rvara=org)
    call shdf5_irec_f(ndims,idims,'sec',rvara=sec)
    do ilandsea = 1, nlandsea
       ilat = int((lat(1,ilandsea)+90.0)*2.0)+1
       ilon = int((lon(1,ilandsea)+180.0)*2.0)+1
       if(apa(ilon,ilat)==apa(ilon,ilat))then
          apa_out(1,ilandsea) = apa(ilon,ilat)
          lab_out(1,ilandsea) = lab(ilon,ilat)
          occ_out(1,ilandsea) = occ(ilon,ilat)
          org_out(1,ilandsea) = org(ilon,ilat)
          sec_out(1,ilandsea) = sec(ilon,ilat)
       elseif(apa(ilon,ilat+1)==apa(ilon,ilat+1))then
          apa_out(1,ilandsea) = apa(ilon,ilat+1)
          lab_out(1,ilandsea) = lab(ilon,ilat+1)
          occ_out(1,ilandsea) = occ(ilon,ilat+1)
          org_out(1,ilandsea) = org(ilon,ilat+1)
          sec_out(1,ilandsea) = sec(ilon,ilat+1)
       elseif(apa(ilon+1,ilat)==apa(ilon+1,ilat))then
          apa_out(1,ilandsea) = apa(ilon+1,ilat)
          lab_out(1,ilandsea) = lab(ilon+1,ilat)
          occ_out(1,ilandsea) = occ(ilon+1,ilat)
          org_out(1,ilandsea) = org(ilon+1,ilat)
          sec_out(1,ilandsea) = sec(ilon+1,ilat)
       elseif(apa(ilon,ilat-1)==apa(ilon,ilat-1))then
          apa_out(1,ilandsea) = apa(ilon,ilat-1)
          lab_out(1,ilandsea) = lab(ilon,ilat-1)
          occ_out(1,ilandsea) = occ(ilon,ilat-1)
          org_out(1,ilandsea) = org(ilon,ilat-1)
          sec_out(1,ilandsea) = sec(ilon,ilat-1)
       elseif(apa(ilon-1,ilat)==apa(ilon-1,ilat))then
          apa_out(1,ilandsea) = apa(ilon-1,ilat)
          lab_out(1,ilandsea) = lab(ilon-1,ilat)
          occ_out(1,ilandsea) = occ(ilon-1,ilat)
          org_out(1,ilandsea) = org(ilon-1,ilat)
          sec_out(1,ilandsea) = sec(ilon-1,ilat)
       elseif(apa(ilon,ilat+2)==apa(ilon,ilat+2))then
          apa_out(1,ilandsea) = apa(ilon,ilat+2)
          lab_out(1,ilandsea) = lab(ilon,ilat+2)
          occ_out(1,ilandsea) = occ(ilon,ilat+2)
          org_out(1,ilandsea) = org(ilon,ilat+2)
          sec_out(1,ilandsea) = sec(ilon,ilat+2)
       elseif(apa(ilon+1,ilat+1)==apa(ilon+1,ilat+1))then
          apa_out(1,ilandsea) = apa(ilon+1,ilat+1)
          lab_out(1,ilandsea) = lab(ilon+1,ilat+1)
          occ_out(1,ilandsea) = occ(ilon+1,ilat+1)
          org_out(1,ilandsea) = org(ilon+1,ilat+1)
          sec_out(1,ilandsea) = sec(ilon+1,ilat+1)
       elseif(apa(ilon+2,ilat)==apa(ilon+2,ilat))then
          apa_out(1,ilandsea) = apa(ilon+2,ilat)
          lab_out(1,ilandsea) = lab(ilon+2,ilat)
          occ_out(1,ilandsea) = occ(ilon+2,ilat)
          org_out(1,ilandsea) = org(ilon+2,ilat)
          sec_out(1,ilandsea) = sec(ilon+2,ilat)
       elseif(apa(ilon+1,ilat-1)==apa(ilon+1,ilat-1))then
          apa_out(1,ilandsea) = apa(ilon+1,ilat-1)
          lab_out(1,ilandsea) = lab(ilon+1,ilat-1)
          occ_out(1,ilandsea) = occ(ilon+1,ilat-1)
          org_out(1,ilandsea) = org(ilon+1,ilat-1)
          sec_out(1,ilandsea) = sec(ilon+1,ilat-1)
       elseif(apa(ilon,ilat-2)==apa(ilon,ilat-2))then
          apa_out(1,ilandsea) = apa(ilon,ilat-2)
          lab_out(1,ilandsea) = lab(ilon,ilat-2)
          occ_out(1,ilandsea) = occ(ilon,ilat-2)
          org_out(1,ilandsea) = org(ilon,ilat-2)
          sec_out(1,ilandsea) = sec(ilon,ilat-2)
       elseif(apa(ilon-1,ilat-1)==apa(ilon-1,ilat-1))then
          apa_out(1,ilandsea) = apa(ilon-1,ilat-1)
          lab_out(1,ilandsea) = lab(ilon-1,ilat-1)
          occ_out(1,ilandsea) = occ(ilon-1,ilat-1)
          org_out(1,ilandsea) = org(ilon-1,ilat-1)
          sec_out(1,ilandsea) = sec(ilon-1,ilat-1)
       elseif(apa(ilon-2,ilat)==apa(ilon-2,ilat))then
          apa_out(1,ilandsea) = apa(ilon-2,ilat)
          lab_out(1,ilandsea) = lab(ilon-2,ilat)
          occ_out(1,ilandsea) = occ(ilon-2,ilat)
          org_out(1,ilandsea) = org(ilon-2,ilat)
          sec_out(1,ilandsea) = sec(ilon-2,ilat)
       elseif(apa(ilon-1,ilat+1)==apa(ilon-1,ilat+1))then
          apa_out(1,ilandsea) = apa(ilon-1,ilat+1)
          lab_out(1,ilandsea) = lab(ilon-1,ilat+1)
          occ_out(1,ilandsea) = occ(ilon-1,ilat+1)
          org_out(1,ilandsea) = org(ilon-1,ilat+1)
          sec_out(1,ilandsea) = sec(ilon-1,ilat+1)
       endif
    enddo

    deallocate(apa)
    deallocate(lab)
    deallocate(occ)
    deallocate(org)
    deallocate(sec)
    call shdf5_close_f()

    close(12)

    return
  end subroutine wise_soil_database

end Module ed_enhanced_soil
