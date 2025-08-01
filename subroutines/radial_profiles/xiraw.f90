!-----------------------------------------------------------------------
function xiraw(re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,mudisk,cosd,gsd)
    ! In: re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts
    ! Out: logxiraw,gsd
    implicit none
    integer         , intent(in)  :: ndelta,npts
    double precision, intent(in)  :: re, spin, h, honr, rmin, mudisk, cosd
    double precision, intent(in)  :: rlp(ndelta),dcosdr(ndelta)
    double precision, intent(out) :: gsd
    integer                       :: kk,get_index
    double precision cosfac,interper,xiraw,dareafac,dglpfacthick,newtex

    !Calculate source to disc blueshift at this radius
    gsd = dglpfacthick(re,spin,h,mudisk, cosd)
    !Find the rlp bin that corresponds to re
    kk = get_index(rlp,ndelta,re,rmin,npts)
    !Interpolate to get |d\cos\delta/dr| at r=re
    cosfac = interper(rlp,dcosdr,ndelta,re,kk)
    !Extrapolate to Newtonian if needs be
    if( kk .eq. npts ) cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
    !Now can do the calculation
    xiraw = gsd**2 * cosfac / dareafac(re,spin) 
    return
end function xiraw  
!-----------------------------------------------------------------------
