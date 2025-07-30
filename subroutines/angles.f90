! !-----------------------------------------------------------------------
!       function dinang(a,r,h,mus)
! ! mus = \cos\delta
! ! Calculates cosine of angle of incidence for a photon hitting the
! ! disk mid-plane from a lamppost source
! ! mui = pi_alpha n^alpha / ( pi_mu u^mu )
! !     = pi_alpha n^alpha / ( -u^t )
! !     = ( ql / r ) / ( u^t )
!         implicit none
!       double precision dinang,a,r,h,mus
!       double precision DeltaL,Delta,ql,Omega
!       !Define useful combinations
!       DeltaL  = h**2 - 2*h + a**2
!       Delta   = r**2 - 2*r + a**2
!       !Carter's constant of incoming photon      
!       ql      = (1.0-mus**2) * (h**2+a**2)**2/DeltaL - a**2
!       ql      = sqrt( ql )
!       !Orbital angular velocity of disk material
!       Omega   = sign(a,1.d0) * 1.0 / ( r**1.5 + abs(a) )
!       ! Calculate mui = (ql/r)/ut
!       dinang = 1-2./r + 4*a*Omega/r
!       dinang = dinang - ( (r**2+a**2)**2 - a**2*Delta )/r**2 * Omega**2
!       dinang = sqrt( dinang ) * ql / r
!       return
!       end
! !-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
      function dinang(a,r,h,mus)
! mus = \cos\delta
! Calculates cosine of angle of incidence for a photon hitting the
! disk mid-plane from a lamppost source
! mui = pi_alpha n^alpha / ( pi_mu u^mu )
!     = pi_alpha n^alpha / ( -u^t )
!     = ( ql / r ) / ( u^t )
      use blcoordinate
      use isco
      implicit none
      double precision dinang,a,r,h,mus,f1234(4),pcros
      double precision Delta,ql,Omega,kk,hh,vr,vt,vphi,kr,velocity(3),lambda
      double precision vt_, vr_, kinematic_energy,angular_momentum
      double precision phi_r,time_r,aff_r,r_coord
      double precision r1,mu1,phi1,t1,sigma1
      double precision r2,mu2,phi2,t2,sigma2
      integer t_r1, t_r2
      velocity = 0.0D0
      !Define useful combinations
      Delta   = r**2 - 2*r + a**2
      !Carter's constant of incoming photon      
      call initialdirection(mus,sqrt(1-mus**2.0),0.d0,0.d0,1.d0,a,h,velocity ,lambda,ql,f1234)
      if(r.gt.risco)then
         !Orbital angular velocity of disk material
         Omega   = sign(a,1.d0) * 1.0 / ( r**1.5 + abs(a) )
         ! Calculate mui = (ql/r)/ut
         dinang = 1-2./r + 4*a*Omega/r
         dinang = dinang - ( (r**2+a**2)**2 - a**2*Delta )/r**2 * Omega**2
         dinang = sqrt( dinang ) * sqrt(ql) / r
      else
         pcros = Pemdisk(f1234,lambda,ql,0.d0,1.d0,a,h,1.d0,0.d0,1.d10,1.d0+sqrt(1.d0-a**2))
         if(pcros.lt.0.d0)then
            dinang=0.d0
         else
            vt_ = vt(r,a)
            vr_ = vr(r)
            kr=sqrt((r**2.0+a**2.0)**2.0-Delta*(ql+a**2.0))/Delta      
            call YNOGK(pcros-0.0001,f1234,lambda,ql,0.d0,1.d0,a,h,1.d0,&
                 r1,mu1,phi1,t1,sigma1, t_r1, t_r2)
            call YNOGK(pcros,f1234,lambda,ql,0.d0,1.d0,a,h,1.d0,&
                              r2,mu2,phi2,t2,sigma2, t_r1, t_r2) 

            if (r2.lt.r1)then
                  kr=-kr
            endif
            dinang = -(sqrt(ql) / r) / (-1.d0 * vt_ + kr * vr_)
      endif
      endif
      return
      end
!-----------------------------------------------------------------------

      
!-----------------------------------------------------------------------
      function demang(a,gdo,mu0,r,alpha,beta)
! Calculates emission angle for the disk mid-plane
        implicit none       
        double precision, intent(in) :: a, gdo, mu0, r, alpha, beta
        double precision             :: mue
        double precision             :: demang
        mue = gdo / r * sqrt( beta**2 + mu0**2 * (alpha**2-a**2) )
        demang = min( mue , 1.d0 )
        ! g = - (1+z)^-1 / [ pe_\mu U^\mu ]
        ! \mue = - pe_\alpha n^\alpha / [ pe_\mu U^\mu ], therefore
        ! \mue = pe_\alpha n^\alpha g|_{z=0}
        return
      end function demang
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
      function arg( z )
      implicit none
      complex z
      real arg
      arg = atan2( aimag(z) , real(z) )
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function darg( z )
      implicit none
      complex (kind=8) z
      double precision darg
      darg = atan2( aimag(z) , real(z) )
      return
      end
!-----------------------------------------------------------------------
