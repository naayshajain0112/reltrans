module rest_frame_mod
   implicit none

contains

   subroutine rest_frame(ear,ne,Gamma,Afe,logne,Cutoff,logxi,thetae,Cp,photar)
   !
   !  Cp : chooses reflection model
   !      -1 xillver      1e15 density and powerlaw illumination  
   !       1 xillverD     high density and powerlaw illumination 
   !       2 xillverDCp   high density and nthcomp  illumination
   !       0 reflionxDCp  reflionx high density and nthcomp  illumination
   !
   !       Last change: Gullo - 2022 Oct

      integer, intent(in) :: ne, Cp
      real, intent(in)    :: ear(0:ne), Gamma, Afe, logne, Cutoff, logxi, thetae
      real, intent(out)   :: photar(ne)
      integer, parameter  :: dim = 6, dimCp = 8
      real                :: xillpar(dim), xillparDCp(dimCp)
      
      if( Cp .ne. 0 ) then
         ! Fill parameter arrays
         xillpar(1) = Gamma     
         xillpar(2) = Afe       
         xillpar(3) = logxi     
         xillpar(4) = Cutoff      
         if( Cp .eq. 1 ) then
            xillpar(4) = logne 
         end if
         xillpar(5) = thetae    
         xillpar(6) = 0.0     
         xillparDCp(1) = Gamma  
         xillparDCp(2) = Afe    
         xillparDCp(3) = logxi  
         xillparDCp(4) = Cutoff   
         xillparDCp(5) = logne 
         xillparDCp(6) = thetae 
         xillparDCp(7) = 0.0  

         call get_xillver(ear, ne, dim, dimCp, xillpar, xillparDCp, Cp, photar)

      else
         ! The model is reflionx
         call normreflionx(ear,ne,Gamma,Afe,logne,Cutoff,logxi,thetae,photar)
      end if

   end subroutine rest_frame

end module rest_frame_mod