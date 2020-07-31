Program Incompressible_vorticity
!
! This is a Serial Benchmarked Two Dimensional Incompressible Viscous Fluid code, 
! using Pseudo-Spectral Method for Spatial Discritisation and Adams-Bashforth Technique for time evolution.
!
! Wave-form /sim exp(i(kx-wt)). Hence derivative w.r.t. x gives ik only (and NOT -ik).
!
!___________________________________________________________________________________________________________________________________________
!
! Equations Solved in this code are in Vorticity (\w) - Stream Function (\psi) Formalism.
!
!\begin{eqnarray}
!&& \frac{\partial \vec{\w}}{\partial t} + \vec{u} \cdot \vec{\nabla} \vec{\w} = \nu \nabla^2 \vec{\w} \\
!&& \vec{\w} = \vec{\nabla} \times \vec{u} \\
!&& \nabla^2 \psi = - \w 
!\end{eqnarray}
!
!___________________________________________________________________________________________________________________________________________
!
implicit none

! Define Grid Size.
!integer ( kind = 4 ), parameter :: Nx = 256
!integer ( kind = 4 ), parameter :: Ny = 256
!integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi = 3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) i,j,t,ifile_write,itime,Ndt,file_no,Nx,Ny,Nh
real ( kind = 8 ) st_time,end_time
real ( kind = 8 ) xmin,xmax,ymin,ymax,Lx,Ly,dx,dy,time,time_min,time_max,dt,Rex,Rey,nu,W0,Az0,m,d,iNxNy, Pi_multpl, pert

real ( kind = 8 ) KE,KEx,KEy,ME,MEx,MEy

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,y,kx,ky
DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: psi,w,w_dum,ux,uy,dw_dx,dw_dy,dpsi_dx,dpsi_dy,rhs_w_old

complex(8), DIMENSION(:,:), ALLOCATABLE ::  wk,wk_dum,wk_new,psik,psik_dum,ukx,uky,ukx_dum,uky_dum
complex(8), DIMENSION(:,:), ALLOCATABLE ::  ikx_psik,iky_psik,ikx_wk,iky_wk,rhs_wk_old,rhs_wk_older,dissip

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW Variables.

namelist /time_param/ time_max, time_min, dt, itime, file_no, ifile_write
namelist /mesh/  Nx,Ny,Pi_multpl,xmin, xmax, ymin, ymax
namelist /phy/ W0,m,d,pert,nu

!===================== USER INPUTS ============================================		
  call CPU_TIME(st_time)
  open(unit=2, file='code.settings', status='unknown')
  open(unit=1, file='input.inp', form='formatted', status='old')
  rewind(1)

  read(1, time_param)
  read(1, mesh)
  read(1, phy)

  close(1)
  
  iNxNy = 1.0d0/(dfloat(Nx*Ny))
  Nh = ( Nx / 2 ) + 1

! System Length.
  if (Pi_multpl == 1)then
    xmin=xmin*pi; xmax=xmax*pi
    ymin=ymin*pi; ymax=ymax*pi
  endif
!  ymin=-pi; ymax=pi
  Lx = xmax-xmin
  Ly = ymax-ymin

  dx = Lx/dfloat(Nx)
  dy = Ly/dfloat(Ny)

  Ndt = (time_max - time_min)/dt
  write(*,*)'Total no of time step:',Ndt
  write(*,*)'file writing frequency',ifile_write,'with every',ifile_write*dt,'time'
  
  write(2,*)'Mesh information:'
  write(2,*)'Nx:',Nx
  write(2,*)'Ny:',Ny
  write(2,*)'Nh',Nh
  write(2,*)'dx:',dx
  write(2,*)'dy:',dy
  write(2,*)'Lx:',Lx
  write(2,*)'xmin:',xmin
  write(2,*)'xmax:',xmax
  write(2,*)'Ly:',Ly
  write(2,*)'ymin:',ymin
  write(2,*)'ymax:',ymax
  write(2,*)'time stepping information:'
  write(2,*)'Total time units:',time_max-time_min,',time step:',dt,',Total time step:',Ndt,',file writing freq:',ifile_write
  close(2)
  
  ALLOCATE ( x(Nx),y(Ny),kx(Nh),ky(Ny) )
  ALLOCATE ( psi(Nx,Ny),w(Nx,Ny),w_dum(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny) )
  ALLOCATE ( dw_dx(Nx,Ny),dw_dy(Nx,Ny),dpsi_dx(Nx,Ny),dpsi_dy(Nx,Ny),rhs_w_old(Nx,Ny) )

  ALLOCATE ( wk(Nh,Ny),wk_dum(Nh,Ny),wk_new(Nh,Ny),psik(Nh,Ny),psik_dum(Nh,Ny) )
  ALLOCATE ( ukx(Nh,Ny),uky(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny),ikx_psik(Nh,Ny),iky_psik(Nh,Ny) )
  ALLOCATE ( ikx_wk(Nh,Ny),iky_wk(Nh,Ny),rhs_wk_old(Nh,Ny),rhs_wk_older(Nh,Ny),dissip(Nh,Ny) )

  call grid(Nx,Ny,xmin,xmax,ymin,ymax,dx,dy,x,y,kx,ky)

! input Generation.
  do i = 1, Nx
    do j = 1, Ny
      ux(i,j) = 0.0d0
      uy(i,j) = 0.0d0
    ! Initial Vorticity Profile.
      if( (dabs(y(j)-1.5d0*pi) .le. d) )then
        w(i,j) = 1.0d0 !W0 /dcosh((y(j)-1.50d0*pi)/d)**2.0-W0/dcosh((y(j)-0.50d0*pi)/d)**2.0
      else if ( (dabs(y(j)-0.5d0*pi) .le. d) )then
        w(i,j) = -1.0d0
      endif
!      w(i,j) = W0/dcosh((y(j)-1.50d0*pi)/d)**2.0 - W0/dcosh((y(j)-0.50d0*pi)/d)**2.0
      w(i,j) = w(i,j)+pert*dcos(m*x(i))
    ! Keep backup for FFTW.
      w_dum(i,j) = w(i,j)
    end do ! j
  end do ! i


!===================== INITIAL TIME DATA ===============================================

! Move to Spectral Space.
  call real2spectral(Nx,Ny,Nh,w_dum,wk)

  do j = 1,Ny
    do i = 1,Nh!/2
      if (i == 1 .and. j == 1) then
        psik(i,j) = (0.0d0,0.0d0)
      else
        psik(i,j) = wk(i,j)/( kx(i)**2 + ky(j)**2 ) 
      endif
      iky_psik(i,j) = + (0.0d0,1.0d0) * ky(j) * psik(i,j) 
      ikx_psik(i,j) = + (0.0d0,1.0d0) * kx(i) * psik(i,j)
      ukx(i,j) = + iky_psik(i,j)
      uky(i,j) = - ikx_psik(i,j)

      wk_dum(i,j)=wk(i,j)
      psik_dum(i,j)=psik(i,j)
    enddo
  enddo
  
! Move to Real Space.
  call spectral2real(Nx,Ny,Nh,psik,psi)
  call spectral2real(Nx,Ny,Nh,ukx,ux)
  call spectral2real(Nx,Ny,Nh,uky,uy)
  
  do i = 1,Nx
    do j = 1,Ny
      psi(i,j) = psi(i,j)*iNxNy
      ux(i,j) = ux(i,j)*iNxNy
      uy(i,j) = uy(i,j)*iNxNy
    end do ! j
  end do ! i
  
  write(1,*)'#time   total_KE   KEy'
9 FORMAT( 12(2x,f22.16) )
!======================= MAIN PROGRAM =====================================================
!$acc data copyin(w,psi,ux,uy,wk), create(ikx_wk,iky_wk,ikx_psik,iky_psik) &
!$acc& create(dw_dx,dw_dy,dpsi_dx,dpsi_dy,rhs_w_old,rhs_wk_old,rhs_wk_older) &
!$acc& create(wk_new,wk_dum,psik_dum,KEx,KEy)
file_no=0
do itime = 0,Ndt
  time=dfloat(itime)*dt
  
  KEx = 0.0d0 ; MEx=0.0d0
  KEy = 0.0d0 ; MEy=0.0d0
  !$acc update device(KEx,KEy)
!  KE = 0.0d0

!$acc parallel loop collapse(2), reduction(+:KEx,KEy)
! energy calculation.   
  do i = 1,Nx
    do j = 1,Ny
      KEx = KEx + (ux(i,j)**2)*0.5d0
      KEy = KEy + (uy(i,j)**2)*0.5d0
    end do ! j
  end do ! i
!$acc end parallel 

!$acc update host(KEx,KEy)
  KEx=KEx*iNxNy
  KEy=KEy*iNxNy
!  KE=KEx+KEy
  
  write(1,9) time,KEx,KEy,MEx,MEy

  if (mod(itime,ifile_write) == 0.0) then
  !$acc update host(ux,uy,psi,w)
    write(10+file_no,*)'#x  y  ux  uy  strm_funcn   w   at time=',itime
     do i=1,Nx
       do j=1,Ny
        write(10+file_no,9) x(i),y(j),ux(i,j),uy(i,j),psi(i,j),w(i,j)
       enddo
       write(10+file_no,*)
     enddo
     close (10+file_no)
     call CPU_TIME(end_time)
     write(*,*)'time step :  ',itime,'  file name :fort.',10+file_no,'time elapsed:',end_time-st_time,'sec.'
     file_no=file_no+1
  endif

!$acc parallel loop collapse(2)
  do j = 1,Ny
    do i = 1,Nh!/2
      iky_wk(i,j) = + (0.0d0,1.0d0) * ky(j) * wk(i,j) 
      ikx_wk(i,j) = + (0.0d0,1.0d0) * kx(i) * wk(i,j)
    enddo
  enddo
!$acc end parallel  
  
! Move to Real Space.
  call spectral2real(Nx,Ny,Nh,ikx_wk,dw_dx)
  call spectral2real(Nx,Ny,Nh,iky_wk,dw_dy)
  
! FFTW Normalisations and Evaluation of Non-Linear Terms.  
!$acc parallel loop collapse(2)
  do j = 1,Ny
    do i = 1,Nx
      dpsi_dx(i,j) = -uy(i,j)
      dpsi_dy(i,j) =  ux(i,j)
      dw_dx(i,j)   =  dw_dx(i,j)*iNxNy
      dw_dy(i,j)   =  dw_dy(i,j)*iNxNy
      rhs_w_old(i,j)  =  (dpsi_dx(i,j)*dw_dy(i,j) - dpsi_dy(i,j)*dw_dx(i,j)) !&

!      rhs_w_old(i,j)  = - (ux(i,j)*dw_dx(i,j) + uy(i,j)*dw_dy(i,j)) 
    end do ! j
  end do ! i
!$acc end parallel

! Move to Spectral Space.
  call real2spectral(Nx,Ny,Nh,rhs_w_old,rhs_wk_old)

!$acc parallel loop collapse(2)
! Evaluation of Nonlinear Term.  
  do j = 1,Ny
    do i = 1,Nh  
      ! De - Aliazing Technique - 2/3 Truncation...
      if (dsqrt(kx(i)*kx(i) + ky(j)*ky(j)) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then    
        rhs_wk_old(i,j) = 0.0d0
      endif
      !rhs calculation
      rhs_wk_old(i,j) = rhs_wk_old(i,j) - nu*(kx(i)*kx(i) + ky(j)*ky(j))*wk(i,j)
      
      !Adam_Bashforth step
      wk_new(i,j) = wk(i,j) + ( (3.0d0/2.0d0)*rhs_wk_old(i,j) - (1.0d0/2.0d0)*rhs_wk_older(i,j) )*dt
      
      !data copying for next time step
      rhs_wk_older(i,j) = rhs_wk_old(i,j)
      wk(i,j) = wk_new(i,j)
      
      !poison solver, derivatives and field calculation
      if (i == 1 .and. j == 1) then
        psik(i,j) = (0.0d0,0.0d0)
      else
        psik(i,j) = wk(i,j)/( kx(i)**2 + ky(j)**2 ) 
      endif
      iky_psik(i,j) = + (0.0d0,1.0d0) * ky(j) * psik(i,j) 
      ikx_psik(i,j) = + (0.0d0,1.0d0) * kx(i) * psik(i,j)
      ukx(i,j) = + iky_psik(i,j)
      uky(i,j) = - ikx_psik(i,j)

      wk_dum(i,j)=wk(i,j)
      psik_dum(i,j)=psik(i,j)
    enddo
  enddo
 !$acc end parallel 

! Move to Real Space.
  call spectral2real(Nx,Ny,Nh,psik_dum,psi)
  call spectral2real(Nx,Ny,Nh,wk_dum,w)
  call spectral2real(Nx,Ny,Nh,ukx,ux)
  call spectral2real(Nx,Ny,Nh,uky,uy)
  
!$acc parallel loop collapse(2)  
  do i = 1,Nx
    do j = 1,Ny
      w(i,j) = w(i,j)*iNxNy
      psi(i,j) = psi(i,j)*iNxNy
      ux(i,j) = ux(i,j)*iNxNy
      uy(i,j) = uy(i,j)*iNxNy
    end do ! j
  end do ! i
!$acc end parallel
  
end do ! time
!$acc end data
 close(1)
contains

!=======================================================================================
subroutine grid(Nx,Ny,xmin,xmax,ymin,ymax,dx,dy,x,y,kx,ky)
  implicit none
  integer i,j,Nx,Ny
  real (kind=8) Lx,Ly,xmin,xmax,ymin,ymax,dx,dy,x(Nx),y(Ny),kx(Nx),ky(Ny)
  
  Lx=xmax-xmin
  Ly=ymax-ymin
  
  dx = Lx/dfloat(Nx)
  dy = Ly/dfloat(Ny)
  
  do i=1,Nx
    x(i)=xmin+real(i-1)*dx    
  enddo
  do i=1,Nx/2+1
    kx(i) = 2.0d0*pi*dfloat(i-1)/Lx
  enddo
  
  do j=1,Ny/2
    ky(j) = 2.0d0*pi*dfloat(j-1)/Ly
    y(j)=ymin+real(j-1)*dy
  enddo
  do j=Ny/2+1,Ny
    ky(j) = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    y(j)=ymin+real(j-1)*dy
  enddo
  
  return
end subroutine
!=======================================================================================
subroutine real2spectral(Nx,Ny,Nh,w_dum,wk)
  implicit none
  integer Nx,Ny,Nh
  real(kind=8) w_dum(Nx,Ny)
  complex(kind=8) wk(Nh,Ny)
  
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, w_dum, wk, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)
  
  return
end subroutine
!=======================================================================================
subroutine spectral2real(Nx,Ny,Nh,wk,w)
  implicit none
  integer Nx,Ny,Nh
  real(kind=8) w(Nx,Ny)
  complex(kind=8) wk(Nh,Ny)
  
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, wk, w, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
  return
end subroutine
!=======================================================================================

end program  Incompressible_vorticity
