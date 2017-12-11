subroutine impedance_driver(tmstps,Period,rho,mu,r_root,r_min,y11,y12,y21,y22,Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV)
  use f90_tools
  use new_match 
  implicit none

  integer, intent(in)    :: tmstps
  real(lng), intent(in)  :: Period,rho,mu,r_root,r_min,Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV
  real(lng), intent(out) :: y11(tmstps),y12(tmstps),y21(tmstps),y22(tmstps)
  integer :: j
!write(*,*)'Impedancedriver started' !mpaun commented out
!write(*,*)'timesteps',tmstps !mpaun
  do j = 1, tmstps
    y11(j) = 0.0
    y12(j) = 0.0
    y21(j) = 0.0
    y22(j) = 0.0
!write(*,*)j,' Impedance initiated'
  end do
!write(*,*)'New match called'
  call impedance (tmstps,Period,rho,mu,r_root,r_min,y11,y12,y21,y22,Lr,q,g,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV)
end subroutine impedance_driver
