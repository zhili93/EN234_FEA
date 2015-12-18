!     Subroutines for basic 3D shell elements



!==========================SUBROUTINE el_linelast_3dbeam ==============================
subroutine el_linelast_3Dbeam(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only:  dNvdx => vol_avg_shape_function_derivatives_3D            !average volume derivative
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,np_nodes,i,j

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,24),B_equal(6,20)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu,h, D44, D11, D12              ! Material properties
    real (prec)  ::  x_slave(3,8), T(24,20), theta(2,4)
    real (prec)  ::  e1(3,1),e2(3,1),e3(3,1),R(6,6),a1(3,1),b1(3,1)
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))


    np_nodes =8
    n_points = 8


    call initialize_integration_points(n_points, np_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    h= element_properties(1)
    E = element_properties(2)
    xnu = element_properties(3)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:2,1:2) = d12
    D(1,1) = d11
    D(2,2) = d11
    !D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

   !===================calculate coords of slave nodes and T matrix===========

       theta(1,1:n_nodes)=dof_total(4:5*n_nodes-1:5)+dof_increment(4:5*n_nodes-1:5)
       theta(2,1:n_nodes)=dof_total(5:5*n_nodes:5)+dof_increment(5:5*n_nodes:5)
      x_slave = 0.d0
      x_slave(1,4)=x(1,1)-0.5d0*h*cos(theta(1,1))*sin(theta(2,1))
      x_slave(2,4)=x(2,1)+0.5d0*h*sin(theta(1,1))
      x_slave(3,4)=x(3,1)+0.5d0*h*cos(theta(1,1))*cos(theta(2,1))
      x_slave(1,1)=x(1,1)+0.5d0*h*cos(theta(1,1))*sin(theta(2,1))
      x_slave(2,1)=x(2,1)-0.5d0*h*sin(theta(1,1))
      x_slave(3,1)=x(3,1)-0.5d0*h*cos(theta(1,1))*cos(theta(2,1))

      x_slave(1,3)=x(1,2)-0.5d0*h*cos(theta(1,2))*sin(theta(2,2))
      x_slave(2,3)=x(2,2)+0.5d0*h*sin(theta(1,2))
      x_slave(3,3)=x(3,2)+0.5d0*h*cos(theta(1,2))*cos(theta(2,2))
      x_slave(1,2)=x(1,2)+0.5d0*h*cos(theta(1,2))*sin(theta(2,2))
      x_slave(2,2)=x(2,2)-0.5d0*h*sin(theta(1,2))
      x_slave(3,2)=x(3,2)-0.5d0*h*cos(theta(1,2))*cos(theta(2,2))

      x_slave(1,8)=x(1,4)-0.5d0*h*cos(theta(1,4))*sin(theta(2,4))
      x_slave(2,8)=x(2,4)+0.5d0*h*sin(theta(1,4))
      x_slave(3,8)=x(3,4)+0.5d0*h*cos(theta(1,4))*cos(theta(2,4))
      x_slave(1,5)=x(1,4)+0.5d0*h*cos(theta(1,4))*sin(theta(2,4))
      x_slave(2,5)=x(2,4)-0.5d0*h*sin(theta(1,4))
      x_slave(3,5)=x(3,4)-0.5d0*h*cos(theta(1,4))*cos(theta(2,4))

      x_slave(1,7)=x(1,3)-0.5d0*h*cos(theta(1,3))*sin(theta(2,3))
      x_slave(2,7)=x(2,3)+0.5d0*h*sin(theta(1,3))
      x_slave(3,7)=x(3,3)+0.5d0*h*cos(theta(1,3))*cos(theta(2,3))
      x_slave(1,6)=x(1,3)+0.5d0*h*cos(theta(1,3))*sin(theta(2,3))
      x_slave(2,6)=x(2,3)-0.5d0*h*sin(theta(1,3))
      x_slave(3,6)=x(3,3)-0.5d0*h*cos(theta(1,3))*cos(theta(2,3))
      !  write(IOW,*) 'dof number'
      !  write(IOW,*)  length_dof_array
      !  write(IOW,*) 'master'
      !!  write(IOW,*)  x
      !  write(IOW,*) 'slave'
      !  write(IOW,*)  x_slave

      T=0.d0

      T(1,1)=1
      T(1,5)=x_slave(3,1)-x(3,1)
      T(2,2)=1
      T(2,4)=-(x_slave(3,1)-x(3,1))
      T(3,3)=1
      T(3,4)=x_slave(2,1)-x(2,1)
      T(3,5)=-(x_slave(1,1)-x(1,1))
      T(10,1)=1
      T(10,5)=x_slave(3,4)-x(3,1)
      T(11,2)=1
      T(11,4)=-(x_slave(3,4)-x(3,1))
      T(12,3)=1
      T(12,4)=x_slave(2,4)-x(2,1)
      T(12,5)=-(x_slave(1,4)-x(1,1))

      T(4,6)=1
      T(4,10)=x_slave(3,2)-x(3,2)
      T(5,7)=1
      T(5,9)=-(x_slave(3,2)-x(3,2))
      T(6,8)=1
      T(6,9)=x_slave(2,2)-x(2,2)
      T(6,10)=-(x_slave(1,2)-x(1,2))
      T(7,6)=1
      T(7,10)=x_slave(3,3)-x(3,2)
      T(8,7)=1
      T(8,9)=-(x_slave(3,3)-x(3,2))
      T(9,8)=1
      T(9,9)=x_slave(2,3)-x(2,2)
      T(9,10)=-(x_slave(1,3)-x(1,2))

      T(13,16)=1
      T(13,20)=x_slave(3,5)-x(3,4)
      T(14,17)=1
      T(14,19)=-(x_slave(3,5)-x(3,4))
      T(15,18)=1
      T(15,19)=x_slave(2,5)-x(2,4)
      T(15,20)=-(x_slave(1,5)-x(1,4))
      T(22,16)=1
      T(22,20)=x_slave(3,8)-x(3,4)
      T(23,17)=1
      T(23,19)=-(x_slave(3,8)-x(3,4))
      T(24,18)=1
      T(24,19)=x_slave(2,8)-x(2,4)
      T(24,20)=-(x_slave(1,8)-x(1,4))

      T(16,11)=1
      T(16,15)=x_slave(3,6)-x(3,3)
      T(17,12)=1
      T(17,14)=-(x_slave(3,6)-x(3,3))
      T(18,13)=1
      T(18,14)=x_slave(2,6)-x(2,3)
      T(18,15)=-(x_slave(1,6)-x(1,3))
      T(19,11)=1
      T(19,15)=x_slave(3,7)-x(3,3)
      T(20,12)=1
      T(20,14)=-(x_slave(3,7)-x(3,3))
      T(21,13)=1
      T(21,14)=x_slave(2,7)-x(2,3)
      T(21,15)=-(x_slave(1,7)-x(1,3))

  !================== calculate stiffness matix====================


    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),np_nodes,N,dNdxi)
        dxdxi = matmul(x_slave(1:3,1:np_nodes),dNdxi(1:np_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:np_nodes,1:3) = matmul(dNdxi(1:np_nodes,1:3),dxidx)

        B = 0.d0
        B(1,1:3*np_nodes-2:3) = dNdx(1:np_nodes,1)
        B(2,2:3*np_nodes-1:3) = dNdx(1:np_nodes,2)
        B(3,3:3*np_nodes:3)   = dNdx(1:np_nodes,3)
        B(4,1:3*np_nodes-2:3) = dNdx(1:np_nodes,2)
        B(4,2:3*np_nodes-1:3) = dNdx(1:np_nodes,1)
        B(5,1:3*np_nodes-2:3) = dNdx(1:np_nodes,3)
        B(5,3:3*np_nodes:3)   = dNdx(1:np_nodes,1)
        B(6,2:3*np_nodes-1:3) = dNdx(1:np_nodes,3)
        B(6,3:3*np_nodes:3)   = dNdx(1:np_nodes,2)

        !=========calculate rotation matrix R===========

        e1(1:3,1)=matmul(x_slave(1:3,1:8),dNdxi(1:8,3))
        e2(1:3,1)=matmul(x_slave(1:3,1:8),dNdxi(1:8,1))
        e3(1,1)=e1(2,1)*e2(3,1)-e1(3,1)*e2(2,1)
        e3(2,1)=-(e1(1,1)*e2(3,1)-e1(3,1)*e2(1,1))
        e3(3,1)=e1(1,1)*e2(2,1)-e1(2,1)*e2(1,1)
        e3=e3/sqrt(e3(1,1)**2+e3(2,1)**2+e3(3,1)**2)

        a1=(e1/sqrt(e1(1,1)**2+e1(2,1)**2+e1(3,1)**2))+(e2/sqrt(e2(1,1)**2+e2(2,1)**2+e2(3,1)**2))
        b1(1,1)=e3(2,1)*a1(3,1)-e3(3,1)*a1(2,1)
        b1(2,1)=-(e3(1,1)*a1(3,1)-e3(3,1)*a1(1,1))
        b1(3,1)=e3(1,1)*a1(2,1)-e3(2,1)*a1(1,1)

        e1=a1-b1
        e1=e1/sqrt(e1(1,1)**2+e1(2,1)**2+e1(3,1)**2)
        e2=a1+b1
        e2=e2/sqrt(e2(1,1)**2+e2(2,1)**2+e2(3,1)**2)

        R(1,1)=e1(1,1)**2
        R(1,2)=e1(2,1)**2
        R(1,3)=e1(3,1)**2
        R(1,4)=e1(1,1)*e1(2,1)
        R(1,5)=e1(1,1)*e1(3,1)
        R(1,6)=e1(2,1)*e1(3,1)
        R(2,1)=e2(1,1)**2
        R(2,2)=e2(2,1)**2
        R(2,3)=e2(3,1)**2
        R(2,4)=e2(1,1)*e2(2,1)
        R(2,5)=e2(1,1)*e2(3,1)
        R(2,6)=e2(2,1)*e2(3,1)
        R(3,1)=e3(1,1)**2
        R(3,2)=e3(2,1)**2
        R(3,3)=e3(3,1)**2
        R(3,4)=e3(1,1)*e3(2,1)
        R(3,5)=e3(1,1)*e3(3,1)
        R(3,6)=e3(2,1)*e3(3,1)
        R(4,1)=2.d0*e1(1,1)*e2(1,1)
        R(4,2)=2.d0*e1(2,1)*e2(2,1)
        R(4,3)=2.d0*e1(3,1)*e2(3,1)
        R(4,4)=e1(1,1)*e2(2,1)+e1(2,1)*e2(1,1)
        R(4,5)=e1(1,1)*e2(3,1)+e1(3,1)*e2(1,1)
        R(4,6)=e1(2,1)*e2(3,1)+e1(3,1)*e2(2,1)
        R(5,1)=2.d0*e1(1,1)*e3(1,1)
        R(5,2)=2.d0*e1(2,1)*e3(2,1)
        R(5,3)=2.d0*e1(3,1)*e3(3,1)
        R(5,4)=e1(1,1)*e3(2,1)+e1(2,1)*e3(1,1)
        R(5,5)=e1(1,1)*e3(3,1)+e1(3,1)*e3(1,1)
        R(5,6)=e1(2,1)*e3(3,1)+e1(3,1)*e3(2,1)
        R(6,1)=2.d0*e2(1,1)*e3(1,1)
        R(6,2)=2.d0*e2(2,1)*e3(2,1)
        R(6,3)=2.d0*e2(3,1)*e3(3,1)
        R(6,4)=e2(1,1)*e3(2,1)+e2(2,1)*e3(1,1)
        R(6,5)=e2(1,1)*e3(3,1)+e2(3,1)*e3(1,1)
        R(6,6)=e2(2,1)*e3(3,1)+e2(3,1)*e3(2,1)
       ! write(IOW,*) 'rotation'
        !write(IOW,*)  R


        strain = matmul(B,matmul(T,dof_total))
        dstrain = matmul(B,matmul(T,dof_increment))

        stress = matmul(D,matmul(R,strain+dstrain))
        B_equal=matmul(R,matmul(B,T))
        !write(IOW,*) 'stress'
        !write(IOW,*)  stress
        element_residual(1:5*n_nodes) = element_residual(1:5*n_nodes) - &
        matmul(transpose(B_equal),stress)*determinant*w(kint)

        element_stiffness(1:5*n_nodes,1:5*n_nodes) = element_stiffness(1:5*n_nodes,1:5*n_nodes) &
            + matmul(transpose(B_equal),matmul(D,B_equal))*w(kint)*determinant

    end do

    return
end subroutine el_linelast_3Dbeam

