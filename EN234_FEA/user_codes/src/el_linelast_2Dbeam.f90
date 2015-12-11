!     Subroutines for basic 2d beam elements



!==========================SUBROUTINE el_linelast_2dbeam ==============================
subroutine el_linelast_2dbeam(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only:  dNvdx => vol_avg_shape_function_derivatives_2D             ! average derivatives
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
    integer      :: n_points,kint,np_nodes

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(2)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(2,2)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,8),B_equal(2,6)           ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  h,k,E, xnu            ! Material properties
    real (prec)  ::  T(8,6),theta1,theta2,x_slave(2,4)
    real (prec)  ::  e1(2,1),e2(2,1),R(2,3)

    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/2,length_coord_array/2/))
   ! four nodes parenent element
    np_nodes = 4
    n_points = 4

   ! write(IOW,*) 'coords'
   ! write(IOW,*)  x

    call initialize_integration_points(n_points, np_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    h=element_properties(1)
    k=element_properties(2)
    E = element_properties(3)
    xnu = element_properties(4)

    D(1,1) = E
    D(2,2) = 0.5d0*E/(1.d0+xnu)

    !write(IOW,*) 'material'
    !write(IOW,*)  D
   !==========calculate transform matrix
    x_slave=0.d0
    T=0.d0
    theta1=dof_total(3)+dof_increment(3)
    theta2=dof_total(6)+dof_increment(6)
    ! theta1=0.d0
    !  theta2=0.d0
    write(IOW,*) 'displacement'
    write(IOW,*)  dof_total
    write(IOW,*)  dof_increment
    x_slave(1,1)=x(1,1)+0.5*h*sin(theta1)
    x_slave(2,1)=x(2,1)-0.5*h*cos(theta1)
    x_slave(1,4)=x(1,1)-0.5*h*sin(theta1)
    x_slave(2,4)=x(2,1)+0.5*h*cos(theta1)
    x_slave(1,2)=x(1,2)+0.5*h*sin(theta2)
    x_slave(2,2)=x(2,2)-0.5*h*cos(theta2)
    x_slave(1,3)=x(1,2)-0.5*h*sin(theta2)
    x_slave(2,3)=x(2,2)+0.5*h*cos(theta2)

   ! do kint=1,4
   !    if (kint<3)then
   !    T((kint-1)*2+1,1)=1.d0
    !   T((kint-1)*2+1,3)=x(2,1)-x_slave(2,kint)
    !   T((kint-1)*2+2,2)=1.d0
    !   T((kint-1)*2+2,3)=-x(1,1)+x_slave(1,kint)
    !   else
    !   T((kint-1)*2+1,4)=1.d0
    !   T((kint-1)*2+1,6)=x(2,2)-x_slave(2,kint)
    !   T((kint-1)*2+2,5)=1.d0
    !   T((kint-1)*2+2,6)=-x(1,2)+x_slave(1,kint)
    !   end if
   ! end do
   T(1,1)=1.d0
   T(1,3)=x(2,1)-x_slave(2,1)
   T(2,2)=1.d0
   T(2,3)=-x(1,1)+x_slave(1,1)

   T(7,1)=1.d0
   T(7,3)=x(2,1)-x_slave(2,4)
   T(8,2)=1.d0
   T(8,3)=-x(1,1)+x_slave(1,4)

   T(3,4)=1.d0
   T(3,6)=x(2,2)-x_slave(2,2)
   T(4,5)=1.d0
   T(4,6)=-x(1,2)+x_slave(1,2)

   T(5,4)=1.d0
   T(5,6)=x(2,2)-x_slave(2,3)
   T(6,5)=1.d0
   T(6,6)=-x(1,2)+x_slave(1,3)
   ! write(IOW,*) 'slave'
   ! write(IOW,*)  x_slave
   ! write(IOW,*) 'transform'
   ! write(IOW,*)  T
   !write(IOW,*) 'coord'
   !write(IOW,*)  dof_total+dof_increment
   !================== calculate stiffness matix====================

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),np_nodes,N,dNdxi)
        dxdxi = matmul(x_slave(1:2,1:np_nodes),dNdxi(1:np_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:np_nodes,1:2) = matmul(dNdxi(1:np_nodes,1:2),dxidx)
       !==========calculate B matrix for parent element============
        B = 0.d0
        B(1,1:2*np_nodes-1:2) = dNdx(1:np_nodes,1)
        B(2,2:2*np_nodes:2) = dNdx(1:np_nodes,2)
        B(3,1:2*np_nodes-1:2) = dNdx(1:np_nodes,2)
        B(3,2:2*np_nodes:2) = dNdx(1:np_nodes,1)

        !=========calculate rotation matrix R===========

        e1(1:2,1)=matmul(x_slave(1:2,1:4),dNdxi(1:4,1))
        e1=e1/(sqrt(e1(1,1)**2+e1(2,1)**2))
        e2(1,1)=-e1(2,1)
        e2(2,1)=e1(1,1)
        R=0.d0
        R(1,1)=e1(1,1)**2
        R(1,2)=e2(2,1)**2
        R(1,3)=e1(1,1)*e1(2,1)
        R(2,1)=e1(1,1)*e2(1,1)
        R(2,2)=e1(2,1)*e2(2,1)
        R(2,3)=0.5d0*(e1(1,1)*e2(2,1)+e1(2,1)*e2(1,1))

        !==========calculate stiffness and residual======

        strain = matmul(B,matmul(T,dof_total))
        dstrain = matmul(B,matmul(T,dof_increment))
      
        stress = matmul(D,matmul(R,strain+dstrain))
        B_equal=matmul(R,matmul(B,T))
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - &
        matmul(transpose(B_equal),stress)*w(kint)*determinant*k

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B_equal),matmul(D,B_equal))*w(kint)*determinant*k
    !write(IOW,*) 'stress'
    !write(IOW,*)  stress
    !write(IOW,*) 'e1'
    !write(IOW,*)  e1
    !write(IOW,*) 'rotation'
    !write(IOW,*)  R

    end do
   ! write(IOW,*) 'element_residual'
   ! write(IOW,*)  element_residual
    !write(IOW,*) 'stiffness'
   ! write(IOW,*)  element_stiffness
   ! write(IOW,*) 'dof'
   ! write(IOW,*)  dof_total+dof_increment
  
    return
end subroutine el_linelast_2dbeam


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_2dbeam_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
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
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )

    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_linelast_2dbeam_dynamic


!==========================SUBROUTINE fieldvars_linelast_2dbasic ==============================
subroutine fieldvars_linelast_2dbeam(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only:  dNvdx => vol_avg_shape_function_derivatives_2D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

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

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,np_nodes

    real (prec)  ::  sdev(2)                           ! Deviatoric stress
    real (prec)  ::  p, smises                          ! Pressure and Mises stress
    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(2)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(2,2)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,length_dof_array),B_equal(2,6)           ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  h,k1,E, xnu            ! Material properties
    real (prec)  ::  T(8,6),theta1,theta2,x_slave(2,4)
    real (prec)  ::  e1(2,1),e2(2,1),R(2,3)
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/2,length_coord_array/2/))

   ! four nodes parenent element
    np_nodes = 4
    n_points = 4

    call initialize_integration_points(n_points, np_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    h=element_properties(1)
    k1=element_properties(2)
    E = element_properties(3)
    xnu = element_properties(4)

    D(1,1) = E
    D(2,2) = 0.5d0*E/(1.d0+xnu)

 !==========calculate transform matrix
    x_slave=0.d0
    T=0.d0
    theta1=dof_total(3)+dof_increment(3)
    theta2=dof_total(6)+dof_increment(6)
    x_slave(1,1)=x(1,1)+0.5*h*sin(theta1)
    x_slave(2,1)=x(2,1)-0.5*h*cos(theta1)
    x_slave(1,2)=x(1,1)-0.5*h*sin(theta1)
    x_slave(2,2)=x(2,1)+0.5*h*cos(theta1)
    x_slave(1,3)=x(1,2)-0.5*h*sin(theta2)
    x_slave(2,3)=x(2,2)+0.5*h*cos(theta2)
    x_slave(1,4)=x(1,2)+0.5*h*sin(theta2)
    x_slave(2,4)=x(2,2)-0.5*h*cos(theta2)

    do kint=1,4
       if (kint<3)then
       T((kint-1)*2+1,1)=1.d0
       T((kint-1)*2+1,3)=x(2,1)-x_slave(2,kint)
       T((kint-1)*2+2,2)=1.d0
       T((kint-1)*2+2,3)=-x(1,1)+x_slave(1,kint)
       else
       T((kint-1)*2+1,4)=1.d0
       T((kint-1)*2+1,6)=x(2,2)-x_slave(2,kint)
       T((kint-1)*2+2,5)=1.d0
       T((kint-1)*2+2,6)=-x(1,2)+x_slave(1,kint)
       end if
    end do

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),np_nodes,N,dNdxi)
        dxdxi = matmul(x_slave(1:2,1:np_nodes),dNdxi(1:np_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:np_nodes,1:2) = matmul(dNdxi(1:np_nodes,1:2),dxidx)

       !==========calculate B matrix for parent element============
        B = 0.d0
        B(1,1:2*np_nodes-1:2) = dNdx(1:np_nodes,1)
        B(2,2:2*np_nodes:2) = dNdx(1:np_nodes,2)
        B(3,1:2*np_nodes-1:2) = dNdx(1:np_nodes,2)
        B(3,2:2*np_nodes:2) = dNdx(1:np_nodes,1)

        !=========calculate rotation matrix R===========

        e1(1:2,1)=matmul(x_slave(1:2,1:4),dNdxi(1:4,2))
        e1=e1/(sqrt(e1(1,1)**2+e1(2,1)**2))
        e2(1,1)=-e1(2,1)
        e2(2,1)=e1(1,1)
        R=0.d0
        R(1,1)=e1(1,1)**2
        R(1,2)=e2(2,1)**2
        R(1,3)=e1(1,1)*e1(2,1)
        R(2,1)=e1(1,1)*e2(1,1)
        R(2,2)=e1(2,1)*e2(2,1)
        R(2,3)=0.5d0*(e1(1,1)*e2(2,1)+e1(2,1)*e2(1,1))

        !==========calculate stiffness and residual======

        strain = matmul(B,matmul(T,dof_total))
        dstrain = matmul(B,matmul(T,dof_increment))
        stress = matmul(D,matmul(R,strain+dstrain))


       ! p = stress(1)
        !sdev = stress
        !sdev(1:2) = sdev(1:2)-p
        !smises = dsqrt( dot_product(sdev(1:2),sdev(1:2)) + 2.d0*dot_product(sdev(3:3),sdev(3:3)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
           else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_linelast_2dbeam

