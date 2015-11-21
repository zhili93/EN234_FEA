!     Subroutines for 3d explicit dynamic elements


!==========================SUBROUTINE explicit dynamic ==============================
subroutine explicit_dynamics_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only : gurson
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
    integer      :: n_points,kint,point_state_variables,ii,failure_count

    real (prec)  ::  depsilon(3,3),dW(3,3),dR(3,3)
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
!   real (prec)  ::                                   ! Material properties
    real (prec)  ::  I(3,3),vel
    real (prec)  ::  df(3,3),f_mid(3,3),f_mid_i(3,3),J,eta,d_eta,dL(3,3),dNdy(20,3),dNvdy(20,3)
    real (prec)  ::  u(3,length_coord_array/3),du(3,length_coord_array/3)
    real (prec)  ::  dL_bar(3,3), temp(3,3),temp_i(3,3), temp_d
    real (prec)  :: el_initial_state_variables(8),el_updated_state_variables(8)
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
    stress=0.d0

    I=0.d0
    I(1,1)=1.d0
    I(2,2)=1.d0
    I(3,3)=1.d0

    du = reshape(dof_increment,(/3,length_coord_array/3/))
    u = reshape(dof_total,(/3,length_coord_array/3/))

    eta=0.d0
    d_eta=0.d0
    dNvdy=0.d0
    vel=0.d0

!============================calculate volume average items==================================
      do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        df(1:3,1:3)=matmul(du(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        f_mid(1:3,1:3)=matmul((u(1:3,1:n_nodes)+0.5d0*du(1:3,1:n_nodes)),dNdx(1:n_nodes,1:3))+I(1:3,1:3)
        call invert_small(f_mid,f_mid_i,J)

        eta=eta+w(kint)*determinant*J
        dL=matmul(df,f_mid_i)
        d_eta=d_eta+w(kint)*determinant*J*(dL(1,1)+dL(2,2)+dL(3,3))

        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),f_mid_i(1:3,1:3))
        dNvdy(1:n_nodes,1:3)=dNvdy(1:n_nodes,1:3)+dNdy(1:n_nodes,1:3)*w(kint)*determinant*J

        vel=vel+w(kint)*determinant
      end do

      eta=eta/vel
      d_eta=d_eta/(eta*vel)
      dNvdy(1:n_nodes,1:3)=dNvdy(1:n_nodes,1:3)/(eta*vel)

!===============calculate element residual=====================

    point_state_variables = 8
    failure_count=0
    !     --  Loop over integration points
      do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        df(1:3,1:3)=matmul(du(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        f_mid(1:3,1:3)=matmul((u(1:3,1:n_nodes)+0.5d0*du(1:3,1:n_nodes)),dNdx(1:n_nodes,1:3))+I(1:3,1:3)
        call invert_small(f_mid,f_mid_i,J)

        dL=matmul(df,f_mid_i)
        dL_bar=dL+I*((d_eta-dL(1,1)-dL(2,2)-dL(3,3))/3.d0)

        depsilon= (dL_bar+transpose(dL_bar))*0.5d0
        dW=(dL_bar-transpose(dL_bar))*0.5d0
        temp=I-dW*0.5d0
        call invert_small(temp,temp_i,temp_d)
        dR=matmul(temp_i,(I+dW*0.5d0))
        !==============calculate B_bar matrix
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),f_mid_i(1:3,1:3))
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        B(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3)-dNdy(1:n_nodes,1)/3.d0+dNvdy(1:n_nodes,1)/3.d0
        B(1,2:3*n_nodes-1:3) = B(1,2:3*n_nodes-1:3)-dNdy(1:n_nodes,2)/3.d0+dNvdy(1:n_nodes,2)/3.d0
        B(1,3:3*n_nodes:3) =   B(1,3:3*n_nodes:3)-dNdy(1:n_nodes,3)/3.d0+dNvdy(1:n_nodes,3)/3.d0
        B(2,1:3*n_nodes-2:3) = B(2,1:3*n_nodes-2:3)-dNdy(1:n_nodes,1)/3.d0+dNvdy(1:n_nodes,1)/3.d0
        B(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3)-dNdy(1:n_nodes,2)/3.d0+dNvdy(1:n_nodes,2)/3.d0
        B(2,3:3*n_nodes:3) =   B(2,3:3*n_nodes:3)-dNdy(1:n_nodes,3)/3.d0+dNvdy(1:n_nodes,3)/3.d0
        B(3,1:3*n_nodes-2:3) = B(3,1:3*n_nodes-2:3)-dNdy(1:n_nodes,1)/3.d0+dNvdy(1:n_nodes,1)/3.d0
        B(3,2:3*n_nodes-1:3) = B(3,2:3*n_nodes-1:3)-dNdy(1:n_nodes,2)/3.d0+dNvdy(1:n_nodes,2)/3.d0
        B(3,3:3*n_nodes:3) =   B(3,3:3*n_nodes:3)-dNdx(1:n_nodes,3)/3.d0+dNvdy(1:n_nodes,3)/3.d0



        el_initial_state_variables=0.d0
        el_updated_state_variables=0.d0
        ii=(kint-1)*point_state_variables+1
        el_initial_state_variables(1:8)=initial_state_variables(ii:(ii+7))

       ! write(IOW,*) el_initial_state_variables
      call gurson(element_properties,n_properties,point_state_variables,&
        el_initial_state_variables,el_updated_state_variables,depsilon,dR,stress)

        updated_state_variables(ii:(ii+7))=el_updated_state_variables(1:8)


        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) &
                                - matmul(transpose(B),stress)*w(kint)*determinant

        if(el_updated_state_variables(8)>element_properties(13)) then
        failure_count=failure_count+1
        end if
       !write(IOW,*) el_updated_state_variables
       !write(IOW,*) dR
    end do
         write(IOW,*) element_deleted
        if(failure_count==n_points) then
         element_deleted =.true.
         end if

    return

end subroutine explicit_dynamics_3d






!==========================SUBROUTINE fieldvars_explicit_dynamics_3d ==============================
subroutine fieldvars_explicit_dynamics_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
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

    integer      :: n_points,kint,k,ii,point_state_variables

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
 !   real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  I(3,3),vel,vf
    real (prec)  ::  df(3,3),f_mid(3,3),f_mid_i(3,3),J,eta,d_eta,dL(3,3),dNdy(20,3),dNvdy(20,3)
    real (prec)  ::  u(3,length_coord_array/3),du(3,length_coord_array/3)
    real (prec)  ::  dL_bar(3,3), temp(3,3),temp_i(3,3), temp_d
    real (prec)  :: el_initial_state_variables(8),el_updated_state_variables(8)

    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio


    !     --  Loop over integration points


    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)



     nodal_fieldvariables = 0.d0


   do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ii=(kint-1)*8+1
        stress(1:6)=updated_state_variables(ii:(ii+5))
        vf=updated_state_variables(ii+7)

       ! write(IOW,*) stress

       ! p = sum(stress(1:3))/3.d0
       ! sdev = stress
       ! sdev(1:3) = sdev(1:3)-p
       ! smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
           ! else if (strcmp(field_variable_names(k),'SMISES',6) ) then
            !    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + vf*N(1:n_nodes)*determinant*w(kint)

            endif
        end do

    end do
    return
end subroutine fieldvars_explicit_dynamics_3d



