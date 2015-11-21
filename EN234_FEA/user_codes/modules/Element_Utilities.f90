module Element_Utilities

    use Types

    real (prec) :: shape_functions_1D(3)
    real (prec) :: shape_functions_2D(9)
    real (prec) :: shape_functions_3D(20)
 
    real (prec) :: shape_function_derivatives_1D(3,1)
    real (prec) :: shape_function_derivatives_2D(9,2)
    real (prec) :: shape_function_derivatives_3D(20,3)
  
    real (prec) :: shape_function_spatial_derivatives_1D(3,1)
    real (prec) :: shape_function_spatial_derivatives_2D(9,2)
    real (prec) :: shape_function_spatial_derivatives_3D(20,3)

    real (prec) :: vol_avg_shape_function_derivatives_1D(3,1)
    real (prec) :: vol_avg_shape_function_derivatives_2D(9,2)
    real (prec) :: vol_avg_shape_function_derivatives_3D(20,3)
  
    real (prec) :: integrationpoints_1D(6,1)
    real (prec) :: integrationpoints_2D(2,9)
    real (prec) :: integrationpoints_3D(3,64)
  
    real (prec) :: integrationweights_1D(6)
    real (prec) :: integrationweights_2D(9)
    real (prec) :: integrationweights_3D(64)
  
    real (prec) :: jacobian_2D(2,2)
    real (prec) :: jacobian_3D(3,3)
  
  



contains

    subroutine gurson (element_properties,n_properties,n_state_variables,initial_state_variables, &
                   updated_state_variables,dstrain,dRot,stress1)

        use Types
        use ParamIO
        use Globals, only: TIME,DTIME

        implicit none

        integer, intent( in )     :: n_properties
        integer, intent( in )     :: n_state_variables
        real (prec), intent( in ) :: element_properties(n_properties)
        real (prec), intent( in ) :: initial_state_variables(8)
        real (prec), intent( in ) :: dstrain(3,3)
        real (prec), intent( in ) :: dRot(3,3)
        real (prec), intent( out ):: stress1(6)
        real (prec), intent( out ):: updated_state_variables(8)



        integer      :: ncount
        real (prec)  :: stress0(6),ematrix,dematrix,vf
        real (prec)  :: E,nv,Y,epsilon_0,m,q1,q2,q3,f_n,epsilon_n,s_n,f_c,f_f,f_f_bar
        real (prec)  :: f_star
        real (prec)  :: p_star,stress0_m(3,3),estress0_m(3,3),I(3,3),destrain(3,3)
        real (prec)  :: S_star(3,3),phi0,phi,sigma_estar
        real (prec)  :: stress1_m(3,3)
        real (prec)  :: deriva(2,2),incre(2,1),resid(2,1),target(2,1)
        real (prec)  :: d_e,d_epsilonv,f1,f2,df1de,df1depsilonv,df2de,df2depsilonv,dt,sigma_star
        real (prec)  :: deriva_i(2,2),deriva_det

       ! initial state variables
        stress0(1:6)=initial_state_variables(1:6)
        ematrix=initial_state_variables(7)
        vf=initial_state_variables(8)
       ! material properties
        E=element_properties(1)
        nv=element_properties(2)
        Y=element_properties(3)
        epsilon_0=element_properties(4)
        m=element_properties(5)
        q1=element_properties(6)
        q2=element_properties(7)
        q3=element_properties(8)
        f_n=element_properties(9)
        epsilon_n =element_properties(10)
        s_n=element_properties(11)
        f_c=element_properties(12)
        f_f=element_properties(13)

        p_star=sum(stress0(1:3))/3.d0+E*(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/(3.d0*(1.d0-2.d0*nv))

        I=0.d0
        I(1,1)=1.d0
        I(2,2)=1.d0
        I(3,3)=1.d0

        stress0_m=0.d0
        stress0_m(1,1)=stress0(1)
        stress0_m(2,2)=stress0(2)
        stress0_m(3,3)=stress0(3)
        stress0_m(1,2)=stress0(4)
        stress0_m(1,3)=stress0(5)
        stress0_m(2,3)=stress0(6)
        stress0_m(2,1)=stress0_m(1,2)
        stress0_m(3,1)=stress0_m(1,3)
        stress0_m(3,2)=stress0_m(2,3)

        estress0_m=stress0_m-(sum(stress0(1:3))/3.d0)*I
        destrain=dstrain-((dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0)*I
        S_star=(E/(1.d0+nv))*destrain+matmul(dRot,matmul(estress0_m,transpose(dRot)))
       !========update state variables=========
 !===============calculate f_star===============================
        if(vf<f_c) then
        f_star=vf
        else if(vf>f_c.AND.vf<f_f) then
        f_f_bar=(q1+sqrt(q1*q1-q3))/q3
        f_star=f_c+((f_f_bar-f_c)*(vf-f_c))/(f_f-f_c)
        else if(vf>f_f) then
        f_f_bar=(q1+sqrt(q1*q1-q3))/q3
        f_star=f_f_bar
        end if

        sigma_estar=sqrt((3.d0*(S_star(1,1)**2+S_star(2,2)**2+S_star(3,3)**2+ &
        2.d0*(S_star(1,2)**2+S_star(2,3)**2+S_star(1,3)**2)))/2.d0)
        phi0=(sigma_estar/Y)**2+2.d0*q1*f_star*cosh(3.d0*q2*p_star/(2.d0*Y))-(1+q3*f_star*f_star)

        if(phi0<0.0000000000001d0) then
        stress1_m=S_star+p_star*I
        stress1(1)=stress1_m(1,1)
        stress1(2)=stress1_m(2,2)
        stress1(3)=stress1_m(3,3)
        stress1(4)=stress1_m(1,2)
        stress1(5)=stress1_m(1,3)
        stress1(6)=stress1_m(2,3)

       ! write(IOW,*) dematrix, ematrix

        updated_state_variables(1:6)=stress1(1:6)
        updated_state_variables(7)=ematrix
        updated_state_variables(8)=vf
        else
       !  write(IOW,*) 'wrong'
        phi=sqrt((sigma_estar/Y)**2+2.d0*q1*f_star*cosh(3.d0*q2*p_star/(2.d0*Y))-(1+q3*f_star*f_star))
        !==========newton raphson========
        deriva=0.d0
        incre(1,1)=1.d0
        incre(2,1)=1.d0
        resid=0.d0
        target=0.d0
        ncount=1
        dt=DTIME
        sigma_star=sigma_estar

        do while(abs(incre(1,1))+abs(incre(2,1))>0.000000001.AND.ncount<30)

        d_e=target(1,1)
        d_epsilonv=target(2,1)

       f1 = (d_e*(- (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 - &
       (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star +&
        (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) - &
        ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
         3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv +&
          2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
          3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0))/(dt*epsilon_0) - ((((sigma_star - &
          (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2.d0/(Y**2 )- (f_star**2)*q3 + &
          2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) -&
           1.d0)**(0.5d0))**m*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/((Y**2)*(- (f_star**2)*q3 +& !???m blan
            (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
            + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))


       f2 =(d_epsilonv*(- (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 &
        - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
         + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) - &
         ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
         3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
         2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
          3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0))/(dt*epsilon_0) - (3.d0*f_star*q1*q2*sinh((q2*(3.d0*p_star &
          + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*(((sigma_star - &
          (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) &
          - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
          (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) &
          - 1.d0)**(0.5d0))**m)/(2.d0*Y*(- (f_star**2)*q3 + (sigma_star - &
          (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) + &
          2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv &
          - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))

        df1de =(- (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 -&
         (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
         2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
         3.d0)))/(2.d0*Y)) + 1.d0)) - ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + &
         (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 &
         - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
         + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0)/(dt*epsilon_0) &
         + (d_e*((6.d0*E*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/((Y**4)*(2.d0*nv + 2.d0)*((f_star**2)*q3 &
         - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
         2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) &
         + (6.d0*E*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**3)/((Y**6)*(2.d0*nv + 2.d0)*((f_star**2)*q3 &
         - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
         + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)**2) + &
         (3.d0*E*(f_star**2)*(q1**2)*(q2**2)*(sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv &
         - 3.d0)))/(2.d0*Y))**2)*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + &  ! add a blancket? sinh
         2.d0)))/((Y**4)*(2.d0*nv + 2.d0)*((f_star**2)*q3 - (sigma_star - &
         (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
         (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)**2)))/(2.d0*dt*epsilon_0*(- (sigma_star &
         - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 - (sigma_star - &
         (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star +&
          (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) - &
          ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
          3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - (sigma_star - &
          (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
          (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0)) + (3.d0*E*(((sigma_star &
          - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
          + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**m)/((Y**2)*(2.d0*nv +&
           2.d0)*((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - (f_star**2)*q3 + &
           2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - &
           1.d0)**(0.5d0)) - (3.d0*E*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
           (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
           3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**m*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
           2.d0))**2)/((Y**4)*(2.d0*nv + 2.d0)*((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) -&
            (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
             3.d0)))/(2.d0*Y)) - 1.d0)**(1.5d0)) - (3.d0*E*m*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
             2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
             (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**(m - 1.d0)*(sigma_star &
             - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2)/((Y**4)*(2.d0*nv + 2.d0)*((f_star**2)*q3 - &
             (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
             (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0))

        df1depsilonv=(3.d0*E*f_star*q1*q2*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
         3.d0)))/(2.d0*Y))*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - (f_star**2)*q3 +&
          2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) -&
           1.d0)**(0.5d0))**m*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/(2.d0*(Y**3)*(6.d0*nv - &
           3.d0)*((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - (f_star**2)*q3 + &
           2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - &
           1.d0)**(1.5d0)) - (d_e*((3.d0*E*(f_star**3)*(q1**3)*(q2**3)*sinh((q2*(3.d0*p_star + &
           (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**3)/(2.d0*(Y**3)*(6.d0*nv - &
           3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
            2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + &
            1.d0)**2) + (3.d0*E*f_star*q1*q2*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
             3.d0)))/(2.d0*Y))*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2)/((Y**5)*(6.d0*nv - &
             3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
             2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) &
             + 1.d0)**2) + (3.d0*E*(f_star**2)*(q1**2)*(q2**3)*cosh((q2*(3.d0*p_star + &
             (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*sinh((q2*(3.d0*p_star + &
             (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)))/(2.d0*(Y**3)*(6.d0*nv - &
             3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
             2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
             3.d0)))/(2.d0*Y)) + 1.d0))))/(2.d0*dt*epsilon_0*(- (sigma_star - &
             (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 - (sigma_star - &
             (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
             (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) - &
             ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
             3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv +&
              2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
              3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0)) + (3.d0*E*f_star*m*q1*q2*sinh((q2*(3.d0*p_star + &
              (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
              2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
              (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**(m - 1.d0)*(sigma_star -&
               (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/(2.d0*(Y**3)*(6.d0*nv - 3.d0)*((f_star**2)*q3 - &
               (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
                + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0))

          df2de =(d_epsilonv*((6.d0*E*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2)))/((Y**4)*(2.d0*nv +&
           2.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
           2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + &
           1.d0)) + (6.d0*E*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**3)/((Y**6)*(2.d0*nv + &
           2.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
           2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + &
           1.d0)**2) + (3.d0*E*(f_star**2)*(q1**2)*(q2**2)*(sinh((q2*(3.d0*p_star + &
           (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**2)*(sigma_star - & ! blan
           (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/((Y**4)*(2.d0*nv + 2.d0)*((f_star**2)*q3 - &
           (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
           + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)**2)))/(2.d0*dt*epsilon_0*(- &
           (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 - (sigma_star - &
           (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
           (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) -&
            ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
            3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
            2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
            3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0)) - (9.d0*E*f_star*q1*q2*sinh((q2*(3.d0*p_star + &
            (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
            2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
            (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**m*(sigma_star - &
            (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/(2.d0*(Y**3)*(2.d0*nv + 2.d0)*(- (f_star**2)*q3 + &
            (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star &
            + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(1.5d0)) - &
            (9.d0*E*f_star*m*q1*q2*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
            3.d0)))/(2.d0*Y))*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
            (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
             3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**(m - 1.d0)*(sigma_star - &
             (3.d0*E*d_e)/(2.d0*nv + 2.d0)))/(2.d0*(Y**3)*(2.d0*nv + 2.d0)*((f_star**2)*q3 - &
             (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
             2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0))


          df2depsilonv =(- (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/((Y**4)*((f_star**2)*q3 -&
           (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
           (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)) - &
           ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
           3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
           2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
           (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0)))**(0.5d0)/(dt*epsilon_0) - &
           (d_epsilonv*((3.d0*E*(f_star**3)*(q1**3)*(q2**3)*sinh((q2*(3.d0*p_star + &
            (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**3)/(2.d0*(Y**3)*(6.d0*nv - &
            3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
            2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
            3.d0)))/(2.d0*Y)) + 1.d0)**2) + (3.d0*E*f_star*q1*q2*sinh((q2*(3.d0*p_star + &
            (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*(sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
            2.d0))**2)/((Y**5)*(6.d0*nv - 3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
            2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - &
            3.d0)))/(2.d0*Y)) + 1.d0)**2) + (3.d0*E*(f_star**2)*(q1**2)*(q2**3)*cosh((q2*(3.d0*p_star + &
            (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*sinh((q2*(3.d0*p_star +&
             (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)))/(2.d0*(Y**3)*(6.d0*nv - &
             3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) -&
              2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + &
              1.d0))))/(2.d0*dt*epsilon_0*(- (sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
              2.d0))**2/((Y**4)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv +&
               2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv -&
                3.d0)))/(2.d0*Y)) + 1.d0)) - ((f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + &
                (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**2)/(2.d0*(Y**2)*((f_star**2)*q3 - &
                (sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - &
                2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) &
                 + 1.d0)))**(0.5d0)) - (9.d0*E*f_star*q1*(q2**2)*cosh((q2*(3.d0*p_star + &
                 (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
                 2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
              (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**m)/(4.d0*(Y**2)*(6.d0*nv &
              - 3.d0)*((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - (f_star**2)*q3 + &
              2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - &
              1.d0)**(0.5d0)) + (9.d0*E*(f_star**2)*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star +&
               (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**2*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv +&
                2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
                (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**m)/(4.d0*(Y**2)*(6.d0*nv &
                - 3.d0)*((sigma_star - (3.d0*E*d_e)/(2.d0*nv + 2.d0))**2/(Y**2) - (f_star**2)*q3 + &
                2.d0*f_star*q1*cosh((q2*(3.d0*p_star + (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - &
                1.d0)**(1.5d0)) + (9.d0*E*(f_star**2)*m*(q1**2)*(q2**2)*sinh((q2*(3.d0*p_star + &
                (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y))**2*(((sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
                2.d0))**2/(Y**2) - (f_star**2)*q3 + 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
                (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) - 1.d0)**(0.5d0))**(m - &
                1.d0))/(4.d0*(Y**2)*(6.d0*nv - 3.d0)*((f_star**2)*q3 - (sigma_star - (3.d0*E*d_e)/(2.d0*nv + &
                2.d0))**2/(Y**2) - 2.d0*f_star*q1*cosh((q2*(3.d0*p_star + &
                (3.d0*E*d_epsilonv)/(6.d0*nv - 3.d0)))/(2.d0*Y)) + 1.d0))


          deriva(1,1)=df1de
          deriva(1,2)=df1depsilonv
          deriva(2,1)=df2de
          deriva(2,2)=df2depsilonv

          resid(1,1)=-f1
          resid(2,1)=-f2

         call invert_small(deriva,deriva_i,deriva_det)
         incre=matmul(deriva_i,resid)
         target=target+incre
         ncount=ncount+1

        enddo

        d_e=target(1,1)
        d_epsilonv=target(2,1)
!================update state variables===========================
       stress1_m=S_star-((d_e*E*3.d0)/(2.d0*(1.d0+nv)*sigma_estar))*S_star+(p_star- &
       (E*d_epsilonv)/(3.d0*(1-2.d0*nv)))*I

        stress1(1)=stress1_m(1,1)
        stress1(2)=stress1_m(2,2)
        stress1(3)=stress1_m(3,3)
        stress1(4)=stress1_m(1,2)
        stress1(5)=stress1_m(1,3)
        stress1(6)=stress1_m(2,3)

        dematrix =(((sigma_star - (3*E*d_e)/(2*nv + 2))**2/(Y**2*(- f_star**2*q3 + &
        (sigma_star - (3*E*d_e)/(2*nv + 2))**2/Y**2 + 2*f_star*q1*cosh((q2*(3*p_star +&
         (3*E*d_epsilonv)/(6*nv - 3)))/(2*Y)) - 1)**(1/2)) + (f_star*q1*q2*sinh((q2*(3*p_star &
         + (3*E*d_epsilonv)/(6*nv - 3)))/(2*Y))*(p_star + (E*d_epsilonv)/(6*nv - 3)))/(2*Y*(-&
          f_star**2*q3 + (sigma_star - (3*E*d_e)/(2*nv + 2))**2/Y**2 + 2*f_star*q1*cosh((q2*(3*p_star &
          + (3*E*d_epsilonv)/(6*nv - 3)))/(2*Y)) - 1)**(1/2)))*(((sigma_star - &
          (3*E*d_e)/(2*nv + 2))**2/Y**2 - f_star**2*q3 + 2*f_star*q1*cosh((q2*(3*p_star + &
          (3*E*d_epsilonv)/(6*nv - 3)))/(2*Y)) - 1)**(1/2))**m)/(- (sigma_star -&
           (3*E*d_e)/(2*nv + 2))**2/(Y**4*(f_star**2*q3 - (sigma_star - (3*E*d_e)/(2*nv + &
           2))**2/Y**2 - 2*f_star*q1*cosh((q2*(3*p_star + (3*E*d_epsilonv)/(6*nv - 3)))/(2*Y)) +&
            1)) - (f_star**2*q1**2*q2**2*sinh((q2*(3*p_star + (3*E*d_epsilonv)/(6*nv - &
            3)))/(2*Y))**2)/(2*Y**2*(f_star**2*q3 - (sigma_star - (3*E*d_e)/(2*nv + 2))**2/Y**2 -&
             2*f_star*q1*cosh((q2*(3*p_star + (3*E*d_epsilonv)/(6*nv - 3)))/(2*Y)) + 1)))**(1/2)



        dematrix=dematrix*dt*epsilon_0/(1.d0-vf)

        vf=1.d0+(vf-1.d0)*(exp(-d_epsilonv))+((f_n*dematrix)/(s_n*sqrt(2.d0*pi)))* &
          exp(-0.5d0*(((ematrix-epsilon_n)/s_n)**2))

       ematrix=ematrix+dematrix



        updated_state_variables(1:6)=stress1(1:6)
        updated_state_variables(7)=ematrix
        updated_state_variables(8)=vf

        end if

    write(IOW,*) ncount

    end subroutine gurson




     subroutine calculate_stress_D(strain,element_properties,stress,D)
        use Types
        use ParamIO
        implicit none

        real (prec), intent(in)    :: strain(:)
        real (prec), intent(in)    :: element_properties(:)
        real (prec), intent(out)   :: stress(:)
        real (prec), intent(out)   :: D(:,:)
!==============local variables================
        integer      :: i,j
        real (prec)  :: sigma_0,e_0,n,k
        real (prec)  :: ekk,ee,sigma_e,E_t,E_s,strainmag
        real (prec)  :: estrain(6)
        real (prec)  :: temp1(6,6),temp2(6,6),e_dyadic_e(6,6)

       !read material parameters
        sigma_0=element_properties(1)
        e_0=element_properties(2)
        n=element_properties(3)
        k=element_properties(4)
 !=================calculate stress===============
   strainmag = dot_product(strain,strain)
 if(strainmag==0.d0) then

        stress=0.d0
 !       write(IOW,*) stress
 else
        !calculate deviatoric strain
        estrain=strain
        ekk=strain(1)+strain(2)+strain(3)
        estrain(1)=estrain(1)-ekk/3.d0
        estrain(2)=estrain(2)-ekk/3.d0
        estrain(3)=estrain(3)-ekk/3.d0
        estrain(4)=estrain(4)/2.d0
        estrain(5)=estrain(5)/2.d0
        estrain(6)=estrain(6)/2.d0
        !calculate von-misers effective strain
        ee=0.d0
        ee=estrain(1)*estrain(1)+estrain(2)*estrain(2)+estrain(3)*estrain(3) &
        +2.d0*(estrain(4)*estrain(4)+estrain(5)*estrain(5)+estrain(6)*estrain(6))
        ee=dsqrt(ee*2.d0/3.d0)
        !calculate σe,Es =σe /ee ,  Et = dσe /dee
        if(ee<e_0) then
        sigma_e=sigma_0*(dsqrt((1.d0+n*n)/((n-1.d0)*(n-1.d0))-(n/(n-1.d0)-ee/e_0)*(n/(n-1.d0)-ee/e_0))-1.d0/(n-1.d0))
        else
        !sigma_e=sigma_0*((ee/e_0)**(1/n))
        sigma_e=sigma_0*exp(log(ee/e_0)/n)
        endif
     !calculate σij
     stress=0.d0
     stress(1:6)=2.d0*sigma_e*estrain(1:6)/(3.d0*ee)
     stress(1:3)=stress(1:3)+k*ekk
 endif
 !===============calculate tangent stiffness matrix [D]====================

 temp1=0.d0
 temp1(1,1)=2
 temp1(2,2)=2
 temp1(3,3)=2
 temp1(4,4)=1
 temp1(5,5)=1
 temp1(6,6)=1
 temp2=0.d0
 temp2(1:3,1:3)=1
 !write(IOW,*) temp1, temp2

 if(strain(1)==0.d0 .and.strain(2)==0.d0 .and.strain(3)==0.d0 .and. &
   strain(4)==0.d0 .and.strain(5)==0.d0 .and.strain(6)==0.d0 ) then

   ee=0.d0
   E_s=(sigma_0/e_0)*(n/(n-1.d0)-ee/e_0)/sqrt((1.d0+n*n)/((n-1.d0)*(n-1.d0))-(n/(n-1.d0)-ee/e_0)*(n/(n-1.d0)-ee/e_0))
   D(1:6,1:6)=(E_s/3.d0)*temp1(1:6,1:6)+(k-2.d0*E_s/9.d0)*temp2(1:6,1:6)

 else
     E_s=sigma_e/ee
     if(ee<e_0) then
     E_t=(sigma_0/e_0)*(n/(n-1.d0)-ee/e_0)/sqrt((1.d0+n*n)/((n-1.d0)*(n-1.d0))-(n/(n-1.d0)-ee/e_0)*(n/(n-1.d0)-ee/e_0))
     else
     !E_t=(sigma_0/e_0)*((ee/e_0)**(1/n-1))/n
     E_t=(sigma_0/e_0)*exp(log(ee/e_0)*(1/n-1.d0))/n
     endif
     e_dyadic_e = spread(estrain,dim=2,ncopies=6)*spread(estrain,dim=1,ncopies=6)
     D(1:6,1:6)=(4.d0*(E_t-E_s)/(9.d0*ee*ee))*e_dyadic_e(1:6,1:6)+(E_s/3.d0)*temp1(1:6,1:6)+(k-2.d0*E_s/9.d0)*temp2(1:6,1:6)

 endif

 !  write(IOW,*) D

    end subroutine calculate_stress_D





    function cross_product(a,b)           ! Compute cross product of two 3 dimensional vectors

        use Types
        implicit none

        real (prec), intent(in)  :: a(3)
        real (prec), intent(in)  :: b(3)

        real (prec) :: cross_product(3)

        cross_product(1) = a(2)*b(3)-a(3)*b(2)
        cross_product(2) = b(1)*a(3)-a(1)*b(3)
        cross_product(3) = a(1)*b(2)-a(2)*b(1)

    end function cross_product

    subroutine invert_small(A,A_inverse,determinant)  ! Compute the inverse and determinant of a 3x3 or 2x2 matrix
        use Types
        use ParamIO
        implicit none

        real (prec), intent(in)    :: A(:,:)
        real (prec), intent(out)   :: determinant
        real (prec), intent(out)   :: A_inverse(:,:)

        real (prec) :: cofactor(3,3)
!
!       A (3x3) or (2x2) input matrix
!       A_inverse (3x3) or (2x2) inverse
!       determinant - determinant of matrix
!
        if (size(A)==4) then
            determinant = A(1,1)*A(2,2)-A(2,1)*A(1,2)
            A_inverse(1,1) = A(2,2)
            A_inverse(2,2) = A(1,1)
            A_inverse(1,2) = -A(1,2)
            A_inverse(2,1) = -A(2,1)
            IF (determinant==0.d0) THEN
                write(IOW,*) ' Error in element utility invert'
                write(IOW,*) ' A 2x2 matrix has a zero determinant'
                stop
            endif
            A_inverse = A_inverse/determinant
        else if (size(A)==9) then

            determinant =   A(1,1)*A(2,2)*A(3,3)  &
                - A(1,1)*A(2,3)*A(3,2)  &
                - A(1,2)*A(2,1)*A(3,3)  &
                + A(1,2)*A(2,3)*A(3,1)  &
                + A(1,3)*A(2,1)*A(3,2)  &
                - A(1,3)*A(2,2)*A(3,1)

            IF (determinant==0.d0) THEN
                write(IOW,*) ' Error in element utility invert'
                write(IOW,*) ' A 3x3 matrix has a zero determinant'
                stop
            endif
            COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
            COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
            COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
            COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
            COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
            COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
            COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
            COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

            A_inverse = transpose(COFACTOR) / determinant
        endif

    end subroutine invert_small

    subroutine eigenvecs33(A,eigenvalues,eigenvectors)     ! Compute eigenvals and eigenvectors of symmetric 3x3 matrix

        use Types
        implicit none

        real (prec), intent(in)  ::  A(3,3)                ! Input matrix
        real (prec), intent(out) :: eigenvalues(3)         ! Eigenvalues
        real (prec), intent(out) :: eigenvectors(3,3)      ! ith eigenvector is eigenvectors(1:3,i)

        real (prec) :: B(3,3)
        real (prec) :: C(3,3)
        real (prec) :: D(3,3)
        real (prec) :: p1,p2
        real (prec) :: p,q,r
        real (prec) :: phi
        real (prec) :: evnorm
        real (prec) :: tol

        integer :: i

        eigenvectors = 0.d0
        p1 = A(1,2)*A(1,2) + A(1,3)*A(1,3) + A(2,3)*A(2,3)
        tol = 1.d-08*(A(1,1)*A(1,1) + A(2,2)*A(2,2) + A(3,3)*A(3,3) + p1)
        if (p1 == 0.d0) then
           eigenvalues(1) = A(1,1)
           eigenvalues(2) = A(2,2)
           eigenvalues(3) = A(3,3)

           eigenvectors(1,1) = 1.d0
           eigenvectors(2,2) = 1.d0
           eigenvectors(3,3) = 1.d0
        else
          q = (A(1,1)+A(2,2)+A(3,3))/3.d0
          p2 = (A(1,1) - q)*(A(1,1)-q) + (A(2,2) - q)*(A(2,2) - q) + (A(3,3) - q)*(A(3,3) - q) + 2.d0 * p1
          p = dsqrt(p2 / 6.d0)
          B = (1.d0 / p) * (A - q * eye3_d)      ! eye3_d is the double precision 3x3 identity matrix
          r =   0.5d0*( B(1,1)*B(2,2)*B(3,3)  &
                      - B(1,1)*B(2,3)*B(3,2)  &
                      - B(1,2)*B(2,1)*B(3,3)  &
                      + B(1,2)*B(2,3)*B(3,1)  &
                      + B(1,3)*B(2,1)*B(3,2)  &
                      - B(1,3)*B(2,2)*B(3,1) )

           if (r < -1.d0) then
              phi = PI_D / 3.d0
           else if (r > 1.d0) then
              phi = 0.d0
           else
              phi = dacos(r) / 3.d0
           end if

           ! the eigenvalues satisfy eig3 <= eig2 <= eig1
           eigenvalues(1) = q + 2.d0 * p * dcos(phi)
           eigenvalues(3) = q + 2.d0 * p * dcos(phi + (2.d0*PI_D/3.d0))
           eigenvalues(2) = 3.d0 * q - eigenvalues(1) - eigenvalues(3)

           do i = 1,3
              B = A - eigenvalues(i)*eye3_D
              C = A - eigenvalues(mod(i,3)+1)*eye3_D
              D = matmul(B,C)
              eigenvectors(1:3,mod(i+1,3)+1) = matmul(D,(/1.d0,1.d0,1.d0/))
           end do
           do i = 1,3
              evnorm = dsqrt(dot_product(eigenvectors(1:3,i),eigenvectors(1:3,i)))
              if (evnorm>tol) then
                 eigenvectors(1:3,i) = eigenvectors(1:3,i)/evnorm
              else
                 eigenvectors(1:3,i) = cross_product(eigenvectors(1:3,mod(i,3)+1),eigenvectors(1:3,mod(i+1,3)+1))
                 evnorm = dsqrt(dot_product(eigenvectors(1:3,i),eigenvectors(1:3,i)))
                 eigenvectors(1:3,i) = eigenvectors(1:3,i)/evnorm
              endif
            end do
         end if
    end subroutine eigenvecs33

    subroutine calculate_principalvals(symvec,principalstress)     ! Compute principal stresses or strains from a stress or strain vector
        use Types
        use ParamIO
        implicit none

        real (prec), intent(in)  ::  symvec(6)   ! symvec contains (A(1,1),A(2,2),A(3,3),A(1,2),A(1,3),A(2,3))
        real (prec), intent(out) :: principalstress(3)
       ! real (prec) ::  principalvals33(3)

        real (prec) :: B(6)
        real (prec) :: p1,p2
        real (prec) :: p,q,r
        real (prec) :: phi

        principalstress = 0.d0
        p1 = symvec(4)*symvec(4)+symvec(5)*symvec(5)+symvec(6)*symvec(6)
        if (p1 == 0.d0) then    ! Vector is diagonal - sort the eigenvalues
           principalstress(1) = maxval(symvec(1:3))
           principalstress(3) = minval(symvec(1:3))
           principalstress(2) = sum(symvec(1:3)) - principalstress(1) - principalstress(2)
        else
          q = sum(symvec(1:3))/3.d0
          p2 = (symvec(1) - q)*(symvec(1)-q) + (symvec(2) - q)*(symvec(2) - q) + (symvec(3) - q)*(symvec(3) - q) + 2.d0 * p1
          p = dsqrt(p2 / 6.d0)
          B = symvec/p
          B(1:3) = B(1:3) - q/p
          r = 0.5d0*(  B(1)*B(2)*B(3)  &
               - B(1)*B(6)*B(6)  &
               - B(4)*B(4)*B(3)  &
               + B(4)*B(6)*B(5)  &
               + B(5)*B(4)*B(6)  &
               - B(5)*B(2)*B(5) )

           if (r < -1.d0) then
              phi = PI_D / 3.d0
           else if (r > 1.d0) then
              phi = 0.d0
           else
              phi = dacos(r) / 3.d0
           end if

           ! Principal values ordered from largest to smallest.
           principalstress(1) = q + 2.d0 * p * cos(phi)
           principalstress(3) = q + 2.d0 * p * cos(phi + (2.d0*PI_D/3.d0))
           principalstress(2) = 3.d0 * q - principalstress(1) - principalstress(3)
        endif

     end subroutine calculate_principalvals

     function rotatesymvec(symvec,R)     ! Apply a rigid rotation to a symmetric 3x3 stress or strain vector

         use Types
         implicit none

         real (prec), intent (in) :: symvec(6)       ! Stress or strain vector stored as [s11,s22,s33,s12,s13,s23]
         real (prec), intent (in) :: R(3,3)          ! Rotation tensor stored as 3x3 matrix

         real (prec) :: rotatesymvec(6)              ! R*symvec*R^T, stored as [s11,s22,s33,s12,s13,s23]

         real (prec) :: temp1(3,3)
         real (prec) :: temp2(3,3)

         temp1(1,1) = symvec(1)
         temp1(2,2) = symvec(2)
         temp1(3,3) = symvec(3)
         temp1(1,2) = symvec(4)
         temp1(1,3) = symvec(5)
         temp1(2,3) = symvec(6)
         temp1(2,1) = temp1(1,2)
         temp1(3,1) = temp1(1,3)
         temp1(3,2) = temp1(2,3)

         temp2 = matmul(R,matmul(temp1,transpose(R)))

         rotatesymvec(1) = temp2(1,1)
         rotatesymvec(2) = temp2(2,2)
         rotatesymvec(3) = temp2(3,3)
         rotatesymvec(4) = temp2(1,2)
         rotatesymvec(5) = temp2(1,3)
         rotatesymvec(6) = temp2(2,3)

     end function rotatesymvec

     function sqrtM33(A)                                     ! Square root of a symmetric 3x3 matrix

        use Types
        implicit none

        real (prec), intent(in) :: A(3,3)
        real (prec) :: sqrtM33(3,3)

        real (prec) :: D(3,3)
        real (prec) :: V(3,3)
        real (prec) :: eig(3)

        call eigenvecs33(A,eig,V)

        D = 0.d0
        D(1,1) = dsqrt(eig(1))
        D(2,2) = dsqrt(eig(2))
        D(3,3) = dsqrt(eig(3))

        sqrtM33 = matmul(V,matmul(D,transpose(V)))

     end function sqrtM33

    subroutine polardecomp(A,V,U,R)                ! Compute left and right polar decompositions of a 3x3 matrix A

       use Types
       implicit none

       real (prec), intent (in)   :: A(3,3)
       real (prec), intent (out)  :: V(3,3)
       real (prec), intent (out)  :: U(3,3)
       real (prec), intent (out)  :: R(3,3)

       real (prec) :: det

       !  Decompose A into A=RU=VR  where U,V are symmetric and R is orthogonal

       R = matmul(A,transpose(A))    ! R is just temporary variable here
       V = sqrtM33(R)                ! V= sqrt(A*A^T)
       call invert_small(V,U,det)         ! U is just temporary variable here
       R = matmul(U,A)               ! R = V^-1*A
       U = matmul(transpose(R),A)    ! U = R^T*A

    end subroutine polardecomp


    subroutine initialize_integration_points(n_points, n_nodes, xi, w)
        use Types
        use ParamIO
        implicit none

        integer, intent( in )       :: n_nodes             ! No. nodes on element
        integer, intent( in )       :: n_points            ! No. integration points
        real( prec ), intent( out ) :: xi(:,:)             ! xi(i,a) ith coord of ath integration point
        real( prec ), intent( out ) :: w(:)                ! Integration weight

        real (prec) :: x1D(4),w1D(4), cn, w1, w2, w11, w12, w22
        integer :: i,j,k,n
  
        if (size(xi,1)==6) then                        ! Integration points for 1d elements
  
            select case ( n_points )
                case (2)
                    xi(1,1) = .5773502691896257D+00
                    xi(2,1) = -.5773502691896257D+00
                    w(1) = .1000000000000000D+01
                    w(2) = .1000000000000000D+01
                    return
                case (3)
                    xi(3,1) = 0.7745966692414834D+00
                    xi(2,1) = .0000000000000000D+00
                    xi(3,1) = -.7745966692414834D+00
                    w(1) = .5555555555555556D+00
                    w(2) = .8888888888888888D+00
                    w(3) = .5555555555555556D+00
                    return
                case (4)
                    xi(1,1) = .8611363115940526D+00
                    xi(2,1) = .3399810435848563D+00
                    xi(3,1) = -.3399810435848563D+00
                    xi(4,1) = -.8611363115940526D+00
                    w(1) = .3478548451374538D+00
                    w(2) = .6521451548625461D+00
                    w(3) = .6521451548625461D+00
                    w(4) = .3478548451374538D+00
                    return
                case (5)
                    xi(1,1) = .9061798459386640D+00
                    xi(2,1) = .5384693101056831D+00
                    xi(3,1) = .0000000000000000D+00
                    xi(4,1) = -.5384693101056831D+00
                    xi(5,1) = -.9061798459386640D+00
                    w(1) = .2369268850561891D+00
                    w(2) = .4786286704993665D+00
                    w(3) = .5688888888888889D+00
                    w(4) = .4786286704993665D+00
                    w(5) = .2369268850561891D+00
                    return
                case (6)
                    xi(1,1) = .9324695142031521D+00
                    xi(2,1) = .6612093864662645D+00
                    xi(3,1) = .2386191860831969D+00
                    xi(4,1) = -.2386191860831969D+00
                    xi(5,1) = -.6612093864662645D+00
                    xi(6,1) = -.9324695142031521D+00
                    w(1) = .1713244923791703D+00
                    w(2) = .3607615730481386D+00
                    w(3) = .4679139345726910D+00
                    w(4) = .4679139345726910D+00
                    w(5) = .3607615730481386D+00
                    w(6) = .1713244923791703D+00
                    return
                case DEFAULT
                    write(IOW,*) ' Error in subroutine initialize_integration_points'
                    write(IOW,*) ' Invalide number of integration points for a 1D integration scheme'
                    write(IOW,*) ' n_points must be between 1 and 6'
                    stop
            end select
        else if (size(xi,1)==2) then                   ! Integration points for 2d elements
            if ( n_nodes>9) then
                write (IOW, 99001) n_nodes
                stop
            end if
            if ( n_points>9 ) then
                write (IOW, 99002) n_points
                stop
            end if
            if ( n_points==1 ) then
                if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
                    xi(1, 1) = 0.D0
                    xi(2, 1) = 0.D0
                    w(1) = 4.D0
                else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
                    xi(1, 1) = 1.D0/3.D0
                    xi(2, 1) = 1.D0/3.D0
                    w(1) = 1.D0/2.D0
                end if
            else if ( n_points==3 ) then
                xi(1, 1) = 0.5D0
                xi(2, 1) = 0.5D0
                w(1) = 1.D0/6.D0
                xi(1, 2) = 0.D0
                xi(2, 2) = 0.5D0
                w(2) = w(1)
                xi(1, 3) = 0.5D0
                xi(2, 3) = 0.D0
                w(3) = w(1)
            else if ( n_points==4 ) then
                if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
                    !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
                    !     43
                    !     12
                    cn = 0.5773502691896260D0
                    xi(1, 1) = -cn
                    xi(1, 2) = cn
                    xi(1, 3) = cn
                    xi(1, 4) = -cn
                    xi(2, 1) = -cn
                    xi(2, 2) = -cn
                    xi(2, 3) = cn
                    xi(2, 4) = cn
                    w(1) = 1.D0
                    w(2) = 1.D0
                    w(3) = 1.D0
                    w(4) = 1.D0
                else if ( n_nodes==3 .or. n_nodes==6 ) then
                    !     xi integration points for triangle
                    xi(1, 1) = 1.D0/3.D0
                    xi(2, 1) = xi(1, 1)
                    w(1) = -27.D0/96.D0
                    xi(1, 2) = 0.6D0
                    xi(2, 2) = 0.2D0
                    w(2) = 25.D0/96.D0
                    xi(1, 3) = 0.2D0
                    xi(2, 3) = 0.6D0
                    w(3) = w(2)
                    xi(1, 4) = 0.2D0
                    xi(2, 4) = 0.2D0
                    w(4) = w(2)
                end if

            else if ( n_points==7 ) then
                ! Quintic integration for triangle
                xi(1,1) = 1.d0/3.d0
                xi(2,1) = xi(1,1)
                w(1) = 0.1125d0
                xi(1,2) = 0.0597158717d0
                xi(2,2) = 0.4701420641d0
                w(2) = 0.0661970763d0
                xi(1,3) = xi(2,2)
                xi(2,3) = xi(1,2)
                w(3) = w(2)
                xi(1,4) = xi(2,2)
                xi(2,4) = xi(2,2)
                w(4) = w(2)
                xi(1,5) = 0.7974269853d0
                xi(2,5) = 0.1012865073d0
                w(5) = 0.0629695902d0
                xi(1,6) = xi(2,5)
                xi(2,6) = xi(1,5)
                w(6) = w(5)
                xi(1,7) = xi(2,5)
                xi(2,7) = xi(2,5)
                w(7) = w(5)
            else if ( n_points==9 ) then
                !     3X3 GAUSS INTEGRATION POINTS IN PSI-ETA COORDINATES
                !     789
                !     456
                !     123
                cn = 0.7745966692414830D0
                xi(1, 1) = -cn
                xi(1, 2) = 0.D0
                xi(1, 3) = cn
                xi(1, 4) = -cn
                xi(1, 5) = 0.D0
                xi(1, 6) = cn
                xi(1, 7) = -cn
                xi(1, 8) = 0.D0
                xi(1, 9) = cn
                xi(2, 1) = -cn
                xi(2, 2) = -cn
                xi(2, 3) = -cn
                xi(2, 4) = 0.D0
                xi(2, 5) = 0.D0
                xi(2, 6) = 0.D0
                xi(2, 7) = cn
                xi(2, 8) = cn
                xi(2, 9) = cn
                w1 = 0.5555555555555560D0
                w2 = 0.8888888888888890D0
                w11 = w1*w1
                w12 = w1*w2
                w22 = w2*w2
                w(1) = w11
                w(2) = w12
                w(3) = w11
                w(4) = w12
                w(5) = w22
                w(6) = w12
                w(7) = w11
                w(8) = w12
                w(9) = w11
            end if

99001       format ( // ' *** ERROR ***'/  &
                '  Number of nodes on 2d element is greater ',  &
                'than 9 '/'  No. nodes = ', I5)
99002       format ( // ' *** ERROR ***'/,  &
                '  Number of int. pts on element is greater ',  &
                'than 9 '/'  No. int pts. = ', I5)
  
   
        else if (size(xi,1)==3) then                  ! Integration points for 3d elements
            if (n_nodes  == 4.or.n_nodes ==10 ) then
                if (n_points == 1) then
                    xi(1,1) = 0.25
                    xi(2,1) = 0.25
                    xi(3,1) = 0.25
                    w(1) = 1.D0/6.D0
                else if (n_points == 4) then
                    xi(1,1) = 0.58541020
                    xi(2,1) = 0.13819660
                    xi(3,1) = xi(2,1)
                    xi(1,2) = xi(2,1)
                    xi(2,2) = xi(1,1)
                    xi(3,2) = xi(2,1)
                    xi(1,3) = xi(2,1)
                    xi(2,3) = xi(2,1)
                    xi(3,3) = xi(1,1)
                    xi(1,4) = xi(2,1)
                    xi(2,4) = xi(2,1)
                    xi(3,4) = xi(2,1)
                    w(1:4) = 1.D0/24.D0
                else if (n_points == 5) then
                    xi(1,1) = 0.25d0
                    xi(2,1) = 0.25d0
                    xi(3,1) = 0.25d0
                    xi(1,2) = 0.5d0
                    xi(2,2) = 1.d0/6.d0
                    xi(3,2) = 1.d0/6.d0
                    xi(1,3) = 1.d0/6.d0
                    xi(2,3) = 0.5d0
                    xi(3,3) = 1.d0/6.d0
                    xi(1,4) = 1.d0/6.d0
                    xi(2,4) = 1.d0/6.d0
                    xi(3,4) = 0.5d0
                    xi(1,5) = 1.d0/6.d0
                    xi(2,5) = 1.d0/6.d0
                    xi(3,5) = 1.d0/6.d0
                    w(1) = -4.d0/30.d0
                    w(2:5) = 3.d0/40.d0
                else
                    write(IOW,*) ' Incorrect number of integration points for tetrahedral element '
                    write(IOW, *) ' called with ',n_points
                    stop
                endif
            else if ( n_nodes == 8 .or. n_nodes == 20 ) then
                if (n_points == 1) then
                    xi(1,1) = 0.D0
                    xi(2,1) = 0.D0
                    xi(3,1) = 0.D0
                    w(1) = 8.D0
                else if (n_points == 8) then
                    x1D(1) = -0.5773502692
                    x1D(2) =  0.5773502692
                    do k = 1,2
                        do j = 1,2
                            do i = 1,2
                                n = 4*(k-1) + 2*(j-1) + i
                                xi(1,n) = x1D(i)
                                xi(2,n) = x1D(j)
                                xi(3,n) = x1D(k)
                            end do
                        end do
                    end do
                    w(1:8) = 1.D0
                else if (n_points == 27) then
                    x1D(1) = -0.7745966692
                    x1D(2) = 0.
                    x1D(3) = 0.7745966692
                    w1D(1) = 0.5555555555D0
                    w1D(2) = 0.888888888D0
                    w1D(3) = 0.55555555555D0
                    do k = 1,3
                        do j = 1,3
                            do i = 1,3
                                n = 9*(k-1) + 3*(j-1) + i
                                xi(1,n) = x1D(i)
                                xi(2,n) = x1D(j)
                                xi(3,n) = x1D(k)
                                w(n) = w1D(i)*w1D(j)*w1D(k)
                            end do
                        end do
                    end do
                else if (n_points == 64) then
                    x1D(1) = .8611363115940526D+00
                    x1D(2) = .3399810435848563D+00
                    x1D(3) = -.3399810435848563D+00
                    x1D(4) = -.8611363115940526D+00
                    w1D(1) = .3478548451374538D+00
                    w1D(2) = .6521451548625461D+00
                    w1D(3) = .6521451548625461D+00
                    w1D(4) = .3478548451374538D+00
                    do k = 1,4
                        do j = 1,4
                            do i = 1,4
                                n = 16*(k-1) + 4*(j-1) + i
                                xi(1,n) = x1D(i)
                                xi(2,n) = x1D(j)
                                xi(3,n) = x1D(k)
                                w(n) = w1D(i)*w1D(j)*w1D(k)
                            end do
                        end do
                    end do
                endif
            endif
        endif
    end subroutine initialize_integration_points
  


    subroutine calculate_shapefunctions(xi,n_nodes,f,df)
        use Types
        implicit none
    
        integer, intent(in)   :: n_nodes
        real (prec), intent(in)  :: xi(:)         ! Integration point coordinates
        real (prec), intent(out) :: f(:)          ! Shape functions
        real (prec), intent(out) :: df(:,:)       ! Shape function derivatives
    
        real (prec) :: z, g1, g2, g3, h1, h2, h3, dg1, dg2, dh1, dh2, dzdp,dzdq, xi4
    
        if (size(df)==3) then    ! 1D shape functions
    
            if (n_nodes==2) then
                f(1) = 0.5d0*(1.d0-xi(1))
                f(2) = 0.5d0*(1.d0+xi(1))
                df(1,1) = -0.5d0
                df(2,1) =  0.5d0
            else if (n_nodes==3) then
                f(1) = -0.5*xi(1)*(1.-xi(1))
                f(2) =  0.5*xi(1)*(1.+xi(1))
                f(3) = (1.-xi(1))*(1.+xi(1))
                df(1,1) = -0.5+xi(1)
                df(2,1) =  0.5+xi(1)
                df(3,1) = -2.d0*xi(1)
            endif
        else if (size(df)==18) then  !2D shape functions
            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2
       
            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
            end if
        else if (size(df)==60) then  !3D shape functions
            if (n_nodes == 4) then
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = xi(3)
                f(4) = 1.-xi(1)-xi(2)-xi(3)
                df(1,1) = 1.
                df(2,2) = 1.
                df(3,3) = 1.
                df(4,1) = -1.
                df(4,2) = -1.
                df(4,3) = -1.
            else if (n_nodes == 10) then
                xi4 = 1.-xi(1)-xi(2)-xi(3)
                f(1) = (2.*xi(1)-1.)*xi(1)
                f(2) = (2.*xi(2)-1.)*xi(2)
                f(3) = (2.*xi(3)-1.)*xi(3)
                f(4) = (2.*xi4-1.)*xi4
                f(5) = 4.*xi(1)*xi(2)
                f(6) = 4.*xi(2)*xi(3)
                f(7) = 4.*xi(3)*xi(1)
                f(8) = 4.*xi(1)*xi4
                f(9) = 4.*xi(2)*xi4
                f(10) = 4.*xi(3)*xi4
                df(1,1) = (4.*xi(1)-1.)
                df(2,2) = (4.*xi(2)-1.)
                df(3,3) = (4.*xi(3)-1.)
                df(4,1) = -(4.*xi4-1.)
                df(4,2) = -(4.*xi4-1.)
                df(4,3) = -(4.*xi4-1.)
                df(5,1) = 4.*xi(2)
                df(5,2) = 4.*xi(1)
                df(6,2) = 4.*xi(3)
                df(6,3) = 4.*xi(2)
                df(7,1) = 4.*xi(3)
                df(7,3) = 4.*xi(1)
                df(8,1) = 4.*(xi4-xi(1))
                df(8,2) = -4.*xi(1)
                df(8,3) = -4.*xi(1)
                df(9,1) = -4.*xi(2)
                df(9,2) = 4.*(xi4-xi(2))
                df(9,3) = -4.*xi(2)
                df(10,1) = -4.*xi(3)*xi4
                df(10,2) = -4.*xi(3)
                df(10,3) = 4.*(xi4-xi(3))
            else if (n_nodes == 8) then
                f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.
                f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.
                f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.
                f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.
                f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.
                f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.
                f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.
                f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.
                df(1,1) = -(1.-xi(2))*(1.-xi(3))/8.
                df(1,2) = -(1.-xi(1))*(1.-xi(3))/8.
                df(1,3) = -(1.-xi(1))*(1.-xi(2))/8.
                df(2,1) = (1.-xi(2))*(1.-xi(3))/8.
                df(2,2) = -(1.+xi(1))*(1.-xi(3))/8.
                df(2,3) = -(1.+xi(1))*(1.-xi(2))/8.
                df(3,1) = (1.+xi(2))*(1.-xi(3))/8.
                df(3,2) = (1.+xi(1))*(1.-xi(3))/8.
                df(3,3) = -(1.+xi(1))*(1.+xi(2))/8.
                df(4,1) = -(1.+xi(2))*(1.-xi(3))/8.
                df(4,2) = (1.-xi(1))*(1.-xi(3))/8.
                df(4,3) = -(1.-xi(1))*(1.+xi(2))/8.
                df(5,1) = -(1.-xi(2))*(1.+xi(3))/8.
                df(5,2) = -(1.-xi(1))*(1.+xi(3))/8.
                df(5,3) = (1.-xi(1))*(1.-xi(2))/8.
                df(6,1) = (1.-xi(2))*(1.+xi(3))/8.
                df(6,2) = -(1.+xi(1))*(1.+xi(3))/8.
                df(6,3) = (1.+xi(1))*(1.-xi(2))/8.
                df(7,1) = (1.+xi(2))*(1.+xi(3))/8.
                df(7,2) = (1.+xi(1))*(1.+xi(3))/8.
                df(7,3) = (1.+xi(1))*(1.+xi(2))/8.
                df(8,1) = -(1.+xi(2))*(1.+xi(3))/8.
                df(8,2) = (1.-xi(1))*(1.+xi(3))/8.
                df(8,3) = (1.-xi(1))*(1.+xi(2))/8.
            else if (n_nodes == 20) then
                f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.
                f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.
                f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.
                f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.
                f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.
                f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.
                f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.
                f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.
                f(9)  = (1.-xi(1)**2.)*(1.-xi(2))*(1.-xi(3))/4.
                f(10) = (1.+xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
                f(11) = (1.-xi(1)**2.)*(1.+xi(2))*(1.-xi(3))/4.
                f(12) = (1.-xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
                f(13) = (1.-xi(1)**2.)*(1.-xi(2))*(1.+xi(3))/4.
                f(14) = (1.+xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
                f(15) = (1.-xi(1)**2.)*(1.+xi(2))*(1.+xi(3))/4.
                f(16) = (1.-xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
                f(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
                f(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
                f(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
                f(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
                df(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
                df(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
                df(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
 
                df(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
                df(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
                df(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
 
                df(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
                df(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
                df(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)-(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
 
                df(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
                df(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
                df(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
                df(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
                df(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
                df(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)+(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
                df(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
                df(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
                df(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
                df(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                df(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                df(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                df(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                df(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                df(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                df(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
                df(9,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
                df(9,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
                df(10,1)  = (1.-xi(2)**2.)*(1.-xi(3))/4.
                df(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.
                df(10,3)  = -(1.-xi(2)**2.)*(1.+xi(1))/4.
                df(11,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
                df(11,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
                df(11,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
                df(12,1)  = -(1.-xi(2)**2.)*(1.-xi(3))/4.
                df(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.
                df(12,3)  = -(1.-xi(2)**2.)*(1.-xi(1))/4.
                df(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.
                df(13,2)  = -(1.-xi(1)**2.)*(1.+xi(3))/4.
                df(13,3)  = (1.-xi(1)**2.)*(1.-xi(2))/4.
                df(14,1)  = (1.-xi(2)**2.)*(1.+xi(3))/4.
                df(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.
                df(14,3)  = (1.-xi(2)**2.)*(1.+xi(1))/4.
                df(15,1)  = 2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.
                df(15,2)  = (1.-xi(1)**2.)*(1.+xi(3))/4.
                df(15,3)  = (1.-xi(1)**2.)*(1.+xi(2))/4.
                df(16,1)  = -(1.-xi(2)**2.)*(1.+xi(3))/4.
                df(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.
                df(16,3)  = (1.-xi(2)**2.)*(1.-xi(1))/4.
                df(17,1) = -(1.-xi(2))*(1.-xi(3)**2.)/4.
                df(17,2) = -(1.-xi(1))*(1.-xi(3)**2.)/4.
                df(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.
                df(18,1) = (1.-xi(2))*(1.-xi(3)**2.)/4.
                df(18,2) = -(1.+xi(1))*(1.-xi(3)**2.)/4.
                df(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.
                df(19,1) = (1.+xi(2))*(1.-xi(3)**2.)/4.
                df(19,2) = (1.+xi(1))*(1.-xi(3)**2.)/4.
                df(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.
                df(20,1) = -(1.+xi(2))*(1.-xi(3)**2.)/4.
                df(20,2) = (1.-xi(1))*(1.-xi(3)**2.)/4.
                df(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.
            endif
   
        endif
   
    end subroutine calculate_shapefunctions
   



    subroutine facenodes(ndims,nelnodes,face,list,nfacenodes)
        use Types
        use ParamIO
        implicit none

        integer, intent (in)      :: ndims
        integer, intent (in)      :: nelnodes
        integer, intent (in)      :: face
        integer, intent (out)     :: list(*)
        integer, intent (out)     :: nfacenodes

        integer :: i3(3), i4(4)
        !
        !        Subroutine to return list of nodes on an element face for standard 2D and 3D solid elements
        !
        if (ndims==2) then
            i3(1:3) = [2,3,1]
            i4(1:4) = [2,3,4,1]

            if (nelnodes == 3) then
                nfacenodes = 2
                list(1) = face
                list(2) = i3(face)

            else if (nelnodes == 4) then
                nfacenodes = 2
                list(1) = face
                list(2) = i4(face)
            else if (nelnodes == 6) then
                nfacenodes = 3
                list(1) = face
                list(2) = i3(face)
                list(3) = face+3
            else if (nelnodes == 8) then
                nfacenodes = 3
                list(1) = face
                list(2) = i4(face)
                list(3) = face+4
            else if (nelnodes == 9) then
                nfacenodes = 3
                if (face==1) list(1:3) = (/1,3,2/)
                if (face==2) list(1:3) = (/3,9,6/)
                if (face==3) list(1:3) = (/9,7,8/)
                if (face==4) list(1:3) = (/7,1,4/)
            endif
  
        else if (ndims==3) then

            if (nelnodes == 4) then
                nfacenodes = 3
                if   (face == 1) list(1:3) = [1,2,3]
                if (face == 2) list(1:3) = [1,4,2]
                if (face == 3) list(1:3) = [2,4,3]
                if (face == 4) list(1:3) = [3,4,1]
            else if (nelnodes == 10) then
                nfacenodes = 6
                if   (face == 1) list(1:6) = [1,2,3,5,6,7]
                if (face == 2) list(1:6) = [1,4,2,8,9,5]
                if (face == 3) list(1:6) = [2,4,3,9,10,6]
                if (face == 4) list(1:6) = [3,4,1,10,8,7]
            else if (nelnodes == 8) then
                nfacenodes = 4
                if (face==1) list(1:4) = [1,2,3,4]
                if (face==2) list(1:4) = [5,8,7,6]
                if (face==3) list(1:4) = [1,5,6,2]
                if (face==4) list(1:4) = [2,6,7,3]
                if (face==5) list(1:4) = [3,7,8,4]
                if (face==6) list(1:4) = [4,8,5,1]
            else  if (nelnodes == 20) then
                nfacenodes = 8
                if (face == 1) list(1:8) = [1,2,3,4,9,10,11,12]
                if (face == 2) list(1:8) = [5,8,7,6,16,15,14,13]
                if (face == 3) list(1:8) = [1,5,6,2,17,13,18,9]
                if (face == 4) list(1:8) = [2,6,7,3,18,14,19,10]
                if (face == 5) list(1:8) = [3,7,8,4,19,15,6,11]
                if (face == 6) list(1:8) = [4,8,5,1,20,16,17,12]
            endif
        endif
    end subroutine facenodes



end module
