program fehr

    use variables
    implicit none
    integer :: tj, ts, clock, done
    character(len=1) :: cr = char(13) ! Carriage return character

    ! Print header
    write(6,"(a50,//,23x,a4,23x,//,9x,a31,9x,//,a50,/)") repeat("#",50),&
    "FEHR","Ehrenfest correlation functions",repeat("#",50)

    ! Initialize random number generator
    call system_clock(clock)
    call RLUXGO(4,clock,0,0)

    ! Read input file
    call read_input()

    ! Allocate arrays
    call allocate_arrays()

    ! Set bath frequencies and coupling coefficients
    call bath_properties()

    ! Monte-carlo loop
    write(6, '("Running trajectories:")')
    do tj = 1, ntraj

        ! Sample initial conditions
        call sample_electronic(tj)
        call sample_nuclear(tj)

        ! Calculate time-t operators (t=0)
        call calculate_obs(1)

        ! Molecular dynamics loop
        do ts = 1, tsteps

            ! Make a single step
            if ( intgt == "diagonalise" .or. intgt == "diag" ) then
                call step_diag()
            else if ( intgt == "velocityverlet" .or. intgt == "vv") then
                call step_vverlet()
            else
                write(6,*) "ERROR: integrator type must be either ",&
                           "'diagonalise'/'diag' or 'velocityverlet'/'vv'!"
                stop
            end if

            ! Calculate time-t operators
            call calculate_obs(ts+1)

        end do

        ! Trajectory progress bar
        if ( tj > 1 .and. mod(tj,(ntraj/100)) == 0 ) then
            done = floor(1.d2*tj/dble(ntraj))
            write(6, '(2x,"- Finished: ",i3,"%",1x,52a)', advance="no") done, &
            "[", repeat("#", done/2), repeat(" ", 50-done/2), "]", cr
            flush(6)
        end if

    end do

    ! Average and output observables
    call average_obs()

    ! Deallocate arrays
    call deallocate_arrays()

end program fehr


subroutine read_input()

    use variables
    implicit none
    integer :: iost=1
    character(len=80) :: dum

    open(11, file="input", action="read")
    read(11, *) dum, init_state
    read(11, *) dum, F
    read(11, *) dum, ntraj
    read(11, *) dum, tsteps
    read(11, *) dum, dt
    read(11, *) dum, intgt
    close(11)

    write(6, '("> Input parameters parsed:")')
    write(6, '(2x,a16,1x,i8)') "- Initial state:", init_state
    write(6, '(2x,a12,5x,i8)') "- Bath DoFs:", F
    write(6, '(2x,a15,2x,i8)') "- Trajectories:", ntraj
    write(6, '(2x,a12,5x,i8)') "- Timesteps:", tsteps
    write(6, '(2x,a14,6x,f5.3)') "- TS duration:", dt
    write(6, '()')

end subroutine read_input


! Allocate arrays and zero observable arrays
subroutine allocate_arrays()

    use variables
    implicit none
    double precision :: evec(S,S), eval(S)

    ! LAPACK work array
    allocate( work(1) )
    call dsyev("V", "U", S, evec, S, eval, work, -1, info)
    lenwork = int(work(1))
    deallocate(work)
    allocate( work(lenwork) )

    ! Standard arrays
    allocate( c(F) )
    allocate( xn(F) )
    allocate( pn(F) )
    allocate( G0(F) )
    allocate( omega(F) )
    allocate( rho(tsteps+1,S,S) )

    ! Zeroing
    rho(:,:,:) = 0.d0

end subroutine allocate_arrays


! Deallocate arrays
subroutine deallocate_arrays()

    use variables
    implicit none

    ! LAPACK work array
    deallocate( work )

    ! Standard arrays
    deallocate( c )
    deallocate( xn )
    deallocate( pn )
    deallocate( G0 )
    deallocate( rho )
    deallocate( omega )

end subroutine deallocate_arrays


! Set bath frequencies and coupling coefficients and output both
subroutine bath_properties()

    use variables
    implicit none
    integer :: i,j
    double precision :: d

    ! Set position of two-level system in the cavity
    d = L/2.d0
    
    open(11, file="freq.out", status="unknown", action="write")
    do i = 1,F
        if ( d == L/2.d0 ) then ! If TLS @ of cavity, use cancellation of modes
            omega(i) = pi * sol * dble(2 * i - 1) / L 
            c(i) = mu * omega(i) * dsqrt(2.d0/eps0/L) * (-1.d0)**(i+1)
        else
            omega(i) = pi * sol * dble(i) / L
            c(i) = mu * omega(i) * dsqrt(2.d0/eps0/L) * sin(omega(i)/sol * d)
        end if
        write(11, *) omega(i), c(i)
    end do
    close(11)

end subroutine bath_properties


! Sample bath variables based on initial density matrix
subroutine sample_nuclear(tj)

    use variables
    implicit none
    integer :: i
    integer, intent(in) :: tj
    double precision :: r1, r2, xn_stdev, pn_stdev

    if ( tj == 1 ) then
        write(6, '(2x,"- Field initial conditions from vacuum state")')
    end if

    do i = 1,F
        xn_stdev = dsqrt( 1.d0 / (2.d0 * omega(i)) ) 
        pn_stdev = dsqrt( omega(i) / 2.d0 )
        call gauss(r1, r2)
        xn(i) = xn_stdev * r1
        pn(i) = pn_stdev * r2
    end do

end subroutine sample_nuclear


! Sample mapping variables based on initial density matrix
subroutine sample_electronic(tj)

    use variables
    implicit none
    integer :: i
    real(4) :: ran(1)
    integer, intent(in) :: tj

    if ( tj == 1 ) then
        select case ( init_state )
            case ( 11 )
            write(6, '(2x,"- Subsystem initial conditions: |1><1|")')
            
            case ( 22 )
            write(6, '(2x,"- Subsystem initial conditions: |2><2|")')
        end select
    end if

    call RANLUX(ran, 1)
    if ( init_state == 11 ) then
            XE(1) = dsqrt(2.d0) * cos(ran(1)*2.d0*pi)
            PE(1) = dsqrt(2.d0) * sin(ran(1)*2.d0*pi)
            XE(2) = 0.d0
            PE(2) = 0.d0
    else if ( init_state == 22 ) then
            XE(2) = dsqrt(2.d0) * cos(ran(1)*2.d0*pi)
            PE(2) = dsqrt(2.d0) * sin(ran(1)*2.d0*pi)
            XE(1) = 0.d0
            PE(1) = 0.d0
    end if

end subroutine sample_electronic


! Calucate observables
subroutine calculate_obs(ts)

    use variables
    implicit none
    integer :: i
    integer, intent(in) :: ts
    double precision :: sigx, sigy

    ! Populations
    do i = 1, S
        rho(ts,i,i) = rho(ts,i,i) + 0.5d0 * (XE(i)**2 + PE(i)**2)
    end do

    ! Coherences
    sigx = XE(1)*XE(2) + PE(1)*PE(2)
    sigy = XE(1)*PE(2) - PE(1)*XE(2)
    rho(ts,1,2) = rho(ts,1,2) + (sigx - eye*sigy)
    rho(ts,2,1) = rho(ts,2,1) + (sigx + eye*sigy)


end subroutine calculate_obs


! Average and output observables
subroutine average_obs()

    use variables
    implicit none
    integer :: i
    character(len=120) :: fmt

    write(6, '("Averaging and outputting observables:")')

    open(11, file="rho.out", status="unknown", action="write")
    write(fmt,*) "(f10.4,6(2x,ES13.5))"
    rho(:,:,:) = rho(:,:,:) / dble(ntraj)
    do i = 1,tsteps+1
        write(11,fmt) dble(i-1)*dt, &
                      real(rho(i,1,1)), real(rho(i,1,2)), aimag(rho(i,1,2)), &
                      real(rho(i,2,1)), aimag(rho(i,2,1)), real(rho(i,2,2))
    end do
    close(11)

end subroutine average_obs


! Makes a single trajectory step using velocity verlet
subroutine step_vverlet()

    use variables
    implicit none
    integer :: i,j
    double precision :: hdt, qdt

    hdt = 0.5d0*dt
    qdt = 0.5d0*hdt

    call potential_force()

    ! Half step in nuclear momenta
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        pn(i) = pn(i) - hdt * (c(i)*(XE(1)*XE(2) + PE(1)*PE(2)))
    end do

    ! Half step in mapping momenta
    PE = PE - hdt * matmul(V, XE)

    ! Full step in nuclear positions
    do i = 1,F
        xn(i) = xn(i) + dt * pn(i)
    end do

    ! Full step in mapping positions
    call potential_force()
    XE = XE + dt * matmul(V, PE)

    ! Half step in mapping momenta
    call potential_force()
    PE = PE - hdt * matmul(V, XE)

    ! Half step in nuclear momenta
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        pn(i) = pn(i) - hdt * (c(i)*(XE(1)*XE(2) + PE(1)*PE(2)))
    end do

end subroutine step_vverlet


! Makes a single trajectory step using diagonalisation of the potential matrix
subroutine step_diag

    use variables
    implicit none
    integer :: i,j
    double complex :: propagator(S) 
    double precision :: hdt, qdt, eval(S), evec(S,S), XEU(S), PEU(S)

    hdt = 0.5d0*dt
    qdt = 0.5d0*hdt

    call potential_force()

    ! Half step in mapping variables
    evec(:,:) = V(:,:)
    call dsyev("V", "U", S, evec, S, eval, work, lenwork, info)
    XEU = matmul(XE, evec)
    PEU = matmul(PE, evec)
    propagator = exp(-eye * hdt * eval) * dcmplx(XEU, PEU)
    XE = matmul(evec, real(propagator))
    PE = matmul(evec, aimag(propagator))

    ! Half step in nuclear momenta
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        pn(i) = pn(i) - hdt * (c(i)*(XE(1)*XE(2) + PE(1)*PE(2)))
    end do

    ! Full step in nuclear positions
    do i = 1,F
        xn(i) = xn(i) + dt * pn(i)
    end do

    ! Half step in nuclear momenta
    call potential_force()
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        pn(i) = pn(i) - hdt * (c(i)*(XE(1)*XE(2) + PE(1)*PE(2)))
    end do

    ! Half step in mapping variables
    evec(:,:) = V(:,:)
    call dsyev("V", "U", S, evec, S, eval, work, lenwork, info)
    XEU = matmul(XE, evec)
    PEU = matmul(PE, evec)
    propagator = exp(-eye * hdt * eval) * dcmplx(XEU, PEU)
    XE = matmul(evec, real(propagator))
    PE = matmul(evec, aimag(propagator))

end subroutine step_diag


! Calculates potential energy and force
subroutine potential_force()

    use variables
    implicit none
    integer :: i,j
    double precision :: tr, trace

    ! State-independent potential and force
    V0 = 0.d0
    do i = 1,F
        V0 = V0 + 0.5d0 * omega(i)**2 * xn(i)**2
        G0(i) = omega(i)**2 * xn(i)
    end do

    ! Potential energy matrix and force tensor
    V(:,:) = 0.d0
    V(1,1) = eps(1)
    V(2,2) = eps(2)
    do i = 1,F
        V(1,2) = V(1,2) + c(i) * xn(i)
        V(2,1) = V(2,1) + c(i) * xn(i)
    end do

    ! Shift trace of v. note that g is already traceless
    tr = trace(V,S)/dble(S)
    V0 = V0 + tr
    do i = 1,S
        V(i,i) = V(i,i) - tr
    end do

end subroutine potential_force


! Returns two numbers sampled from the standard normal distribution
! which are obtained via the Box-Muller transform of a RANLUX uniform call
subroutine gauss(r1, r2)

    use variables
    implicit none
    real(4) :: yield(2)
    double precision, intent(out) :: r1, r2

    call RANLUX(yield, 2)

    r1 = dsqrt(-2.d0 * log(yield(1))) * cos(2.d0 * pi * yield(2))
    r2 = dsqrt(-2.d0 * log(yield(1))) * sin(2.d0 * pi * yield(2))

end subroutine gauss


! Calculates the trace of a square matrix
double precision function trace(A,D)

    implicit none
    integer :: i
    integer, intent(in) :: D
    double precision, intent(in) :: A(D,D)
    
    trace = 0.d0

    do i = 1, D
        trace = trace + A(i,i)
    end do

end function trace
