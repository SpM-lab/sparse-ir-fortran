program main
    use sparse_ir
    use sparse_ir_io
    use sparse_ir_preset
    implicit none

    call test_fit()
    call test_under_fitting()
    call test_over_fitting()
    call test_fermion(.true.)
    call test_boson(.true.)
    call test_fermion_dlr(.true.)
    call test_boson_dlr(.true.)
    call test_fermion(.false.)
    call test_boson(.false.)
    call test_fermion_dlr(.false.)
    call test_boson_dlr(.false.)

    contains

    subroutine test_fit()
        complex(kind(0d0)) :: a(2, 2), y(1, 2), x(1, 2), y_reconst(1, 2)
        type(DecomposedMatrix) :: dm
        a(1, 1) = 2.d0
        a(2, 1) = 0.1d0
        a(1, 2) = 0.d0
        a(2, 2) = 1.d0

        dm = decompose(a, 1d-20)

        y(1, 1) = 0.2d0
        y(1, 2) = 0.1d0

        call fit_impl(y, dm, x)

        y_reconst = transpose(matmul(dm%a, transpose(x)))
        if (maxval(abs(y - y_reconst)) > 1e-12) then
            stop "y and y_reconst do not match!"
        end if
        !write(*, *) y
        !write(*, *) y_reconst
        call finalize_dmat(dm)
    end subroutine

    subroutine test_over_fitting()
        integer, parameter :: n=1, m=2
        complex(kind(0d0)) :: a(n, m), y(1, n), x(1, m), y_reconst(1, n)
        type(DecomposedMatrix) :: dm
        a(1, 1) = 2.d0
        a(1, 2) = 1.d0

        dm = decompose(a, 1d-10)

        y(1, 1) = 0.2d0

        call fit_impl(y, dm, x)

        y_reconst = transpose(matmul(dm%a, transpose(x)))
        if (maxval(abs(y - y_reconst)) > 1e-12) then
            stop "y and y_reconst do not match!"
        end if
        call finalize_dmat(dm)
    end subroutine


    subroutine test_under_fitting()
        integer, parameter :: n=2, m=1
        complex(kind(0d0)) :: a(n, m), y(1, n), x(1, m), y_reconst(1, n)
        type(DecomposedMatrix) :: dm
        a(1, 1) = 1.2d0
        a(2, 1) = 1.d0

        dm = decompose(a, 1d-10)

        y(1, 1) = 1.2d0
        y(1, 2) = 1.0d0

        call fit_impl(y, dm, x)
        !write(*,*) "inv_s: ", dm%inv_s
        !write(*,*) "x: ", x

        ! transpose(x): (m, 1)
        ! matmul(dm%a, transpose(x)): (n, m) * (m, 1) = (n, 1)
        y_reconst = transpose(matmul(dm%a, transpose(x)))
        !write(*,*) "y: ", y_reconst
        if (maxval(abs(y - y_reconst)) > 1e-12) then
            stop "y and y_reconst do not match!"
        end if
        call finalize_dmat(dm)
    end subroutine


    ! fermion
    subroutine test_fermion(preset)
        logical, intent(in) :: preset
        type(IR) :: ir_obj
        integer, parameter :: ndigit = 10, nlambda = 4
        double precision, parameter :: lambda = 1.d1 ** nlambda
        double precision, parameter :: wmax = 1.d0
        double precision :: PI

        double precision, parameter :: beta = lambda/wmax, omega0 = 1.d0/beta
        double precision, parameter :: eps = 1.d-1**ndigit

        complex(kind(0d0)),allocatable :: giv(:,:), gl_matsu(:, :), gl_tau(:, :), gtau(:, :), &
            gtau_reconst(:, :), giv_reconst(:, :)
        integer n, t

        PI =4.D0*DATAN(1.D0)

        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta)
            close(99)
        end if

        if (abs(ir_obj%beta - beta) > 1d-10) then
            stop "beta does not match"
        end if
        if (abs(ir_obj%wmax - wmax) > 1d-10) then
            stop "wmax does not match"
        end if

        ! With ω0 = 1/β,
        !   G(iv) = 1/(iv - ω0),
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        allocate(giv(1, ir_obj%nfreq_f))
        allocate(gtau(1, ir_obj%ntau))
        allocate(gl_matsu(1, ir_obj%size))
        allocate(gl_tau(1, ir_obj%size))
        allocate(gtau_reconst(1, ir_obj%ntau))
        allocate(giv_reconst(1, ir_obj%nfreq_f))

        ! From Matsubara
        do n = 1, ir_obj%nfreq_f
            giv(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_f(n)/beta, kind(0d0)) - omega0)
        end do
        call fit_matsubara_f(ir_obj, giv, gl_matsu)

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 + exp(-beta * omega0))
        end do
        call fit_tau(ir_obj, gtau, gl_tau)

        !do l = 1, ir_obj%size
            !write(*,*) real(gl_matsu(1,l)), real(gl_tau(1,l))
        !end do
        if (maxval(abs(gl_matsu - gl_tau)) > 1d2*eps) then
            stop "gl_matsu and gl_tau do not match!"
        end if

        call evaluate_matsubara_f(ir_obj, gl_matsu, giv_reconst)
        if (maxval(abs(giv - giv_reconst)) > 1d2*eps) then
            write(*, *) "C"
            stop "giv do not match!"
        end if

        call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
        if (maxval(abs(gtau - gtau_reconst)) > 1d2*eps) then
            stop "gtau do not match!"
        end if

        deallocate(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)

        call finalize_ir(ir_obj)
    end subroutine


    ! boson
    subroutine test_boson(preset)
        logical, intent(in) :: preset
        type(IR) :: ir_obj
        integer, parameter :: ndigit = 10, nlambda = 4
        double precision, parameter :: lambda = 1.d1 ** nlambda
        double precision, parameter :: wmax = 1.d0
        double precision :: PI

        double precision, parameter :: beta = lambda/wmax, omega0 = 1.d0/beta
        double precision, parameter :: eps = 1.d-1**ndigit

        complex(kind(0d0)),allocatable :: giv(:,:), gl_matsu(:, :), gl_tau(:, :), gtau(:, :), &
            gtau_reconst(:, :), giv_reconst(:, :)
        integer n, t

        PI=4.D0*DATAN(1.D0)

        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta)
            close(99)
        end if

        if (abs(ir_obj%beta - beta) > 1d-10) then
            stop "beta does not match"
        end if
        if (abs(ir_obj%wmax - wmax) > 1d-10) then
            stop "wmax does not match"
        end if

        ! With ω0 = 1/β,
        !   G(iv) = 1/(iv - ω0),
        !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
        allocate(giv(1, ir_obj%nfreq_b))
        allocate(gtau(1, ir_obj%ntau))
        allocate(gl_matsu(1, ir_obj%size))
        allocate(gl_tau(1, ir_obj%size))
        allocate(gtau_reconst(1, ir_obj%ntau))
        allocate(giv_reconst(1, ir_obj%nfreq_b))

        ! From Matsubara
        do n = 1, ir_obj%nfreq_b
            giv(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_b(n)/beta, kind(0d0)) - omega0)
        end do
        call fit_matsubara_b(ir_obj, giv, gl_matsu)

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 - exp(-beta * omega0))
        end do
        call fit_tau(ir_obj, gtau, gl_tau)

        !do l = 1, ir_obj%size
            !write(*, *) l, real(gl_tau(1, l)), real(gl_matsu(1, l))
        !end do

        if (maxval(abs(gl_matsu - gl_tau)) > 1d2*eps) then
            stop "gl_matsu and gl_tau do not match!"
        end if

        call evaluate_matsubara_b(ir_obj, gl_matsu, giv_reconst)
        if (maxval(abs(giv - giv_reconst)) > 1d2*eps) then
            write(*, *) "A"
            stop "gtau do not match!"
        end if

        call evaluate_tau(ir_obj, gl_tau, gtau_reconst)
        if (maxval(abs(gtau - gtau_reconst)) > 1d2*eps) then
            stop "gtau do not match!"
        end if

        deallocate(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        
        call finalize_ir(ir_obj)
    end subroutine

    ! fermion
    subroutine test_fermion_dlr(preset)
        logical, intent(in) :: preset
        type(IR) :: ir_obj
        integer, parameter :: ndigit = 10, nlambda = 4
        double precision, parameter :: lambda = 1.d1 ** nlambda
        double precision, parameter :: wmax = 1.d0
        double precision :: PI

        double precision, parameter :: beta = lambda/wmax, omega0 = 1.d0/beta
        double precision, parameter :: eps = 1.d-1**ndigit
        integer, parameter :: ntau_dlr = 200, nfreq_dlr = 200

        complex(kind(0d0)),allocatable :: giv_smpl(:,:), gl_matsu(:, :), gl_tau(:, :), gtau_smpl(:, :), &
            gtau_reconst(:, :), giv_reconst(:, :), g_dlr(:, :), giv_ref(:,:), gtau_ref(:,:)
        integer, allocatable :: freq(:) 
        double precision, allocatable :: tau(:)
        integer n, t

        PI =4.D0*DATAN(1.D0)

        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta)
            close(99)
        end if

        if (abs(ir_obj%beta - beta) > 1d-10) then
            stop "beta does not match"
        end if
        if (abs(ir_obj%wmax - wmax) > 1d-10) then
            stop "wmax does not match"
        end if

        ! With ω0 = 1/β,
        !   G(iv) = 1/(iv - ω0),
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        allocate(giv_smpl(1, ir_obj%nfreq_f))
        allocate(gtau_smpl(1, ir_obj%ntau))
        allocate(gl_matsu(1, ir_obj%size))
        allocate(gl_tau(1, ir_obj%size))
        allocate(g_dlr(1, ir_obj%nomega))
        allocate(giv_ref(1, nfreq_dlr))
        allocate(gtau_ref(1, ntau_dlr))
        allocate(gtau_reconst(1, ntau_dlr))
        allocate(giv_reconst(1, nfreq_dlr))
        allocate(freq(nfreq_dlr))
        allocate(tau(ntau_dlr))

        ! From Matsubara
        do n = 1, ir_obj%nfreq_f
            giv_smpl(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_f(n)/beta, kind(0d0)) - omega0)
        end do
        call fit_matsubara_f(ir_obj, giv_smpl, gl_matsu)

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau_smpl(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 + exp(-beta * omega0))
        end do
        call fit_tau(ir_obj, gtau_smpl, gl_tau)

        !do l = 1, ir_obj%size
            !write(*,*) real(gl_matsu(1,l)), real(gl_tau(1,l))
        !end do
        if (maxval(abs(gl_matsu - gl_tau)) > 1d2*eps) then
            stop "gl_matsu and gl_tau do not match!"
        end if

        call to_dlr (ir_obj, gl_matsu, g_dlr) 

        do n = 1, nfreq_dlr
            freq(n) = -nfreq_dlr + 2 * (n) - 1
            giv_ref(1, n) = 1.d0/(cmplx(0d0, PI*freq(n)/beta, kind(0d0)) - omega0)
        end do

        call evaluate_matsubara_f_from_dlr(ir_obj, freq, g_dlr, giv_reconst)
        if (maxval(abs(giv_ref - giv_reconst)) > 1d3*eps) then
            write(*, *) "AAAA"
            stop "giv do not match!"
        end if

        do n = 1, ntau_dlr
            tau(n) = beta * DBLE(n) / DBLE(ntau_dlr + 1)
            gtau_ref(1, n) = - exp(-tau(n) * omega0)/(1.d0 + exp(-beta * omega0))
        end do

        call evaluate_tau_from_dlr(ir_obj, tau, g_dlr, gtau_reconst)
        if (maxval(abs(gtau_ref - gtau_reconst)) > 1d3*eps) then
            stop "gtau do not match!"
        end if

        deallocate(giv_smpl, gtau_smpl, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        deallocate(giv_ref, gtau_ref, g_dlr, freq, tau)
        
        call finalize_ir(ir_obj)
    end subroutine

    ! boson
    subroutine test_boson_dlr(preset)
        logical, intent(in) :: preset
        type(IR) :: ir_obj
        integer, parameter :: ndigit = 10, nlambda = 4
        double precision, parameter :: lambda = 1.d1 ** nlambda
        double precision, parameter :: wmax = 1.d0
        double precision :: PI

        double precision, parameter :: beta = lambda/wmax, omega0 = 1.d0/beta
        double precision, parameter :: eps = 1.d-1**ndigit
        integer, parameter :: ntau_dlr = 200, nfreq_dlr = 200

        complex(kind(0d0)),allocatable :: giv_smpl(:,:), gl_matsu(:, :), gl_tau(:, :), gtau_smpl(:, :), &
            gtau_reconst(:, :), giv_reconst(:, :), g_dlr(:, :), giv_ref(:,:), gtau_ref(:,:)
        integer, allocatable :: freq(:) 
        double precision, allocatable :: tau(:)
        integer n, t

        PI =4.D0*DATAN(1.D0)

        if (preset) then
            ir_obj = mk_ir_preset(nlambda, ndigit, beta)
        else
            open(99, file='ir_nlambda4_ndigit10.dat', status='old')
            ir_obj = read_ir(99, beta)
            close(99)
        end if

        if (abs(ir_obj%beta - beta) > 1d-10) then
            stop "beta does not match"
        end if
        if (abs(ir_obj%wmax - wmax) > 1d-10) then
            stop "wmax does not match"
        end if

        ! With ω0 = 1/β,
        !   G(iv) = 1/(iv - ω0),
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        allocate(giv_smpl(1, ir_obj%nfreq_b))
        allocate(gtau_smpl(1, ir_obj%ntau))
        allocate(gl_matsu(1, ir_obj%size))
        allocate(gl_tau(1, ir_obj%size))
        allocate(g_dlr(1, ir_obj%nomega))
        allocate(giv_ref(1, nfreq_dlr))
        allocate(gtau_ref(1, ntau_dlr))
        allocate(gtau_reconst(1, ntau_dlr))
        allocate(giv_reconst(1, nfreq_dlr))
        allocate(freq(nfreq_dlr))
        allocate(tau(ntau_dlr))

        ! From Matsubara
        do n = 1, ir_obj%nfreq_b
            giv_smpl(1, n) = 1.d0/(cmplx(0d0, PI*ir_obj%freq_b(n)/beta, kind(0d0)) - omega0)
        end do
        call fit_matsubara_b(ir_obj, giv_smpl, gl_matsu)

        ! From tau
        !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
        do t = 1, ir_obj%ntau
            gtau_smpl(1, t) = - exp(-ir_obj%tau(t) * omega0)/(1.d0 - exp(-beta * omega0))
        end do
        call fit_tau(ir_obj, gtau_smpl, gl_tau)

        !do l = 1, ir_obj%size
            !write(*,*) real(gl_matsu(1,l)), real(gl_tau(1,l))
        !end do
        if (maxval(abs(gl_matsu - gl_tau)) > 1d2*eps) then
            stop "gl_matsu and gl_tau do not match!"
        end if

        call to_dlr (ir_obj, gl_matsu, g_dlr) 

        do n = 1, nfreq_dlr
            freq(n) = -nfreq_dlr + 2 * (n)
            giv_ref(1, n) = 1.d0/(cmplx(0d0, PI*freq(n)/beta, kind(0d0)) - omega0)
        end do

        call evaluate_matsubara_b_from_dlr(ir_obj, freq, g_dlr, giv_reconst)
        if (maxval(abs(giv_ref - giv_reconst)) > 1d3*eps) then
            write(*, *) "B"
            stop "giv do not match!"
        end if

        do n = 1, ntau_dlr
            tau(n) = beta * DBLE(n) / DBLE(ntau_dlr + 1)
            gtau_ref(1, n) = - exp(-tau(n) * omega0)/(1.d0 - exp(-beta * omega0))
        end do

        call evaluate_tau_from_dlr(ir_obj, tau, g_dlr, gtau_reconst)
        if (maxval(abs(gtau_ref - gtau_reconst)) > 1d3*eps) then
            stop "gtau do not match!"
        end if

        deallocate(giv_smpl, gtau_smpl, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
        deallocate(giv_ref, gtau_ref, g_dlr, freq, tau)
        
        call finalize_ir(ir_obj)
    end subroutine

end program
