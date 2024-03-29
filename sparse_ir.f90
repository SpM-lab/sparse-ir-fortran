module sparse_ir
    implicit none

    ! Matrix decomposed in SVD for fitting
    type DecomposedMatrix
        complex(kind(0d0)), allocatable :: a(:, :) ! Original matrix
        double precision, allocatable :: inv_s_dl(:) ! Inverse of dimensionless singular values
        double precision, allocatable :: inv_s(:) ! Inverse of singular values
        complex(kind(0d0)), allocatable :: ut(:, :), v(:, :)
        integer :: m, n, ns
    end type

    ! Sampling points, basis functions
    type IR
        integer :: size, ntau, nfreq_f, nfreq_b, nomega
        double precision :: beta, lambda, wmax, eps, eps_svd
        double precision, allocatable :: s(:), tau(:), x(:)
        double precision, allocatable :: omega(:), y(:)
        integer, allocatable :: freq_f(:), freq_b(:)
        complex(kind(0d0)), allocatable :: u_data(:,:)
        complex(kind(0d0)), allocatable :: uhat_f_data(:,:), uhat_b_data(:,:)
        complex(kind(0d0)), allocatable :: v_data(:,:), dlr_data(:,:)
        type(DecomposedMatrix) :: u
        type(DecomposedMatrix) :: uhat_f, uhat_b
        type(DecomposedMatrix) :: dlr
    end type

    contains

    subroutine init_ir(obj, beta, lambda, eps, s, x, freq_f, freq_b, u, uhat_f, uhat_b, y, v, dlr, eps_svd)
        type(IR), intent(inout) :: obj
        double precision, intent(in) :: beta, lambda, eps, s(:), x(:), y(:), eps_svd
        complex(kind(0d0)), intent(in) :: u(:,:), uhat_f(:, :), uhat_b(:, :), v(:, :), dlr(:, :)
        integer, intent(in) :: freq_f(:), freq_b(:)

        if (allocated(obj%x)) then
            stop 'IR%x is already allocated. You should call finalize_ir before recalling init_ir.'
        end if

        obj%size = size(s)
        obj%ntau = size(x)
        obj%nfreq_f = size(freq_f)
        obj%nfreq_b = size(freq_b)
        obj%lambda = lambda
        obj%eps = eps
        obj%eps_svd = eps_svd
        obj%nomega = size(y)

        allocate(obj%x(obj%ntau))
        obj%x = x

        allocate(obj%tau(obj%ntau))

        allocate(obj%s(obj%size))
        obj%s = sqrt(5.0d-1*obj%lambda) * s

        allocate(obj%freq_f(obj%nfreq_f))
        obj%freq_f = freq_f

        allocate(obj%freq_b(obj%nfreq_b))
        obj%freq_b = freq_b

        allocate(obj%u_data(obj%ntau, obj%size))
        obj%u_data = u

        allocate(obj%uhat_f_data(obj%nfreq_f, obj%size))
        obj%uhat_f_data = uhat_f

        allocate(obj%uhat_b_data(obj%nfreq_b, obj%size))
        obj%uhat_b_data = uhat_b

        allocate(obj%y(obj%nomega))
        obj%y = y

        allocate(obj%omega(obj%nomega))

        allocate(obj%v_data(obj%nomega, obj%size))
        obj%v_data = v

        allocate(obj%dlr_data(obj%size, obj%nomega))
        obj%dlr_data = transpose(dlr)

        obj%u = decompose(obj%u_data, obj%eps_svd)
        obj%uhat_f = decompose(obj%uhat_f_data, obj%eps_svd)
        obj%uhat_b = decompose(obj%uhat_b_data, obj%eps_svd)
        obj%dlr = decompose2(obj%dlr_data, obj%eps_svd)

        ! Here we define basis sets for the input value of beta. 
        call set_beta(obj, beta)
    end subroutine

    subroutine set_beta(obj, beta)
        type(IR), intent(inout) :: obj
        double precision, intent(in) :: beta

        obj%beta = beta
        obj%wmax = obj%lambda / beta

        obj%tau = 5.0d-1 * beta * (obj%x + 1.d0)
        obj%omega = obj%y * obj%wmax

        obj%u%a(:, :) = sqrt(2.0d0/beta)*obj%u_data(:, :)
        obj%uhat_f%a(:, :) = sqrt(beta) * obj%uhat_f_data(:, :)
        obj%uhat_b%a(:, :) = sqrt(beta) * obj%uhat_b_data(:, :)
        obj%dlr%a(:, :) = sqrt(5.0d-1*beta)*obj%dlr_data(:, :)

        obj%u%inv_s(:) = sqrt(5.0d-1*beta) * obj%u%inv_s_dl(:)
        obj%uhat_f%inv_s(:) = (1.0d0 / sqrt(beta)) * obj%uhat_f%inv_s_dl(:)
        obj%uhat_b%inv_s(:) = (1.0d0 / sqrt(beta)) * obj%uhat_b%inv_s_dl(:)
        obj%dlr%inv_s(:) = sqrt(2.0d0 / beta) * obj%dlr%inv_s_dl(:)

    end subroutine

    subroutine finalize_ir(obj)
        type(IR) :: obj
    
        if (allocated(obj%x)) deallocate(obj%x)
        if (allocated(obj%tau)) deallocate(obj%tau)
        if (allocated(obj%s)) deallocate(obj%s)
        if (allocated(obj%freq_f)) deallocate(obj%freq_f)
        if (allocated(obj%freq_b)) deallocate(obj%freq_b)
        if (allocated(obj%u_data)) deallocate(obj%u_data)
        if (allocated(obj%uhat_f_data)) deallocate(obj%uhat_f_data)
        if (allocated(obj%uhat_b_data)) deallocate(obj%uhat_b_data)
        if (allocated(obj%y)) deallocate(obj%y)
        if (allocated(obj%omega)) deallocate(obj%omega)
        if (allocated(obj%v_data)) deallocate(obj%v_data)
        if (allocated(obj%dlr_data)) deallocate(obj%dlr_data)
    
        call finalize_dmat(obj%u)
        call finalize_dmat(obj%uhat_f)
        call finalize_dmat(obj%uhat_b)
        call finalize_dmat(obj%dlr)
    end subroutine

    subroutine finalize_dmat(dmat)
        type(DecomposedMatrix) :: dmat
    
        if (allocated(dmat%a)) deallocate(dmat%a)
        if (allocated(dmat%inv_s)) deallocate(dmat%inv_s)
        if (allocated(dmat%inv_s)) deallocate(dmat%inv_s_dl)
        if (allocated(dmat%ut)) deallocate(dmat%ut)
        if (allocated(dmat%v)) deallocate(dmat%v)
    end subroutine

    ! SVD of matrix a. Singular values smaller than eps * the largest one are dropped.
    function decompose(a, eps) result(dmat)
        complex(kind(0d0)), intent(in) :: a(:, :)
        double precision, intent(in) :: eps

        integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns
        complex(kind(0d0)), allocatable :: a_copy(:, :), u(:, :), &
            vt(:, :), work(:)
        double precision, allocatable :: rwork(:), s(:)
        integer, allocatable :: iwork(:)
        type(DecomposedMatrix)::dmat

        if (allocated(dmat%a)) then
            stop 'DMAT%a is already allocated. You should call finalize_dmat before recalling decompose.'
        end if

        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        lda = m
        ldu = m
        ldvt = n
        lwork = mn*mn + 3*mn

        allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork((5*mn+7)*mn), iwork(8*mn))

        a_copy(1:m, 1:n) = a(1:m, 1:n)
        call zgesdd('S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)

        if (info /= 0) then
            stop 'Failure in ZGESDD.'
        end if

        ! Number of relevant singular values s(i)/s(1) >= eps
        ns = 0
        do i = 1, mn
            if (s(i)/s(1) < eps) then
                exit
            end if
            ns = ns + 1
        end do

        allocate(dmat%a(m, n))
        allocate(dmat%inv_s_dl(ns))
        allocate(dmat%inv_s(ns))
        allocate(dmat%ut(ns, m))
        allocate(dmat%v(n, ns))

        ! dmat%a temporarily stores the same data of input a
        dmat%a = a
        dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
        ! inv_s temporarily stores the same data of inv_s_dl
        dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
        dmat%ut(1:ns, 1:m) = conjg(transpose(u(1:m, 1:ns)))
        dmat%v(1:n, 1:ns) = conjg(transpose(vt(1:ns, 1:n)))
        dmat%m = size(a, 1)
        dmat%n = size(a, 2)
        dmat%ns = ns

        deallocate(work, a_copy, s, u, vt, rwork, iwork)
    end function

    ! SVD of matrix a. Singular values smaller than eps * the largest one are dropped.
    function decompose2(a, eps) result(dmat)
        complex(kind(0d0)), intent(in) :: a(:, :)
        double precision, intent(in) :: eps
    
        integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns
        complex(kind(0d0)), allocatable :: a_copy(:, :), u(:, :), &
            vt(:, :), work(:)
        double precision, allocatable :: rwork(:), s(:)
        type(DecomposedMatrix)::dmat

        if (allocated(dmat%a)) then
            stop 'DMAT%a is already allocated. You should call finalize_dmat before recalling decompose2.'
        end if
    
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        lda = m
        ldu = m
        ldvt = n
        lwork = 2*mn + m + n
    
        allocate(work(lwork), a_copy(m,n), s(m), u(ldu,m), vt(ldvt,n), rwork(5*n))
    
        a_copy(1:m, 1:n) = a(1:m, 1:n)
        call zgesvd('S', 'S', m, n, a_copy, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
    
        if (info /= 0) then
            stop 'Failure in ZGESVD.'
        end if
    
        ! Number of relevant singular values s(i)/s(1) >= eps
        ns = 0
        do i = 1, mn
            if (s(i)/s(1) < eps) then
                exit
            end if
            ns = ns + 1
        end do
    
        allocate(dmat%a(m, n))
        allocate(dmat%inv_s_dl(ns))
        allocate(dmat%inv_s(ns))
        allocate(dmat%ut(ns, m))
        allocate(dmat%v(n, ns))

        ! dmat%a temporarily stores the same data of input a
        dmat%a = a
        dmat%inv_s_dl(1:ns) = 1.0D0 / s(1:ns)
        ! inv_s temporarily stores the same data of inv_s_dl
        dmat%inv_s(1:ns) = dmat%inv_s_dl(1:ns)
        dmat%ut(1:ns, 1:m) = conjg(transpose(u(1:m, 1:ns)))
        dmat%v(1:n, 1:ns) = conjg(transpose(vt(1:ns, 1:n)))
        dmat%m = size(a, 1)
        dmat%n = size(a, 2)
        dmat%ns = ns
    
        deallocate(work, a_copy, s, u, vt, rwork)
    end function

    subroutine fit_matsubara_f(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        call fit_impl(arr, obj%uhat_f, res)
    end subroutine

    subroutine fit_matsubara_b(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        call fit_impl(arr, obj%uhat_b, res)
    end subroutine

    subroutine fit_tau(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        call fit_impl(arr, obj%u, res)
    end subroutine

    subroutine evaluate_tau(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), PARAMETER :: cone  = (1.0d0, 0.0d0)
        complex(kind(0d0)), PARAMETER :: czero  = (0.0d0, 0.0d0)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%u%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%u%m) stop 'wrong number of columns of output array.'
        !
        CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
                   l2, obj%u%a, m, czero, res(:, :), l2)
        ! res = matmul(arr, transpose(obj%u%a))
    end subroutine

    subroutine evaluate_matsubara_f(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), PARAMETER :: cone  = (1.0d0, 0.0d0)
        complex(kind(0d0)), PARAMETER :: czero  = (0.0d0, 0.0d0)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%uhat_f%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%uhat_f%m) stop 'wrong number of columns of output array.'
        !
        CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
                   l2, obj%uhat_f%a, m, czero, res(:, :), l2)
        ! res = matmul(arr, transpose(obj%uhat_f%a))
    end subroutine

    subroutine evaluate_matsubara_b(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), PARAMETER :: cone  = (1.0d0, 0.0d0)
        complex(kind(0d0)), PARAMETER :: czero  = (0.0d0, 0.0d0)
        integer :: m, n, l1, l2
        !
        l1 = size(arr, 1)
        n = size(arr, 2)
        l2 = size(res, 1)
        m = size(res, 2)
        IF (l1 .NE. l2) stop 'wrong number of rows of input array.'
        IF (n .NE. obj%uhat_b%n) stop 'wrong number of columns of input array.'
        IF (m .NE. obj%uhat_b%m) stop 'wrong number of columns of output array.'
        !
        CALL ZGEMM('n', 't', l1, m, n, cone, arr(:,:), &
                   l2, obj%uhat_b%a, m, czero, res(:, :), l2)
        ! res = matmul(arr, transpose(obj%uhat_b%a))
    end subroutine

    ! Implementation of fit
    subroutine fit_impl(arr, mat, res)
        complex(kind(0d0)), intent (in) :: arr(:, :)
        type(DecomposedMatrix), intent(in) :: mat
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), PARAMETER :: cone  = (1.0d0, 0.0d0)
        complex(kind(0d0)), PARAMETER :: czero  = (0.0d0, 0.0d0)
        complex(kind(0d0)), allocatable :: ut_arr(:, :)

        integer :: nb, m, n, ns, i, j

        ! ut(ns, m)
        ! v(n, ns)
        ! arr(nb, m)
        ! mat(m, n)
        ! ut_arr(ns, nb)
        ! res(nb, n)
        nb = size(arr, 1)
        ns = mat%ns
        m = mat%m
        n = mat%n
        allocate(ut_arr(ns, nb))

        if (size(res, 1) /= nb .or. size(res, 2) /= n) then
            stop 'Invalid size of output array'
        end if

        !ut(ns, m) * arr(nb, m) -> ut_arr(ns, nb)
        ut_arr = 0.0
        call zgemm("n", "t", ns, nb, m, cone, mat%ut, ns, arr, nb, czero, ut_arr, ns)
        do j = 1, ns
            do i = 1, nb
                ut_arr(j, i) = ut_arr(j, i) * mat%inv_s(j)
            end do
        end do

        ! ut_arr(ns, nb) * v(n, ns) -> (nb, n)
        res = 0.0
        call zgemm("t", "t", nb, n, ns, cone, ut_arr, ns, mat%v, n, czero, res, nb)

        deallocate(ut_arr)
    end subroutine

    subroutine to_dlr(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        call fit_impl(arr, obj%dlr, res)
    end subroutine
      
    subroutine evaluate_tau_from_dlr(obj, tau, arr, res)
        type(IR), intent(in) :: obj
        double precision, intent(in) :: tau(:)
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        double precision :: kernel
        integer :: ntau, nt, p, t, l1, l2
        !
        ntau = size(tau)
        nt = size(res, 2)
        IF (ntau .NE. nt) stop 'wrong number of columns of output array.'
        l1 = size(arr, 1)
        l2 = size(res, 1)
        IF (l1 .NE. l2) stop 'wrong number of rows of output array.'
        !
        do t = 1, ntau
          if ((tau(t) < 0.0d0) .OR. (obj%beta < tau(t))) then
            stop 'tau must be in [0, beta].'
          end if
        end do
        !
        res(:, :) = (0d0, 0d0)
        !
        do t = 1, ntau
            do p = 1, obj%nomega
                if (obj%omega(p) < 0) then
                    if ( (obj%beta - tau(t)) * obj%omega(p) < -100.d0) then
                        kernel = 0.0d0
                    else if (obj%beta * obj%omega(p) < -30.d0) then
                        kernel = exp( (obj%beta - tau(t)) * obj%omega(p))
                    else
                        kernel = exp( (obj%beta - tau(t)) * obj%omega(p)) / (exp(obj%beta * obj%omega(p)) + 1.0d0) 
                    end if
                else
                    if (tau(t) * obj%omega(p) > 100.d0) then
                        kernel = 0.0d0
                    else if (obj%beta * obj%omega(p) > 30.d0) then
                        kernel = exp(- tau(t) * obj%omega(p))
                    else
                        kernel = exp(- tau(t) * obj%omega(p)) / (1.0d0 + exp(- obj%beta * obj%omega(p))) 
                    end if
                end if
                res(:, t) = res(:, t) - kernel * arr(:, p)
            end do
        end do
        !
    end subroutine
      
    subroutine evaluate_matsubara_f_from_dlr(obj, freq, arr, res)
        type(IR), intent(in) :: obj
        integer, intent(in) :: freq(:)
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), PARAMETER :: ci  = (0.0d0, 1.0d0)
        double precision :: PI
        complex(kind(0d0)) :: cfreq
        complex(kind(0d0)) :: kernel
        integer :: nfreq, nf, p, n, l1, l2
        !
        PI =4.D0*DATAN(1.D0)
        !
        nfreq = size(freq)
        nf = size(res, 2)
        IF (nfreq .NE. nf) stop 'wrong number of columns of output array.'
        l1 = size(arr, 1)
        l2 = size(res, 1)
        IF (l1 .NE. l2) stop 'wrong number of rows of output array.'
        !
        do n = 1, nfreq
          if (MOD(freq(n), 2) == 0) stop 'one of input integers is not odd.'
        end do
        !
        res(:, :) = (0d0, 0d0)
        !
        do n = 1, nfreq
            cfreq = ci * pi * REAL(freq(n), kind(0d0)) / obj%beta
            do p = 1, obj%nomega
                kernel = 1.0d0 / (cfreq - obj%omega(p))
                res(:, n) = res(:, n) + kernel * arr(:, p)
            end do
        end do
        !
    end subroutine

    subroutine evaluate_matsubara_b_from_dlr(obj, freq, arr, res)
        type(IR), intent(in) :: obj
        integer, intent(in) :: freq(:)
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        complex(kind(0d0)), PARAMETER :: ci  = (0.0d0, 1.0d0)
        double precision :: PI
        complex(kind(0d0)) :: cfreq
        complex(kind(0d0)) :: kernel
        integer :: nfreq, nf, p, n, l1, l2
        !
        PI =4.D0*DATAN(1.D0)
        !
        nfreq = size(freq)
        nf = size(res, 2)
        IF (nfreq .NE. nf) stop 'wrong number of columns of output array.'
        l1 = size(arr, 1)
        l2 = size(res, 1)
        IF (l1 .NE. l2) stop 'wrong number of rows of output array.'
        !
        do n = 1, nfreq
          if (MOD(freq(n), 2) .ne. 0) stop 'one of input integers is not even.'
        end do
        !
        res(:, :) = (0d0, 0d0)
        !
        do n = 1, nfreq
            cfreq = ci * pi * REAL(freq(n), kind(0d0)) / obj%beta
            do p = 1, obj%nomega
                kernel = TANH(5.0d-1 * obj%beta * obj%omega(p)) / (cfreq - obj%omega(p))
                res(:, n) = res(:, n) + kernel * arr(:, p)
            end do
        end do
        !
    end subroutine

end module
