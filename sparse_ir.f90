module sparse_ir
    implicit none

    ! Matrix decomposed in SVD for fitting
    type DecomposedMatrix
        complex(kind(0d0)), pointer :: a(:, :) ! Original matrix
        double precision, pointer :: inv_s(:) ! Inverse of singular values
        complex(kind(0d0)), pointer :: ut(:, :), v(:, :)
        integer :: m, n, ns
    end type

    ! Sampling points, basis functions
    type IR
        integer :: size, ntau, nfreq_f, nfreq_b
        double precision :: beta, wmax, eps
        double precision, pointer :: s(:), tau(:)
        integer, pointer :: freq_f(:), freq_b(:)
        type(DecomposedMatrix) :: u
        type(DecomposedMatrix) :: uhat_f, uhat_b
    end type

    contains

    subroutine init_ir(obj, beta, lambda, eps, s, x, freq_f, freq_b, u, uhat_f, uhat_b, eps_svd)
        type(IR), intent(inout) :: obj
        double precision, intent(in) :: beta, lambda, eps, s(:), x(:), eps_svd
        complex(kind(0d0)), intent(in) :: u(:,:), uhat_f(:, :), uhat_b(:, :)
        integer, intent(in) :: freq_f(:), freq_b(:)

        obj%size = size(s)
        obj%ntau = size(x)
        obj%nfreq_f = size(freq_f)
        obj%nfreq_b = size(freq_b)
        obj%beta = beta
        obj%wmax = lambda/beta
        obj%eps = eps

        allocate(obj%s(obj%size))
        obj%s = sqrt(0.5*lambda) * s

        allocate(obj%tau(obj%ntau))
        obj%tau = 0.5 * beta * (x + 1.d0)

        allocate(obj%freq_f(obj%nfreq_f))
        obj%freq_f = freq_f

        allocate(obj%freq_b(obj%nfreq_b))
        obj%freq_b = freq_b

        obj%u = decompose(sqrt(2/beta)*u, eps_svd)
        obj%uhat_f = decompose(sqrt(beta) * uhat_f, eps_svd)
        obj%uhat_b = decompose(sqrt(beta) * uhat_b, eps_svd)
    end subroutine

    ! SVD of matrix a. Singular values smaller than esp * the largest one are dropped.
    function decompose(a, eps) result(dmat)
        complex(kind(0d0)), intent(in) :: a(:, :)
        double precision, intent(in) :: eps

        integer :: i, info, lda, ldu, ldvt, lwork, m, n, mn, ns
        complex(kind(0d0)), allocatable :: a_copy(:, :), u(:, :), &
            vt(:, :), work(:)
        double precision, allocatable :: rwork(:), s(:)
        integer, allocatable :: iwork(:)
        type(DecomposedMatrix)::dmat

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
        allocate(dmat%inv_s(ns))
        allocate(dmat%ut(ns, m))
        allocate(dmat%v(n, ns))

        dmat%a = a
        dmat%inv_s(1:ns) = 1/s(1:ns)
        dmat%ut(1:ns, 1:m) = conjg(transpose(u(1:m, 1:ns)))
        dmat%v(1:n, 1:ns) = conjg(transpose(vt(1:ns, 1:n)))
        dmat%m = size(a, 1)
        dmat%n = size(a, 2)
        dmat%ns = ns

        deallocate(work, a_copy, s, u, vt, rwork, iwork)
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
        ! TODO: using zgemm
        res = matmul(arr, transpose(obj%u%a))
    end subroutine

    subroutine evaluate_matsubara_f(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        ! TODO: using zgemm
        res = matmul(arr, transpose(obj%uhat_f%a))
    end subroutine

    subroutine evaluate_matsubara_b(obj, arr, res)
        type(IR), intent(in) :: obj
        complex(kind(0d0)), intent (in) :: arr(:, :)
        complex(kind(0d0)), intent(out) :: res(:, :)
        ! TODO: using zgemm
        res = matmul(arr, transpose(obj%uhat_b%a))
    end subroutine

    ! Implementation of fit
    subroutine fit_impl(arr, mat, res)
        complex(kind(0d0)), intent (in) :: arr(:, :)
        type(DecomposedMatrix), intent(in) :: mat
        complex(kind(0d0)), intent(out) :: res(:, :)

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
        call zgemm("n", "t", ns, nb, m, 1.d0, mat%ut, ns, arr, nb, 0.d0, ut_arr, ns)
        do j = 1, ns
            do i = 1, nb
                ut_arr(j, i) = ut_arr(j, i) * mat%inv_s(j)
            end do
        end do

        ! ut_arr(ns, nb) * v(n, ns) -> (nb, n)
        res = 0.0
        call zgemm("t", "t", nb, n, ns, 1.d0, ut_arr, ns, mat%v, n, 0.d0, res, nb)

        deallocate(ut_arr)
    end subroutine


end module
