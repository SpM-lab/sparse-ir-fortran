module sparse_ir_io
    use sparse_ir
    implicit none

    contains

    ! Read sampling points, basis functions
    function read_ir(unit, beta) result(obj)
        integer, intent (in) :: unit
        double precision, intent (in) :: beta

        type(IR) :: obj
        integer :: version
        character(len=100) :: tmp_str

        read(unit,*) tmp_str, version
        if (version == 1) then
            obj = read_v1(unit, beta)
        else
            write(*, *) "Invalid version number", version
            stop "Stopping..."
        end if
    end function


    ! Read sampling points, basis functions (version 1)
    function read_v1(unit, beta) result(obj)
        integer, intent (in) :: unit
        double precision, intent (in) :: beta

        type(IR) :: obj

        character(len=100) :: tmp_str
        integer :: i, l, t, n
        double precision :: rtmp, rtmp2

        double precision :: lambda, eps
        double precision, parameter :: rtol = 1e-20
        integer :: size, ntau, nfreq_f, nfreq_b, nomega
        double precision, allocatable :: s(:), tau(:), omega(:)
        integer, allocatable :: freq_f(:), freq_b(:)
        complex(kind(0d0)), allocatable :: u(:, :)
        complex(kind(0d0)), allocatable :: uhat_f(:, :), uhat_b(:, :)
        complex(kind(0d0)), allocatable :: v(:, :), dlr(:, :)

        read(unit,*) tmp_str, lambda
        read(unit,*) tmp_str, eps

        ! Singular values
        read(unit,*)
        read(unit,*) size
        allocate(s(size))
        do i=1, size
            read(unit, *) s(i)
        end do

        ! Sampling times
        read(unit,*)
        read(unit,*) ntau
        allocate(tau(ntau))
        do i=1, ntau
            read(unit, *) tau(i)
        end do

        ! Basis functions on sampling times
        read(unit,*)
        allocate(u(ntau, size))
        do l = 1, size
            do t = 1, ntau
                read(unit, *) rtmp
                u(t, l) = rtmp
            end do
        end do

        ! Sampling frequencies (F)
        read(unit,*)
        read(unit,*) nfreq_f
        allocate(freq_f(nfreq_f))
        do i=1, nfreq_f
            read(unit, *) freq_f(i)
        end do

        read(unit,*)
        allocate(uhat_f(nfreq_f, size))
        do l = 1, size
            do n = 1, nfreq_f
                read(unit, *) rtmp, rtmp2
                uhat_f(n, l) = cmplx(rtmp, rtmp2, kind(0d0))
            end do
        end do

        ! Sampling frequencies (B)
        read(unit,*)
        read(unit,*) nfreq_b
        allocate(freq_b(nfreq_b))
        do i=1, nfreq_b
            read(unit, *) freq_b(i)
        end do

        read(unit,*)
        allocate(uhat_b(nfreq_b, size))
        do l = 1, size
            do n = 1, nfreq_b
                read(unit, *) rtmp, rtmp2
                uhat_b(n, l) = cmplx(rtmp, rtmp2, kind(0d0))
            end do
        end do

        ! Sampling poles on real frequencies
        read(unit,*)
        read(unit,*) nomega
        allocate(omega(nomega))
        do i=1, nomega
            read(unit, *) omega(i)
        end do

        ! Right singular functions on sampling poles
        read(unit,*)
        allocate(v(nomega, size))
        allocate(dlr(nomega, size))
        do l = 1, size
            do i = 1, nomega
                read(unit, *) rtmp
                v(i, l) = rtmp
                dlr(i, l) = - s(l) * v(i, l)
            end do
        end do

        call init_ir(obj, beta, lambda, eps, s, tau, freq_f, freq_b, u, uhat_f, uhat_b, omega, v, dlr, 1d-16)

        deallocate(u, uhat_f, uhat_b, v, dlr)
    end function

end module
