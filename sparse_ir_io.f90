  !----------------------------------------------------------------------
  MODULE sparse_ir_io
  !----------------------------------------------------------------------
  !
  USE sparse_ir
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  INTEGER, PARAMETER :: DP = KIND(0d0)
  !
  PUBLIC :: read_ir
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  FUNCTION read_ir(unit, beta, positive_only) RESULT(obj)
  !-----------------------------------------------------------------------
  !!
  !! This function calls read_v1 
  !! to read the file including the ir-basis objects.
  !!
  !
  INTEGER, INTENT(IN) :: unit
  !! Unit number
  REAL(KIND = DP), INTENT(IN) :: beta
  !! inverse temperature
  LOGICAL, INTENT(IN), OPTIONAL :: positive_only
  !! if true, take the Matsubara frequencies
  !! only from the positive region
  !
  TYPE(IR) :: obj
  !! to contain all the ir-basis objects from the file
  CHARACTER(LEN = 100) :: tmp_str
  !! dummy for characters
  INTEGER :: version
  !! version number
  !
  READ(unit,*) tmp_str, version
  IF (version == 1) then
    IF ((.NOT. PRESENT(positive_only))) then
      obj = read_v1(unit, beta)
    ELSE
      obj = read_v1(unit, beta, positive_only)
    ENDIF
  ELSE
    WRITE(*, *) "Invalid version number", version
    STOP "Stopping..."
  ENDIF
  !-----------------------------------------------------------------------
  END FUNCTION read_ir
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION read_v1(unit, beta, positive_only) RESULT(obj)
  !-----------------------------------------------------------------------
  !!
  !! This function reads the file to get the ir-basis objects.
  !! (version 1)
  !!
  !
  INTEGER, INTENT(IN) :: unit
  !! Unit number
  REAL(KIND = DP), INTENT(IN) :: beta
  !! inverse temperature
  LOGICAL, INTENT(IN), OPTIONAL :: positive_only
  !! if true, take the Matsubara frequencies
  !! only from the positive region
  !
  REAL(KIND = DP), PARAMETER :: rtol = 1e-20
  !! 
  TYPE(IR) :: obj
  !! to contain all the ir-basis objects from the file
  CHARACTER(LEN = 100) :: tmp_str
  !! dummy for characters
  INTEGER :: i
  !! counter
  INTEGER :: l
  !! counter
  INTEGER :: t
  !! counter
  INTEGER :: n
  !! counter
  INTEGER :: size
  !! total number of IR basis functions (size of s)
  INTEGER :: ntau
  !! total number of sampling points of imaginary time
  INTEGER :: nfreq_f
  !! total number of sampling Matsubara freqs (Fermionic)
  INTEGER :: nfreq_b
  !! total number of sampling Matsubara freqs (Bosonic)
  INTEGER :: nomega
  !! total number of sampling points of real frequency
  INTEGER, ALLOCATABLE :: freq_f(:)
  !! integer part of sampling Matsubara freqs (Fermion)
  INTEGER, ALLOCATABLE :: freq_b(:)
  !! integer part of sampling Matsubara freqs (Boson)
  REAL(KIND = DP) :: rtmp
  !! dummy for real variables
  REAL(KIND = DP) :: rtmp2
  !! dummy for real variables
  REAL(KIND = DP) :: lambda
  !! lambda = 10^{nlambda}, 
  !! which determines maximum sampling point of real frequency
  REAL(KIND = DP) :: eps
  !! eps = 10^{-ndigit}
  REAL(KIND = DP), ALLOCATABLE :: s(:)
  !! singular values
  REAL(KIND = DP), ALLOCATABLE :: tau(:)
  !! sampling points of imaginary time (dimenisonless)
  REAL(KIND = DP), ALLOCATABLE :: omega(:)
  !! sampling points of real frequency (dimenisonless)
  REAL(KIND = DP), ALLOCATABLE :: u(:, :)
  !! dimensionless ir-basis functions of tau
  REAL(KIND = DP), ALLOCATABLE :: v(:, :)
  !! this may be not used after getting dlr
  REAL(KIND = DP), ALLOCATABLE :: dlr(:, :)
  !! change-of-basis matrix from IR basis to DLR basis
  !! dlr(i, l) = - s(l) * v(i, l)
  COMPLEX(KIND = DP), ALLOCATABLE :: uhat_f(:, :)
  !! dimensionless ir-basis functions of Matsubara freqs
  COMPLEX(KIND = DP), ALLOCATABLE :: uhat_b(:, :)
  !! dimensionless ir-basis functions of Matsubara freqs
  !
  READ(unit,*) tmp_str, lambda
  READ(unit,*) tmp_str, eps
  !
  ! Singular values
  READ(unit,*)
  READ(unit,*) size
  ALLOCATE(s(size))
  DO i=1, size
    READ(unit, *) s(i)
  ENDDO
  !
  ! Sampling times
  READ(unit,*)
  READ(unit,*) ntau
  ALLOCATE(tau(ntau))
  DO i=1, ntau
    READ(unit, *) tau(i)
  ENDDO
  !
  ! Basis functions on sampling times
  READ(unit,*)
  ALLOCATE(u(ntau, size))
  DO l = 1, size
    DO t = 1, ntau
      READ(unit, *) rtmp
      u(t, l) = rtmp
    ENDDO
  ENDDO
  !
  ! Sampling frequencies (F)
  READ(unit,*)
  READ(unit,*) nfreq_f
  ALLOCATE(freq_f(nfreq_f))
  DO i=1, nfreq_f
    READ(unit, *) freq_f(i)
  ENDDO
  !
  READ(unit,*)
  ALLOCATE(uhat_f(nfreq_f, size))
  DO l = 1, size
    DO n = 1, nfreq_f
      READ(unit, *) rtmp, rtmp2
      uhat_f(n, l) = CMPLX(rtmp, rtmp2, KIND = DP)
    ENDDO
  ENDDO
  !
  ! Sampling frequencies (B)
  READ(unit,*)
  READ(unit,*) nfreq_b
  ALLOCATE(freq_b(nfreq_b))
  DO i=1, nfreq_b
    READ(unit, *) freq_b(i)
  ENDDO
  !
  READ(unit,*)
  ALLOCATE(uhat_b(nfreq_b, size))
  DO l = 1, size
    DO n = 1, nfreq_b
      READ(unit, *) rtmp, rtmp2
      uhat_b(n, l) = CMPLX(rtmp, rtmp2, KIND = DP)
    ENDDO
  ENDDO
  !
  ! Sampling poles on real frequencies
  READ(unit,*)
  READ(unit,*) nomega
  ALLOCATE(omega(nomega))
  DO i=1, nomega
    READ(unit, *) omega(i)
  ENDDO
  !
  ! Right singular functions on sampling poles
  READ(unit,*)
  ALLOCATE(v(nomega, size))
  ALLOCATE(dlr(nomega, size))
  DO l = 1, size
    DO i = 1, nomega
      READ(unit, *) rtmp
      v(i, l) = rtmp
      dlr(i, l) = - s(l) * v(i, l)
    ENDDO
  ENDDO
  !
  IF ((.NOT. PRESENT(positive_only))) then
    CALL init_ir(obj, beta, lambda, eps, s, tau, freq_f, freq_b, u, uhat_f, uhat_b, omega, v, dlr, 1d-16)
  ELSE
    CALL init_ir(obj, beta, lambda, eps, s, tau, freq_f, freq_b, u, uhat_f, uhat_b, omega, v, dlr, 1d-16, positive_only)
  ENDIF
  !
  DEALLOCATE(s, tau, u, freq_f, uhat_f, freq_b, uhat_b, omega, v, dlr)
  !
  !-----------------------------------------------------------------------
  END FUNCTION read_v1
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  END MODULE sparse_ir_io
  !-----------------------------------------------------------------------
