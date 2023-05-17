  !-----------------------------------------------------------------------
  PROGRAM test
  !-----------------------------------------------------------------------
  USE sparse_ir
  USE sparse_ir_io
  USE sparse_ir_preset
  !
  IMPLICIT NONE
  !
  ! use preset
  CALL test_fermion(.true., .false., .false., .false.)
  CALL test_boson  (.true., .false., .false., .false.)
  ! use the file
  CALL test_fermion(.false., .false., .false., .false.)
  CALL test_boson  (.false., .false., .false., .false.)
  CALL test_fermion(.false., .true., .false., .false.)
  CALL test_boson  (.false., .true., .false., .false.)
  CALL test_fermion(.false., .true., .false., .true.)
  CALL test_boson  (.false., .true., .false., .true.)
  CALL test_fermion(.false., .true., .true., .false.)
  CALL test_boson  (.false., .true., .true., .false.)
  CALL test_fermion(.false., .true., .true., .true.)
  CALL test_boson  (.false., .true., .true., .true.)
  !
  ! use preset
  CALL test_fermion_dlr(.true., .false., .false., .false., .false.) ! positive_only = false
  CALL test_boson_dlr  (.true., .false., .false., .false., .false.) ! positive_only = false
  ! use the file
  CALL test_fermion_dlr(.false., .false., .false., .false., .false.) ! positive_only = false
  CALL test_boson_dlr  (.false., .false., .false., .false., .false.) ! positive_only = false
  CALL test_fermion_dlr(.false., .true., .false., .false., .false.) ! zz zz zz
  CALL test_boson_dlr  (.false., .true., .false., .false., .false.) ! zz zz zz
  !CALL test_fermion_dlr(.false., .true., .false., .true., .false.) ! zd dz zz
  !CALL test_boson_dlr  (.false., .true., .false., .true., .false.) ! zd dz zz
  !CALL test_fermion_dlr(.false., .true., .true., .false., .false.) ! dz zz zd
  !CALL test_boson_dlr  (.false., .true., .true., .false., .false.) ! dz zz zd
  CALL test_fermion_dlr(.false., .true., .true., .true., .false.) ! dd dz zd
  CALL test_boson_dlr  (.false., .true., .true., .true., .false.) ! dd dz zd
  !CALL test_fermion_dlr(.false., .true., .false., .false., .true.) ! zz zd dz
  !CALL test_boson_dlr  (.false., .true., .false., .false., .true.) ! zz zd dz
  CALL test_fermion_dlr(.false., .true., .false., .true., .true.) ! zd dd dz
  CALL test_boson_dlr  (.false., .true., .false., .true., .true.) ! zd dd dz
  CALL test_fermion_dlr(.false., .true., .true., .false., .true.) ! dz zd dd
  CALL test_boson_dlr  (.false., .true., .true., .false., .true.) ! dz zd dd
  !CALL test_fermion_dlr(.false., .true., .true., .true., .true.) ! dd dd dd
  !CALL test_boson_dlr  (.false., .true., .true., .true., .true.) ! dd dd dd
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE test_fermion(preset, positive_only, lflag_gtau, lflag_gl)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs some tests for fermionic function.
  !!
  LOGICAL, INTENT(IN) :: preset
  LOGICAL, INTENT(IN) :: positive_only
  LOGICAL, INTENT(IN) :: lflag_gtau
  LOGICAL, INTENT(IN) :: lflag_gl
  TYPE(IR) :: ir_obj
  INTEGER, PARAMETER :: ndigit = 10, nlambda = 4
  REAL(KIND = DP), PARAMETER :: lambda = 1.d1 ** nlambda
  REAL(KIND = DP), PARAMETER :: wmax = 1.d0
  !
  REAL(KIND = DP), PARAMETER :: beta = lambda/wmax, omega0 = 1.d0/beta
  REAL(KIND = DP), PARAMETER :: eps = 1.d-1**ndigit
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: giv(:,:), gl_matsu(:, :), gl_tau(:, :), gtau(:, :), &
    gtau_reconst(:, :), giv_reconst(:, :)
  REAL(KIND = DP), ALLOCATABLE :: gl_matsu_d(:, :), gl_tau_d(:, :), gtau_d(:, :), &
    gtau_reconst_d(:, :)
  INTEGER n, t, l
  !
  IF (preset) THEN
    ir_obj = mk_ir_preset(nlambda, ndigit, beta, positive_only)
  ELSE
    OPEN(99, file='ir_nlambda4_ndigit10.dat', status='old')
    ir_obj = read_ir(99, beta, positive_only)
    CLOSE(99)
  ENDIF
  !
  IF (abs(ir_obj%beta - beta) > 1d-10) THEN
    WRITE(*,*) "beta does not match"
    STOP 1
  ENDIF
  IF (abs(ir_obj%wmax - wmax) > 1d-10) THEN
    WRITE(*,*) "wmax does not match"
    STOP 1
  ENDIF
  !
  ! With ω0 = 1/β,
  !   G(iv) = 1/(iv - ω0),
  !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
  ALLOCATE(giv(1, ir_obj%nfreq_f))
  ALLOCATE(giv_reconst(1, ir_obj%nfreq_f))
  ALLOCATE(gtau(1, ir_obj%ntau))
  ALLOCATE(gl_matsu(1, ir_obj%size))
  ALLOCATE(gl_tau(1, ir_obj%size))
  ALLOCATE(gtau_reconst(1, ir_obj%ntau))
  ALLOCATE(gtau_d(1, ir_obj%ntau))
  ALLOCATE(gl_matsu_d(1, ir_obj%size))
  ALLOCATE(gl_tau_d(1, ir_obj%size))
  ALLOCATE(gtau_reconst_d(1, ir_obj%ntau))
  !
  ! From Matsubara
  DO n = 1, ir_obj%nfreq_f
    giv(1, n) = 1.d0/(CMPLX(0d0, pi*ir_obj%freq_f(n)/beta, KIND = DP) - omega0)
  ENDDO
  IF (lflag_gl) THEN
    CALL fit_matsubara_f(ir_obj, giv, gl_matsu_d)
  ELSE
    CALL fit_matsubara_f(ir_obj, giv, gl_matsu)
  ENDIF
  !
  ! From tau
  !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
  DO t = 1, ir_obj%ntau
    gtau(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 + EXP(-beta * omega0))
    gtau_d(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 + EXP(-beta * omega0))
  ENDDO
  IF (lflag_gtau .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau_d, gl_tau_d)
  ELSEIF ((.NOT. lflag_gtau) .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau, gl_tau_d)
  ELSEIF (lflag_gtau .AND. (.NOT. lflag_gl)) THEN
    CALL fit_tau(ir_obj, gtau_d, gl_tau)
  ELSE
    CALL fit_tau(ir_obj, gtau, gl_tau)
  ENDIF
  !
  IF (lflag_gl) THEN
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) gl_matsu_d(1,l)
    !  WRITE(*,*) gl_tau_d(1,l)
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu_d - gl_tau_d)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ELSE
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) REAL(gl_matsu(1,l)), AIMAG(gl_matsu(1,l))
    !  WRITE(*,*) REAL(gl_tau(1,l)), AIMAG(gl_tau(1,l))
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu - gl_tau)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ENDIF
  IF (lflag_gl) THEN
    CALL evaluate_matsubara_f(ir_obj, gl_matsu_d, giv_reconst)
  ELSE
    CALL evaluate_matsubara_f(ir_obj, gl_matsu, giv_reconst)
  ENDIF
  !DO n = 1, ir_obj%nfreq_f
  !  WRITE(*,*) giv(1,n)
  !  WRITE(*,*) giv_reconst(1,n)
  !ENDDO
  IF (MAXVAL(ABS(giv - giv_reconst)) > 1d2*eps) THEN
    WRITE(*,*) "giv do not match!"
    STOP 1
  ENDIF
  !
  IF (lflag_gtau .AND. lflag_gl) THEN
    CALL evaluate_tau(ir_obj, gl_tau_d, gtau_reconst_d)
  ELSEIF ((.NOT. lflag_gtau) .AND. lflag_gl) THEN
    CALL evaluate_tau(ir_obj, gl_tau_d, gtau_reconst)
    gtau_reconst_d = REAL(gtau_reconst, KIND = DP)
  ELSEIF (lflag_gtau .AND. (.NOT. lflag_gl)) THEN
    CALL evaluate_tau(ir_obj, gl_tau, gtau_reconst_d)
  ELSE
    CALL evaluate_tau(ir_obj, gl_tau, gtau_reconst)
    gtau_reconst_d = REAL(gtau_reconst, KIND = DP)
  ENDIF
  !DO n = 1, ir_obj%ntau
  !  WRITE(*,*) gtau_d(1,n)
  !  WRITE(*,*) gtau_reconst_d(1,n)
  !ENDDO
  IF (MAXVAL(ABS(gtau_d - gtau_reconst_d)) > 1d2*eps) THEN
    WRITE(*,*) "gtau do not match!"
    STOP 1
  ENDIF 
  !
  !WRITE(*,*) "test_fermion"
  !WRITE(*,*) "preset = ", preset
  !WRITE(*,*) "positive_only = ", positive_only
  !WRITE(*,*) "lflag_gtau = ", lflag_gtau
  !WRITE(*,*) "lflag_gl = ", lflag_gl
  !WRITE(*,*)
  !
  DEALLOCATE(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
  DEALLOCATE(gtau_d, gl_matsu_d, gl_tau_d, gtau_reconst_d)
  !
  CALL finalize_ir(ir_obj)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE test_fermion
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE test_boson(preset, positive_only, lflag_gtau, lflag_gl)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs some tests for bosonic function.
  !!
  LOGICAL, INTENT(IN) :: preset
  LOGICAL, INTENT(IN) :: positive_only
  LOGICAL, INTENT(IN) :: lflag_gtau
  LOGICAL, INTENT(IN) :: lflag_gl
  TYPE(IR) :: ir_obj
  INTEGER, PARAMETER :: ndigit = 10, nlambda = 4
  REAL(KIND = DP), PARAMETER :: lambda = 1.d1 ** nlambda
  REAL(KIND = DP), PARAMETER :: wmax = 1.d0
  !
  REAL(KIND = DP), PARAMETER :: beta = lambda/wmax, omega0 = 1.d0/beta
  REAL(KIND = DP), PARAMETER :: eps = 1.d-1**ndigit
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: giv(:,:), gl_matsu(:, :), gl_tau(:, :), gtau(:, :), &
    gtau_reconst(:, :), giv_reconst(:, :)
  REAL(KIND = DP), ALLOCATABLE :: gl_matsu_d(:, :), gl_tau_d(:, :), gtau_d(:, :), &
    gtau_reconst_d(:, :)
  INTEGER n, t, l
  !
  IF (preset) THEN
    ir_obj = mk_ir_preset(nlambda, ndigit, beta, positive_only)
  ELSE
    OPEN(99, file='ir_nlambda4_ndigit10.dat', status='old')
    ir_obj = read_ir(99, beta, positive_only)
    CLOSE(99)
  ENDIF
  !
  IF (abs(ir_obj%beta - beta) > 1d-10) THEN
    WRITE(*,*) "beta does not match"
    STOP 1
  ENDIF
  IF (abs(ir_obj%wmax - wmax) > 1d-10) THEN
    WRITE(*,*) "wmax does not match"
    STOP 1
  ENDIF
  !
  ! With ω0 = 1/β,
  !   G(iv) = 1/(iv - ω0),
  !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
  ALLOCATE(giv(1, ir_obj%nfreq_b))
  ALLOCATE(giv_reconst(1, ir_obj%nfreq_b))
  ALLOCATE(gtau(1, ir_obj%ntau))
  ALLOCATE(gl_matsu(1, ir_obj%size))
  ALLOCATE(gl_tau(1, ir_obj%size))
  ALLOCATE(gtau_reconst(1, ir_obj%ntau))
  ALLOCATE(gtau_d(1, ir_obj%ntau))
  ALLOCATE(gl_matsu_d(1, ir_obj%size))
  ALLOCATE(gl_tau_d(1, ir_obj%size))
  ALLOCATE(gtau_reconst_d(1, ir_obj%ntau))
  !
  ! From Matsubara
  DO n = 1, ir_obj%nfreq_b
    giv(1, n) = 1.d0/(CMPLX(0d0, pi*ir_obj%freq_b(n)/beta, KIND = DP) - omega0)
  ENDDO
  IF (lflag_gl) THEN
    CALL fit_matsubara_b(ir_obj, giv, gl_matsu_d)
  ELSE
    CALL fit_matsubara_b(ir_obj, giv, gl_matsu)
  ENDIF
  !
  ! From tau
  !   G(τ=0) = - exp(-τ ω0)/(1-exp(-β ω0)),
  DO t = 1, ir_obj%ntau
    gtau(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 - EXP(-beta * omega0))
    gtau_d(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 - EXP(-beta * omega0))
  ENDDO
  IF (lflag_gtau .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau_d, gl_tau_d)
  ELSEIF ((.NOT. lflag_gtau) .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau, gl_tau_d)
  ELSEIF (lflag_gtau .AND. (.NOT. lflag_gl)) THEN
    CALL fit_tau(ir_obj, gtau_d, gl_tau)
  ELSE
    CALL fit_tau(ir_obj, gtau, gl_tau)
  ENDIF 
  !
  IF (lflag_gl) THEN
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) gl_matsu_d(1,l)
    !  WRITE(*,*) gl_tau_d(1,l)
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu_d - gl_tau_d)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ELSE
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) REAL(gl_matsu(1,l)), AIMAG(gl_matsu(1,l))
    !  WRITE(*,*) REAL(gl_tau(1,l)), AIMAG(gl_tau(1,l))
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu - gl_tau)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ENDIF
  IF (lflag_gl) THEN
    CALL evaluate_matsubara_b(ir_obj, gl_matsu_d, giv_reconst)
  ELSE
    CALL evaluate_matsubara_b(ir_obj, gl_matsu, giv_reconst)
  ENDIF
  !DO n = 1, ir_obj%nfreq_f
  !  WRITE(*,*) giv(1,n)
  !  WRITE(*,*) giv_reconst(1,n)
  !ENDDO
  IF (MAXVAL(ABS(giv - giv_reconst)) > 1d2*eps) THEN
    WRITE(*,*) "giv do not match!"
    STOP 1
  ENDIF
  !
  IF (lflag_gtau .AND. lflag_gl) THEN
    CALL evaluate_tau(ir_obj, gl_tau_d, gtau_reconst_d)
  ELSEIF ((.NOT. lflag_gtau) .AND. lflag_gl) THEN
    CALL evaluate_tau(ir_obj, gl_tau_d, gtau_reconst)
    gtau_reconst_d = REAL(gtau_reconst, KIND = DP)
  ELSEIF (lflag_gtau .AND. (.NOT. lflag_gl)) THEN
    CALL evaluate_tau(ir_obj, gl_tau, gtau_reconst_d)
  ELSE
    CALL evaluate_tau(ir_obj, gl_tau, gtau_reconst)
    gtau_reconst_d = REAL(gtau_reconst, KIND = DP)
  ENDIF
  !DO n = 1, ir_obj%ntau
  !  WRITE(*,*) gtau_d(1,n)
  !  WRITE(*,*) gtau_reconst_d(1,n)
  !ENDDO
  IF (MAXVAL(ABS(gtau_d - gtau_reconst_d)) > 1d2*eps) THEN
    WRITE(*,*) "gtau do not match!"
    STOP 1
  ENDIF
  !
  !WRITE(*,*) "test_boson"
  !WRITE(*,*) "preset = ", preset
  !WRITE(*,*) "positive_only = ", positive_only
  !WRITE(*,*) "lflag_gtau = ", lflag_gtau
  !WRITE(*,*) "lflag_gl = ", lflag_gl
  !WRITE(*,*)
  !
  DEALLOCATE(giv, gtau, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
  DEALLOCATE(gtau_d, gl_matsu_d, gl_tau_d, gtau_reconst_d)
  !
  CALL finalize_ir(ir_obj)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE test_boson
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE test_fermion_dlr(preset, positive_only, lflag_gtau, lflag_gl, lflag_gdlr)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs some tests of DLR for fermionic function.
  !!
  LOGICAL, INTENT(IN) :: preset
  LOGICAL, INTENT(IN) :: positive_only
  LOGICAL, INTENT(IN) :: lflag_gtau
  LOGICAL, INTENT(IN) :: lflag_gl
  LOGICAL, INTENT(IN) :: lflag_gdlr
  TYPE(IR) :: ir_obj
  INTEGER, PARAMETER :: ndigit = 10, nlambda = 4
  REAL(KIND = DP), PARAMETER :: lambda = 1.d1 ** nlambda
  REAL(KIND = DP), PARAMETER :: wmax = 1.d0
  !
  REAL(KIND = DP), PARAMETER :: beta = lambda/wmax, omega0 = 1.d0/beta
  REAL(KIND = DP), PARAMETER :: eps = 1.d-1**ndigit
  INTEGER, PARAMETER :: ntau_dlr = 200, nfreq_dlr = 200
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: giv_smpl(:,:), gl_matsu(:, :), gl_tau(:, :), gtau_smpl(:, :), &
    gtau_reconst(:, :), giv_reconst(:, :), g_dlr(:, :), giv_ref(:,:), gtau_ref(:,:)
  REAL(KIND = DP), ALLOCATABLE :: gl_matsu_d(:, :), gl_tau_d(:, :), gtau_smpl_d(:, :), &
    gtau_reconst_d(:, :), g_dlr_d(:, :)
  INTEGER, ALLOCATABLE :: freq(:) 
  REAL(KIND = DP), ALLOCATABLE :: tau(:)
  INTEGER n, t
  !
  IF (preset) THEN
    ir_obj = mk_ir_preset(nlambda, ndigit, beta, positive_only)
  ELSE
    OPEN(99, file='ir_nlambda4_ndigit10.dat', status='old')
    ir_obj = read_ir(99, beta, positive_only)
    CLOSE(99)
  ENDIF
  !
  IF (abs(ir_obj%beta - beta) > 1d-10) THEN
    WRITE(*,*) "beta does not match"
    STOP 1
  ENDIF
  IF (abs(ir_obj%wmax - wmax) > 1d-10) THEN
    WRITE(*,*) "wmax does not match"
    STOP 1
  ENDIF
  !
  ! With ω0 = 1/β,
  !   G(iv) = 1/(iv - ω0),
  !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
  ALLOCATE(giv_smpl(1, ir_obj%nfreq_f))
  ALLOCATE(gtau_smpl(1, ir_obj%ntau))
  ALLOCATE(gl_matsu(1, ir_obj%size))
  ALLOCATE(gl_tau(1, ir_obj%size))
  ALLOCATE(g_dlr(1, ir_obj%nomega))
  ALLOCATE(giv_ref(1, nfreq_dlr))
  ALLOCATE(gtau_ref(1, ntau_dlr))
  ALLOCATE(gtau_reconst(1, ntau_dlr))
  ALLOCATE(giv_reconst(1, nfreq_dlr))
  ALLOCATE(freq(nfreq_dlr))
  ALLOCATE(tau(ntau_dlr))
  ALLOCATE(gtau_smpl_d(1, ir_obj%ntau))
  ALLOCATE(gl_matsu_d(1, ir_obj%size))
  ALLOCATE(gl_tau_d(1, ir_obj%size))
  ALLOCATE(g_dlr_d(1, ir_obj%nomega))
  ALLOCATE(gtau_reconst_d(1, ntau_dlr))
  !
  ! From Matsubara
  DO n = 1, ir_obj%nfreq_f
    giv_smpl(1, n) = 1.d0/(CMPLX(0d0, pi*ir_obj%freq_f(n)/beta, KIND = DP) - omega0)
  ENDDO
  IF (lflag_gl) THEN
    CALL fit_matsubara_f(ir_obj, giv_smpl, gl_matsu_d)
  ELSE
    CALL fit_matsubara_f(ir_obj, giv_smpl, gl_matsu)
  ENDIF
  !
  ! From tau
  !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
  DO t = 1, ir_obj%ntau
    gtau_smpl(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 + EXP(-beta * omega0))
    gtau_smpl_d(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 + EXP(-beta * omega0))
  ENDDO
  IF (lflag_gtau .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau_smpl_d, gl_tau_d)
  ELSEIF ((.NOT. lflag_gtau) .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau_smpl, gl_tau_d)
  ELSEIF (lflag_gtau .AND. (.NOT. lflag_gl)) THEN
    CALL fit_tau(ir_obj, gtau_smpl_d, gl_tau)
  ELSE
    CALL fit_tau(ir_obj, gtau_smpl, gl_tau)
  ENDIF
  !
  IF (lflag_gl) THEN
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) gl_matsu_d(1,l)
    !  WRITE(*,*) gl_tau_d(1,l)
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu_d - gl_tau_d)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ELSE
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) REAL(gl_matsu(1,l)), AIMAG(gl_matsu(1,l))
    !  WRITE(*,*) REAL(gl_tau(1,l)), AIMAG(gl_tau(1,l))
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu - gl_tau)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ENDIF
  !
  IF (lflag_gl .AND. lflag_gdlr) THEN
    CALL to_dlr(ir_obj, gl_matsu_d, g_dlr_d)
  ELSEIF ((.NOT. lflag_gl) .AND. lflag_gdlr) THEN
    CALL to_dlr(ir_obj, gl_matsu, g_dlr_d)
  ELSEIF (lflag_gl .AND. (.NOT. lflag_gdlr)) THEN
    CALL to_dlr(ir_obj, gl_matsu_d, g_dlr)
  ELSE
    CALL to_dlr(ir_obj, gl_matsu, g_dlr)
  ENDIF
  !
  DO n = 1, nfreq_dlr
    freq(n) = -nfreq_dlr + 2 * n - 1
    giv_ref(1, n) = 1.d0/(CMPLX(0d0, pi*freq(n)/beta, KIND = DP) - omega0)
  ENDDO
  !
  IF (lflag_gdlr) THEN
    CALL evaluate_matsubara_f_from_dlr(ir_obj, freq, g_dlr_d, giv_reconst)
  ELSE
    CALL evaluate_matsubara_f_from_dlr(ir_obj, freq, g_dlr, giv_reconst)
  ENDIF
  IF (MAXVAL(ABS(giv_ref - giv_reconst)) > 1d3*eps) THEN
    WRITE(*,*) "giv do not match!"
    STOP 1
  ENDIF
  !
  DO n = 1, ntau_dlr
    tau(n) = beta * DBLE(n) / DBLE(ntau_dlr + 1)
    gtau_ref(1, n) = - EXP(-tau(n) * omega0)/(1.d0 + EXP(-beta * omega0))
  ENDDO
  !
  IF (lflag_gdlr .AND. lflag_gtau) THEN
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr_d, gtau_reconst_d)
    gtau_reconst = CMPLX(gtau_reconst_d, 0.0d0, KIND = DP)
  ELSEIF ((.NOT. lflag_gdlr) .AND. lflag_gtau) THEN
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr, gtau_reconst_d)
    gtau_reconst = CMPLX(gtau_reconst_d, 0.0d0, KIND = DP)
  ELSEIF (lflag_gdlr .AND. (.NOT. lflag_gtau)) THEN
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr_d, gtau_reconst)
  ELSE
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr, gtau_reconst)
  ENDIF
  IF (MAXVAL(ABS(gtau_ref - gtau_reconst)) > 1d3*eps) THEN
    WRITE(*,*) "gtau do not match!"
    STOP 1
  ENDIF
  !
  !WRITE(*,*) "test_fermion_dlr"
  !WRITE(*,*) "preset = ", preset
  !WRITE(*,*) "positive_only = ", positive_only
  !WRITE(*,*) "lflag_gtau = ", lflag_gtau
  !WRITE(*,*) "lflag_gl = ", lflag_gl
  !WRITE(*,*) "lflag_gdlr = ", lflag_gdlr
  !WRITE(*,*)
  !
  DEALLOCATE(giv_smpl, gtau_smpl, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
  DEALLOCATE(giv_ref, gtau_ref, g_dlr, freq, tau)
  DEALLOCATE(gtau_smpl_d, gl_matsu_d, gl_tau_d, gtau_reconst_d, g_dlr_d)
  !
  CALL finalize_ir(ir_obj)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE test_fermion_dlr
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE test_boson_dlr(preset, positive_only, lflag_gtau, lflag_gl, lflag_gdlr)
  !-----------------------------------------------------------------------
  !!
  !! This routine performs some tests of DLR for bosonic function.
  !!
  LOGICAL, INTENT(IN) :: preset
  LOGICAL, INTENT(IN) :: positive_only
  LOGICAL, INTENT(IN) :: lflag_gtau
  LOGICAL, INTENT(IN) :: lflag_gl
  LOGICAL, INTENT(IN) :: lflag_gdlr
  TYPE(IR) :: ir_obj
  INTEGER, PARAMETER :: ndigit = 10, nlambda = 4
  REAL(KIND = DP), PARAMETER :: lambda = 1.d1 ** nlambda
  REAL(KIND = DP), PARAMETER :: wmax = 1.d0
  !
  REAL(KIND = DP), PARAMETER :: beta = lambda/wmax, omega0 = 1.d0/beta
  REAL(KIND = DP), PARAMETER :: eps = 1.d-1**ndigit
  INTEGER, PARAMETER :: ntau_dlr = 200, nfreq_dlr = 200
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: giv_smpl(:,:), gl_matsu(:, :), gl_tau(:, :), gtau_smpl(:, :), &
    gtau_reconst(:, :), giv_reconst(:, :), g_dlr(:, :), giv_ref(:,:), gtau_ref(:,:)
  REAL(KIND = DP), ALLOCATABLE :: gl_matsu_d(:, :), gl_tau_d(:, :), gtau_smpl_d(:, :), &
    gtau_reconst_d(:, :), g_dlr_d(:, :)
  INTEGER, ALLOCATABLE :: freq(:) 
  REAL(KIND = DP), ALLOCATABLE :: tau(:)
  INTEGER n, t
  !
  IF (preset) THEN
    ir_obj = mk_ir_preset(nlambda, ndigit, beta, positive_only)
  ELSE
    OPEN(99, file='ir_nlambda4_ndigit10.dat', status='old')
    ir_obj = read_ir(99, beta, positive_only)
    CLOSE(99)
  ENDIF
  !
  IF (abs(ir_obj%beta - beta) > 1d-10) THEN
    WRITE(*,*) "beta does not match"
    STOP 1
  ENDIF
  IF (abs(ir_obj%wmax - wmax) > 1d-10) THEN
    WRITE(*,*) "wmax does not match"
    STOP 1
  ENDIF
  !
  ! With ω0 = 1/β,
  !   G(iv) = 1/(iv - ω0),
  !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
  ALLOCATE(giv_smpl(1, ir_obj%nfreq_b))
  ALLOCATE(gtau_smpl(1, ir_obj%ntau))
  ALLOCATE(gl_matsu(1, ir_obj%size))
  ALLOCATE(gl_tau(1, ir_obj%size))
  ALLOCATE(g_dlr(1, ir_obj%nomega))
  ALLOCATE(giv_ref(1, nfreq_dlr))
  ALLOCATE(gtau_ref(1, ntau_dlr))
  ALLOCATE(gtau_reconst(1, ntau_dlr))
  ALLOCATE(giv_reconst(1, nfreq_dlr))
  ALLOCATE(freq(nfreq_dlr))
  ALLOCATE(tau(ntau_dlr))
  ALLOCATE(gtau_smpl_d(1, ir_obj%ntau))
  ALLOCATE(gl_matsu_d(1, ir_obj%size))
  ALLOCATE(gl_tau_d(1, ir_obj%size))
  ALLOCATE(g_dlr_d(1, ir_obj%nomega))
  ALLOCATE(gtau_reconst_d(1, ntau_dlr))
  !
  ! From Matsubara
  DO n = 1, ir_obj%nfreq_b
    giv_smpl(1, n) = 1.d0/(CMPLX(0d0, pi*ir_obj%freq_b(n)/beta, KIND = DP) - omega0)
  ENDDO
  IF (lflag_gl) THEN
    CALL fit_matsubara_b(ir_obj, giv_smpl, gl_matsu_d)
  ELSE
    CALL fit_matsubara_b(ir_obj, giv_smpl, gl_matsu)
  ENDIF
  !
  ! From tau
  !   G(τ=0) = - exp(-τ ω0)/(1+exp(-β ω0)),
  DO t = 1, ir_obj%ntau
    gtau_smpl(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 - EXP(-beta * omega0))
    gtau_smpl_d(1, t) = - EXP(-ir_obj%tau(t) * omega0)/(1.d0 - EXP(-beta * omega0))
  ENDDO
  IF (lflag_gtau .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau_smpl_d, gl_tau_d)
  ELSEIF ((.NOT. lflag_gtau) .AND. lflag_gl) THEN
    CALL fit_tau(ir_obj, gtau_smpl, gl_tau_d)
  ELSEIF (lflag_gtau .AND. (.NOT. lflag_gl)) THEN
    CALL fit_tau(ir_obj, gtau_smpl_d, gl_tau)
  ELSE
    CALL fit_tau(ir_obj, gtau_smpl, gl_tau)
  ENDIF
  !
  IF (lflag_gl) THEN
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) gl_matsu_d(1,l)
    !  WRITE(*,*) gl_tau_d(1,l)
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu_d - gl_tau_d)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ELSE
    !DO l = 1, ir_obj%size
    !  WRITE(*,*) REAL(gl_matsu(1,l)), AIMAG(gl_matsu(1,l))
    !  WRITE(*,*) REAL(gl_tau(1,l)), AIMAG(gl_tau(1,l))
    !ENDDO
    IF (MAXVAL(ABS(gl_matsu - gl_tau)) > 1d2*eps) THEN
      WRITE(*,*) "gl_matsu and gl_tau do not match!"
      STOP 1
    ENDIF
  ENDIF
  !
  IF (lflag_gl .AND. lflag_gdlr) THEN
    CALL to_dlr(ir_obj, gl_matsu_d, g_dlr_d)
  ELSEIF ((.NOT. lflag_gl) .AND. lflag_gdlr) THEN
    CALL to_dlr(ir_obj, gl_matsu, g_dlr_d)
  ELSEIF (lflag_gl .AND. (.NOT. lflag_gdlr)) THEN
    CALL to_dlr(ir_obj, gl_matsu_d, g_dlr)
  ELSE
    CALL to_dlr(ir_obj, gl_matsu, g_dlr)
  ENDIF
  !
  DO n = 1, nfreq_dlr
    freq(n) = -nfreq_dlr + 2 * n
    giv_ref(1, n) = 1.d0/(CMPLX(0d0, pi*freq(n)/beta, KIND = DP) - omega0)
  ENDDO
  !
  IF (lflag_gdlr) THEN
    CALL evaluate_matsubara_b_from_dlr(ir_obj, freq, g_dlr_d, giv_reconst)
  ELSE
    CALL evaluate_matsubara_b_from_dlr(ir_obj, freq, g_dlr, giv_reconst)
  ENDIF
  IF (MAXVAL(ABS(giv_ref - giv_reconst)) > 1d3*eps) THEN
    WRITE(*,*) "giv do not match!"
    STOP 1
  ENDIF
  !
  DO n = 1, ntau_dlr
    tau(n) = beta * DBLE(n) / DBLE(ntau_dlr + 1)
    gtau_ref(1, n) = - EXP(-tau(n) * omega0)/(1.d0 - EXP(-beta * omega0))
  ENDDO
  !
  IF (lflag_gdlr .AND. lflag_gtau) THEN
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr_d, gtau_reconst_d)
    gtau_reconst = CMPLX(gtau_reconst_d, 0.0d0, KIND = DP)
  ELSEIF ((.NOT. lflag_gdlr) .AND. lflag_gtau) THEN
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr, gtau_reconst_d)
    gtau_reconst = CMPLX(gtau_reconst_d, 0.0d0, KIND = DP)
  ELSEIF (lflag_gdlr .AND. (.NOT. lflag_gtau)) THEN
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr_d, gtau_reconst)
  ELSE
    CALL evaluate_tau_from_dlr(ir_obj, tau, g_dlr, gtau_reconst)
  ENDIF
  IF (MAXVAL(ABS(gtau_ref - gtau_reconst)) > 1d3*eps) THEN
    WRITE(*,*) "gtau do not match!"
    STOP 1
  ENDIF
  !
  !WRITE(*,*) "test_boson_dlr"
  !WRITE(*,*) "preset = ", preset
  !WRITE(*,*) "positive_only = ", positive_only
  !WRITE(*,*) "lflag_gtau = ", lflag_gtau
  !WRITE(*,*) "lflag_gl = ", lflag_gl
  !WRITE(*,*) "lflag_gdlr = ", lflag_gdlr
  !WRITE(*,*)
  !
  DEALLOCATE(giv_smpl, gtau_smpl, gl_matsu, gl_tau, gtau_reconst, giv_reconst)
  DEALLOCATE(giv_ref, gtau_ref, g_dlr, freq, tau)
  DEALLOCATE(gtau_smpl_d, gl_matsu_d, gl_tau_d, gtau_reconst_d, g_dlr_d)
  !
  CALL finalize_ir(ir_obj)
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE test_boson_dlr
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  END PROGRAM test
  !-----------------------------------------------------------------------