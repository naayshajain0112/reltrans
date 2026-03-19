module lag_freq_mod
    use constants
    use lib_types, only: config_type, model_args_type, arrays_type
    implicit none

    ! gslope and ABslope are physical constants, not free parameters
    real, parameter :: gslope  = 1.0
    real, parameter :: ABslope = 1.0

contains

! ----------------------------------------------------------------------
! Main entry point: coherent lag-frequency calculation
! ----------------------------------------------------------------------
subroutine lag_freq(config, model_args, arrays, absorbx)
    type(config_type),     intent(in)    :: config
    type(model_args_type), intent(inout) :: model_args
    type(arrays_type),     intent(inout) :: arrays
    real,                  intent(in)    :: absorbx(:)

    integer :: Ea1, Ea2, Eb1, Eb2
    integer :: j
    real    :: f
    real    :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb
    real    :: f_ref   ! reference frequency for spectral scaling

    call energy_bounds(config, Ea1, Ea2, Eb1, Eb2)

    ! Compute reference frequency once (midpoint of first bin)
    f_ref = 0.5 * (arrays%fix(1) + arrays%fix(2))   ! fix is 1-based

    ! Apply eta boost to all reflection (m=2..nlp) components
    if (model_args%nlp .gt. 1) then
        call apply_eta_boost(model_args, arrays)
    end if

    do j = 1, config%nf
        f = config%flo * (config%fhi / config%flo)**((real(j) - 0.5) / real(config%nf))

        call accumulate_band(config, model_args, arrays, absorbx, &
                             Ea1, Ea2, j, f, f_ref,               &
                             ReGrawEa, ImGrawEa)

        call accumulate_band(config, model_args, arrays, absorbx, &
                             Eb1, Eb2, j, f, f_ref,               &
                             ReGrawEb, ImGrawEb)

        ! Cross-spectrum: G_AB = G_A* x G_B  (real and imaginary parts)
        arrays%ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
        arrays%ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
    end do
end subroutine lag_freq


! ----------------------------------------------------------------------
! Apply eta boost to the reflection (m >= 2) W arrays in-place
! ----------------------------------------------------------------------
subroutine apply_eta_boost(model_args, arrays)
    type(model_args_type), intent(in)    :: model_args
    type(arrays_type),     intent(inout) :: arrays

    arrays%ReW0(2,:,:) = model_args%eta * arrays%ReW0(2,:,:)
    arrays%ImW0(2,:,:) = model_args%eta * arrays%ImW0(2,:,:)
    arrays%ReW1(2,:,:) = model_args%eta * arrays%ReW1(2,:,:)
    arrays%ImW1(2,:,:) = model_args%eta * arrays%ImW1(2,:,:)
    arrays%ReW2(2,:,:) = model_args%eta * arrays%ReW2(2,:,:)
    arrays%ImW2(2,:,:) = model_args%eta * arrays%ImW2(2,:,:)
    arrays%ReW3(2,:,:) = model_args%eta * arrays%ReW3(2,:,:)
    arrays%ImW3(2,:,:) = model_args%eta * arrays%ImW3(2,:,:)
end subroutine apply_eta_boost


! ----------------------------------------------------------------------
! Accumulate the raw Green's function over one energy band
!
! Arguments:
!   i1, i2   : energy bin range (inclusive)
!   j        : frequency bin index (selects W array slice)
!   f        : current Fourier frequency (Hz)
!   f_ref    : reference frequency for gslope/ABslope scaling
!   ReGraw   : output – real part of accumulated Green's function
!   ImGraw   : output – imaginary part
! ----------------------------------------------------------------------
subroutine accumulate_band(config, model_args, arrays, absorbx, &
                            i1, i2, j, f, f_ref,                &
                            ReGraw, ImGraw)
    type(config_type),     intent(in) :: config
    type(model_args_type), intent(in) :: model_args
    type(arrays_type),     intent(in) :: arrays
    real,                  intent(in) :: absorbx(:)
    integer,               intent(in) :: i1, i2, j
    real,                  intent(in) :: f, f_ref
    real,                  intent(out):: ReGraw, ImGraw

    integer :: i, m
    real    :: E, fac, DelAB_nu, g_nu
    real    :: tau_d, tau_p, phase_d, phase_p
    complex :: W0, W1, W2, W3, Sraw, Stemp
    complex :: cexp_d, cexp_p, cexp_phi

    real    :: f_ratio   ! (f_ref / f) used for spectral scaling

    f_ratio = f_ref / f
    ReGraw  = 0.0
    ImGraw  = 0.0

    do i = i1, i2
        Sraw = (0.0, 0.0)

        ! Mid-bin energy: use bin i and i-1, with a safe lower guard
        if (i .gt. 1) then
            E = 0.5 * (arrays%earx(i) + arrays%earx(i-1))
        else
            ! Edge bin: use half the width from the bin above as an approximation
            E = arrays%earx(i) - 0.5 * (arrays%earx(i+1) - arrays%earx(i))
            E = max(E, arrays%earx(i) * 0.5)   ! floor at half the bin value
        end if

        do m = 1, model_args%nlp
            ! Frequency-dependent phase and amplitude scaling
            DelAB_nu = model_args%DelAB(m) * (f_ratio**ABslope)
            g_nu     = model_args%g(m)     * (f_ratio**gslope)

            ! log-ionisation factor
            fac = log(model_args%gso(m) / ((1.0 + model_args%z) * E))

            ! Phase delays are zero for the primary (m=1) component
            phase_d = 0.0
            phase_p = 0.0
            if (m .gt. 1) then
                tau_d   = model_args%tauso(m) - model_args%tauso(1)
                tau_p   = (model_args%h(m)    - model_args%h(1)) / model_args%beta_p
                phase_d = 2.0 * pi * tau_d * f
                phase_p = 2.0 * pi * tau_p * f
            end if

            cexp_d   = cmplx(cos(phase_d),   sin(phase_d))
            cexp_p   = cmplx(cos(phase_p),   sin(phase_p))
            cexp_phi = cmplx(cos(DelAB_nu),  sin(DelAB_nu))

            W0 = model_args%Anorm * cmplx(arrays%ReW0(m,i,j), arrays%ImW0(m,i,j))
            W1 = model_args%Anorm * cmplx(arrays%ReW1(m,i,j), arrays%ImW1(m,i,j))
            W2 = model_args%Anorm * cmplx(arrays%ReW2(m,i,j), arrays%ImW2(m,i,j))
            W3 = model_args%ionvar * model_args%Anorm &
                                   * cmplx(arrays%ReW3(m,i,j), arrays%ImW3(m,i,j))

            ! Reflection + continuum term (phase-delayed and AB-shifted)
            Stemp = g_nu * cexp_phi * (W1 + W2 + fac * cexp_d * arrays%contx(i,m))
            ! Add direct + ionisation + unshifted continuum
            Stemp = Stemp + W0 + W3 + cexp_d * arrays%contx(i,m)
            ! Sum over reflection components with propagation phase
            Sraw  = Sraw + cexp_p * Stemp
        end do

        ReGraw = ReGraw + real(Sraw)  * absorbx(i)
        ImGraw = ImGraw + aimag(Sraw) * absorbx(i)
    end do

end subroutine accumulate_band


! ----------------------------------------------------------------------
! Return energy-band index bounds from config
! ----------------------------------------------------------------------
subroutine energy_bounds(config, Ea1, Ea2, Eb1, Eb2)
    type(config_type), intent(in)  :: config
    integer,           intent(out) :: Ea1, Ea2, Eb1, Eb2

    Ea1 = config%ie1
    Ea2 = config%ie2
    Eb1 = config%ie3
    Eb2 = config%ie4
end subroutine energy_bounds


! ----------------------------------------------------------------------
! Incoherent (no-cross-term) lag-frequency calculation
! Mirrors lag_freq but accumulates power spectra separately
! ----------------------------------------------------------------------
subroutine lag_freq_nocoh(config, model_args, arrays, absorbx)
    type(config_type),     intent(in)    :: config
    type(model_args_type), intent(inout) :: model_args
    type(arrays_type),     intent(inout) :: arrays
    real,                  intent(in)    :: absorbx(:)

    integer :: Ea1, Ea2, Eb1, Eb2, j
    real    :: f, f_ref
    real    :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb

    call energy_bounds(config, Ea1, Ea2, Eb1, Eb2)

    f_ref = 0.5 * (arrays%fix(1) + arrays%fix(2))

    if (model_args%nlp .gt. 1) then
        call apply_eta_boost(model_args, arrays)
    end if

    do j = 1, config%nf
        f = config%flo * (config%fhi / config%flo)**((real(j) - 0.5) / real(config%nf))

        call accumulate_band(config, model_args, arrays, absorbx, &
                             Ea1, Ea2, j, f, f_ref,               &
                             ReGrawEa, ImGrawEa)

        call accumulate_band(config, model_args, arrays, absorbx, &
                             Eb1, Eb2, j, f, f_ref,               &
                             ReGrawEb, ImGrawEb)

        ! Incoherent: store power spectra, not the cross-spectrum
        arrays%ReGraw(j) = ReGrawEa**2 + ImGrawEa**2   ! |G_A|^2
        arrays%ImGraw(j) = ReGrawEb**2 + ImGrawEb**2   ! |G_B|^2
    end do
end subroutine lag_freq_nocoh

end module lag_freq_mod