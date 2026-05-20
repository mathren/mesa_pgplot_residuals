! pgstar_residuals.inc
!
! HOW TO USE
! ----------
! 1. PASTE the following block into run_star_extras.f90 BEFORE its `contains`:
!
!   integer, parameter :: NR_MAX_MODELS = 200
!   integer, allocatable, save :: nr_zone_buf(:)
!   integer, allocatable, save :: nr_model_buf(:)     ! (NR_MAX_MODELS)
!   integer,              save :: nr_n_stored = 0, nr_n_cells = 0
!   real, allocatable, save :: nr_resid_buf(:,:)   ! (NR_MAX_MODELS, nz)
!   real, allocatable, save :: nr_mass_buf(:,:)
!   real,                 save :: nr_resid_min = -101.0, nr_resid_max = -101.0  ! auto-scale
!
! 2. AFTER `contains` in run_star_extras.f90, add:
!       include 'pgstar_residuals.inc'
!
! 3. In extras_controls, add:
!       s% other_pgstar_plots_info => nr_resid_pgstar_plots_info
!
! ---------------------------------------------------------------------------

! -----------  register plot ----------------------------
   subroutine nr_resid_pgstar_plots_info(id, ierr)
      use star_pgstar
      integer, intent(in)  :: id
      integer, intent(out) :: ierr

      type(pgstar_win_file_data), pointer :: p
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      p => s%pg%pgstar_win_file_ptr(i_Other)
      p%plot              => nr_resid_plot
      p%id                =  i_Other
      p%name              = 'Newton Raphson Residuals - accepted step'
      p%win_flag          = .true.
      p%win_width         =  12.0
      p%win_aspect_ratio  =  0.6
      p%file_flag         = .true.
      p%file_dir          = 'png'         ! folder where to save files
      p%file_prefix       = 'nr_resid_'
      p%file_interval     =  5
      p%file_width        = -1.0
      p%file_aspect_ratio = -1.0

   end subroutine nr_resid_pgstar_plots_info

! -----------  render plot ----------------------------
   subroutine nr_resid_plot(id, device_id, ierr)
      use const_def, only: dp
      integer, intent(in)  :: id, device_id
      integer, intent(out) :: ierr

      type(star_info), pointer :: s
      integer  :: k, neq, col

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      neq = size(s%equ, 1)   ! number of equations per cell

      ! --- allocate history buffer on first call ---
      if (.not. allocated(nr_resid_buf)) then
         allocate(nr_resid_buf(NR_MAX_MODELS, s%nz))
         allocate(nr_model_buf(NR_MAX_MODELS))
         allocate(nr_zone_buf(NR_MAX_MODELS))
         allocate(nr_mass_buf(NR_MAX_MODELS, s%nz))
         nr_resid_buf = -99.0
         nr_mass_buf  =  0.0   ! zero-init so unused cells are safely masked
         nr_n_cells   = s%nz
      end if

      ! --- roll buffer if full ---
      if (nr_n_stored < NR_MAX_MODELS) then
         nr_n_stored = nr_n_stored + 1
         col = nr_n_stored
      else
         nr_model_buf(1:NR_MAX_MODELS-1)      = nr_model_buf(2:NR_MAX_MODELS)
         nr_resid_buf(1:NR_MAX_MODELS-1, :)   = nr_resid_buf(2:NR_MAX_MODELS, :)
         nr_mass_buf(1:NR_MAX_MODELS-1, :)    = nr_mass_buf(2:NR_MAX_MODELS, :)
         col = NR_MAX_MODELS
      end if

      nr_model_buf(col) = s%model_number
      nr_zone_buf(col) = s%nz
      do k = 1, s%nz
         nr_resid_buf(col, k) = safe_log10(maxval(abs(s%equ(1:neq, k))))       ! residuals are cell-centered
         nr_mass_buf(col, k) = s%m(k) - 0.5 * s%dm(k) ! outer face mass coord minus half mass of the cell
      end do

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call nr_resid_render(ierr)
      call pgebuf()

   end subroutine nr_resid_plot

   subroutine nr_resid_render(ierr)
      use pgstar_colors
      use const_def, only: dp, Msun
      integer, intent(out) :: ierr

      real, allocatable :: img(:,:), mass_grid(:)
      real :: tr(6), fg, bg, xlo, xhi, dx
      real :: mlo, mhi, dm, m_j, frac
      integer :: nx, ny, i, j, k, nz_i

      ! ----------------------------------------- colorbar definition
      real :: bright, contra
      real, dimension(9) :: l, r, g, b

      l = (/ &
           0.00, 0.12, 0.25, 0.38, 0.50, &
           0.62, 0.75, 0.88, 1.00 /)

      r = (/ &
           0.00, 0.05, 0.15, 0.32, 0.55, &
           0.78, 0.93, 0.99, 1.00 /)

      g = (/ &
           0.00, 0.02, 0.04, 0.06, 0.10, &
           0.18, 0.32, 0.58, 0.95 /)

      b = (/ &
           0.00, 0.08, 0.18, 0.32, 0.36, &
           0.28, 0.16, 0.08, 0.80 /)

      contra = 1.0
      bright = 0.5

      call pgscir(16,255)
      call pgctab(l, r, g, b, 9, contra, bright)
      ! ----------------------------------------- end colorbar

      ierr = 0
      if (nr_n_stored < 2) return

      nx = nr_n_stored
      ny = nr_n_cells

      allocate(img(nx, ny))
      allocate(mass_grid(ny))

      mhi = -1.0e30
      mlo =  1.0e30
      do i = 1, nx
         nz_i = nr_zone_buf(i)
         do k = 1, nz_i
            if (nr_mass_buf(i,k) > mhi) mhi = nr_mass_buf(i,k)
            ! m(k) is positive; skip any zero-initialised slots
            if (nr_mass_buf(i,k) > 0.0 .and. nr_mass_buf(i,k) < mlo) &
               mlo = nr_mass_buf(i,k)
         end do
      end do
      mhi = mhi / real(Msun)   ! convert to solar masses for display
      mlo = mlo / real(Msun)

      dm = (mhi - mlo) / real(ny - 1)

      do j = 1, ny
         mass_grid(j) = mhi - real(j - 1) * dm
      end do

      ! --- interpolate residuals onto the uniform mass grid ---
      img = -99.0
      do i = 1, nx
         nz_i = nr_zone_buf(i)
         do j = 1, ny
            m_j = mass_grid(j) * real(Msun)   ! back to cgs for comparison
            do k = 1, nz_i - 1
               ! bracket: m(k) >= m_j >= m(k+1)
               if (nr_mass_buf(i,k) >= m_j .and. nr_mass_buf(i,k+1) <= m_j) then
                  if (nr_mass_buf(i,k) > nr_mass_buf(i,k+1)) then
                     frac = (m_j            - nr_mass_buf(i,k+1)) &
                          / (nr_mass_buf(i,k) - nr_mass_buf(i,k+1))
                  else
                     frac = 0.5
                  end if
                  img(i,j) = nr_resid_buf(i,k+1) &
                            + frac * (nr_resid_buf(i,k) - nr_resid_buf(i,k+1))
                  exit
               end if
            end do
         end do
      end do

      ! --- x range centred on model numbers ---
      dx  = max(1.0, real(nr_model_buf(nx) - nr_model_buf(1)) / real(nx - 1))
      xlo = real(nr_model_buf(1))  - 0.5*dx
      xhi = real(nr_model_buf(nx)) + 0.5*dx

      ! --- transformation matrix for pgimag ---
      ! world_y(j) = tr(4) + tr(6)*j
      ! j=1  -> mhi  =>  tr(4) = mhi - tr(6) = mhi + dm
      ! j=ny -> mlo  =>  confirms tr(6) = -dm
      tr(1) = xlo
      tr(2) = (xhi - xlo) / real(nx)
      tr(3) = 0.0
      tr(4) = mhi + dm
      tr(5) = 0.0
      tr(6) = -dm

      ! color scale
      fg = nr_resid_max
      bg = nr_resid_min
      if (fg <= -100) fg = maxval(img)
      if (bg <= -100) bg = minval(img, mask=(img > -19.9))
      if (bg >= fg)   bg = fg - 6.0

      ! --- main panel ---
      ! y-axis: mlo at bottom (centre), mhi at top (surface)
      call pgsvp(0.10, 0.82, 0.12, 0.92)
      call pgswin(xlo, xhi, mlo - 0.5*dm, mhi + 0.5*dm)
      call pgimag(img, nx, ny, 1, nx, 1, ny, fg, bg, tr)

      call pgsci(clr_Foreground)
      call pgsch(1.2)
      call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      call pglab('Model Number', 'Mass (M\d\(2281)\u)', &
                 'Newton Raphson Solver Residuals of accepted step')

      ! --- color wedge ---
      call pgsvp(0.84, 0.90, 0.12, 0.92)
      call pgswin(0.0, 1.0, bg, fg)
      call pgwedg('RI', 0.0, 4.0, bg, fg, 'log10(max(|res|))')

      deallocate(img)
      deallocate(mass_grid)

   end subroutine nr_resid_render
