! pgstar_residuals.inc
!
! HOW TO USE
! ----------
! 1. PASTE the following block into run_star_extras.f90 BEFORE its `contains`:
!
! ! for Kippenhahn style local value of  max across equations of residual
! integer, parameter :: NR_MAX_MODELS = 2000
! integer, allocatable, save :: nr_zone_buf(:)
! integer, allocatable, save :: nr_model_buf(:)
! integer,              save :: nr_n_stored = 0, nr_n_cells = 0
! real, allocatable, save :: nr_resid_buf(:,:), tmp_resid_buf(:,:)
! real, allocatable, save :: nr_mass_buf(:,:), tmp_mass_buf(:,:)
! real,              save :: nr_resid_min = -101.0, nr_resid_max = -101.0  ! auto-scale

! ! for max residual across all cells for each equation
! integer, parameter :: max_resid_hist = 2000
! integer            :: n_resid_hist   = 0
! integer            :: resid_hist_nvar = 0
! integer,       allocatable :: resid_hist_model(:)        ! (max_resid_hist)
! real,          allocatable :: resid_hist_vals(:,:)       ! (nvar, max_resid_hist)
! character(len=32), allocatable :: resid_equ_names(:)
!
! this creates some global variables accessible for the routines below
!
! 2. AFTER `contains` in run_star_extras.f90, add:
!       include 'pgstar_residuals.inc'
!
! 3. In extras_controls:
!    a. to get a "Kippenhahn-style" map of mass coordinate vs. model number
!       colored by the maximum residual (across all equations) for each location:
!
!       s% other_pgstar_plots_info => nr_resid_pgstar_plots_info
!
!    b. to get a plot as a function of model number of the maximum
!       (across mesh points) residual for each structure equation
!
!        s% other_pgstar_plots_info => equ_resid_hist
!
! ---------------------------------------------------------------------------

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

! Kippenhahn-style log10(max(|res|)) across eqs plot
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
     allocate(nr_resid_buf(NR_MAX_MODELS, 5 * s%nz)) ! make 5 times bigger than current nz
     allocate(nr_mass_buf(NR_MAX_MODELS, 5 * s%nz))
     allocate(nr_model_buf(NR_MAX_MODELS))
     allocate(nr_zone_buf(NR_MAX_MODELS))
     nr_resid_buf = -99.0
     nr_mass_buf  =  0.0   ! zero-init so unused cells are safely masked
     nr_n_cells   = s%nz
  else if (s%nz > 0.95 * size(nr_resid_buf, 2)) then
     ! check if remeshing has increased nz beyond 95% of the
     ! current array size, reallocate
     allocate(tmp_resid_buf(NR_MAX_MODELS, 5 * s%nz))
     allocate(tmp_mass_buf(NR_MAX_MODELS, 5  *s%nz))

     tmp_resid_buf = -99.0
     tmp_mass_buf  = 0.0

     ! Copy older compressed history into the newly expanded arrays
     tmp_resid_buf(:, 1:nr_n_cells) = nr_resid_buf(:, 1:nr_n_cells)
     tmp_mass_buf(:, 1:nr_n_cells)  = nr_mass_buf(:, 1:nr_n_cells)

     call move_alloc(tmp_resid_buf, nr_resid_buf)
     call move_alloc(tmp_mass_buf, nr_mass_buf)
     nr_n_cells = s%nz
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


! history of residuals for each structure equation
subroutine equ_resid_hist(id, ierr)
  use star_pgstar
  integer, intent(in)  :: id
  integer, intent(out) :: ierr

  type(pgstar_win_file_data), pointer :: p
  type(star_info), pointer :: s

  ierr = 0
  call star_ptr(id, s, ierr)
  if (ierr /= 0) return

  p => s%pg%pgstar_win_file_ptr(i_Other)
  p%plot              => equ_resid_plot
  p%id                =  i_Other
  p%name              = 'Max Residual per equation across mesh points'
  p%win_flag          = .true.
  p%win_width         =  12.0
  p%win_aspect_ratio  =  0.6
  p%file_flag         = .true.
  p%file_dir          = 'png'         ! folder where to save files
  p%file_prefix       = 'eq_resid_'
  p%file_interval     =  5
  p%file_width        = -1.0
  p%file_aspect_ratio = -1.0

end subroutine equ_resid_hist


subroutine equ_resid_plot(id, device_id, ierr)

  use star_pgstar

  implicit none

  integer, intent(in)  :: id
  integer, intent(in)  :: device_id
  integer, intent(out) :: ierr

  type(star_info), pointer :: s

  integer :: i_eq
  integer :: n, last_model_plotted

  real :: xmin, xmax
  real :: ymin, ymax

  real, dimension(max_resid_hist) :: xvec
  real, dimension(max_resid_hist) :: yvec

  ierr = 0

  call star_ptr(id, s, ierr)
  if (ierr /= 0) return

  if (.not. allocated(resid_hist_model)) then
     resid_hist_nvar = s%nvar_hydro
     allocate(resid_hist_model(max_resid_hist))
     allocate(resid_hist_vals(resid_hist_nvar, max_resid_hist))
     allocate(resid_equ_names(resid_hist_nvar))

     do i_eq = 1, resid_hist_nvar
        resid_equ_names(i_eq) = trim(s%nameofequ(i_eq))
     end do

     n_resid_hist = 0
  end if

  if (s%model_number /= last_model_plotted) then
     last_model_plotted = s%model_number
     if (n_resid_hist < max_resid_hist) then
        n_resid_hist = n_resid_hist + 1
     else
        resid_hist_model(1:max_resid_hist-1) = resid_hist_model(2:max_resid_hist)
        resid_hist_vals(:,1:max_resid_hist-1) = resid_hist_vals(:,2:max_resid_hist)
     end if

     resid_hist_model(n_resid_hist) = s%model_number

     do i_eq = 1, resid_hist_nvar
        resid_hist_vals(i_eq,n_resid_hist) = &
             real(max(1d-40, &
             maxval(abs(s%equ(i_eq,1:s%nz)))))
     end do

  end if

  n = n_resid_hist
  if (n < 2) return ! nothing to plot yet

  xvec(1:n) = real(resid_hist_model(1:n))

  xmin = xvec(1)
  xmax = xvec(n)

  if (xmin >= xmax) xmax = xmin + 1.0

  ymin = minval(resid_hist_vals(:,1:n))
  ymax = maxval(resid_hist_vals(:,1:n))
  ymin = max(ymin, 1e-20)
  ymax = max(ymax, 10.0*ymin)

  call pgslct(device_id)
  call pgsave
  call pgeras
  call pgsvp(0.12, 0.88, 0.12, 0.88)
  call pgswin(xmin, xmax, log10(ymin), log10(ymax))
  call pgscf(1)
  call pgsch(1.0)
  call pgsci(1)
  call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
  call pgmtxt('B',2.5,0.5,0.5,'Model Number')
  call pgmtxt('L',3.0,0.5,0.5, 'log10(max |residual|)')
  call pgmtxt('T',1.0,0.5,0.5,'Maximum structural equation residuals')

  ! draw lines
  do i_eq = 1, resid_hist_nvar
     yvec(1:n) = log10(resid_hist_vals(i_eq,1:n))
     call pgsci(mod(i_eq-1,13)+2)
     call pgline(n, xvec, yvec)
  end do

  ! legend
  call pgsch(0.65)

  do i_eq = 1, resid_hist_nvar
     call pgsci(mod(i_eq-1,13)+2)
     call pgptxt( &
          xmin + (0.55 + 0.06*real(i_eq-1))*(xmax-xmin), &
          log10(ymax) - 0.03*(log10(ymax)-log10(ymin)), &
          0.0, 0.0, &
          trim(resid_equ_names(i_eq)))
  end do

  call pgsci(1)
  call pgunsa

end subroutine equ_resid_plot
