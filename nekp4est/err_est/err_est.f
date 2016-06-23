!> @file err_est.f
!! @ingroup nekp4est
!! @brief Spectral error estimator.
!! @details Set of subroutines calculating error estimator based on
!! variable spectra. Adopted from Catherine Mavriplis code.
!! @author Adam Peplinski
!! @author Nicolas Offermans
!! @date Jun 20, 2016
!=======================================================================
!> @brief Get refinement mark; main interface
!! @param[out]   ref_mark  refinement mark; refine (1), coarsen (-1) or leave unchanged (0)
!! @param[in]    ref_level element refinement level
      subroutine err_est_mark(ref_mark, ref_level)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKP4EST'
      include 'ERR_EST'

!     for tests with advected cone
      include 'GEOM_DEF'
      include 'GEOM'

!     argument list
      integer ref_mark(LELT)
      integer ref_level(LELT)

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl, kl, ll ! loop index

!     min and max values of error estimator
      real eest_min, eest_max

!     work arrays
      real xa(LX1,LY1,LZ1), xb(LX1,LY1,LZ1)

!     initalisation
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./

!     for testing
      logical ifdone
      save ifdone
      data ifdone /.FALSE./

!     for tests with advected cone
!      real xloc, yloc, zloc, phi_loc, om_loc, r_loc, dr
!      save xloc, yloc, zloc, phi_loc, om_loc, r_loc, dr
!      real x, y, z,rmin,rmax
!     simple testing; element list frim the file
!      integer iunit, ierr
!      integer lref_elements
!      parameter (lref_elements=30)
!      integer ref_elements(lref_elements), n_ref_elem
!      save ref_elements, n_ref_elem

!     functions
      real glmax, glmin
!-----------------------------------------------------------------------
      logs = 'Get error estimate.'
      call nekp4est_log(NP4_LP_PRD,logs)

!     reset refinement mark
      call izero(ref_mark,LELT)

!     this is ont the best way of initialisation, but for now I leave it here
      if(.not.ifcalled) then
        ifcalled=.TRUE.
!     initialise error estimator
        call err_est_init

!     for testing
!     refine blocks from element_lst.txt file
!     read the file
!        ierr = 0
!        if(NID.eq.0) then
!           call io_file_freeid(iunit, ierr)
!           if (ierr.eq.0) then
!              open(unit=iunit,file='element_lst.txt',status='old',
!     $        action='read',iostat=ierr)
!              if (ierr.eq.0) then
!                 read(iunit,*) n_ref_elem
!                 if (n_ref_elem.le.lref_elements) then
!                    do il=1,n_ref_elem
!                       read(iunit,*) ref_elements(il)
!                    enddo
!                 else
!                    ierr = 1
!                 endif
!                 close(iunit)
!              endif
!           endif
!        endif
!        call nekp4est_chk_abort(ierr,'ERROR reading element_lst.txt')
!        call bcast(n_ref_elem,ISIZE)
!        call bcast(ref_elements,lref_elements*ISIZE)

!     for tests with advected cone
!        xloc = 0.7 - 0.5
!        yloc = 0.3 - 0.5
!        zloc = 0.5
!        r_loc = sqrt(xloc*xloc +yloc*yloc)
!        xloc = xloc/r_loc
!        yloc = yloc/r_loc
!        phi_loc = atan2(yloc,xloc)
!        om_loc = 1.0!*ATAN(1.0)
!        dr = 0.05/(2.0**(int(PARAM(35))-1))
      endif

!     for tests with advected cone
!     sphere and cylinder xy
!        xloc = 0.5 + r_loc*cos(om_loc*time+ phi_loc)
!        yloc = 0.5 + r_loc*sin(om_loc*time+ phi_loc)
!     cylinder xz
!        xloc = 0.5 + r_loc*cos(om_loc*time+ phi_loc)
!        zloc = 0.5 + r_loc*sin(om_loc*time+ phi_loc)
!        yloc = 0.5

!     for testing
!     refine/derefine whole domain
!      if (.not.ifdone) then
!        ifdone = .TRUE.
!        call ifill(ref_mark,1,LELT)
!      else
!        ifdone = .FALSE.
!        call ifill(ref_mark,-1,LELT)
!      endif
!
!      return

!     for testing
!     refine blocks from element_lst.txt file
!      if (.not.ifdone) then
!        ifdone = .TRUE.
!         do jl=1,NELT
!            do kl=1,n_ref_elem
!                if (LGLEL(jl).eq.ref_elements(kl)) ref_mark(jl) = 1
!            enddo
!         enddo
!      else
!        do jl=1,NELT
!            if (ref_level(jl).gt.0) ref_mark(jl) = -1
!        enddo
!        ifdone = .FALSE.
!      endif
!
!      return



!     Run error estimater for a given set of variables and mark
!     elelments according to refinement/derefinement threshold.
!     In general one should be able to use more than one variable.
!     That is why we define EEST_IFESTV logical array, but we should be
!     carreful and do not overwrite previous refine marks

!     velocity components
      if (IFFLOW) then
!     X component
        if (EEST_IFESTV(1)) then
            call err_est_var(EEST_EST,EEST_SIG,VX,NELV,xa,xb)
            do il=1,NELV
                if(EEST_EST(il).gt.EEST_REFT) then
                    ref_mark(il) = 1
                else if(EEST_EST(il).lt.EEST_DREFT.and.
     $                  ref_mark(il).ne.1) then
                    ref_mark(il) = -1
                endif
            enddo
!     stamp the file
            eest_min = glmin(EEST_EST,NELV)
            eest_max = glmax(EEST_EST,NELV)
            if (NID.eq.0) then
                write(*,*) 'EEST, velocity X'
                write(*,*) 'eest_min = ',eest_min
                write(*,*) 'eest_max = ',eest_max
            endif
        endif

!     Y component
        if (EEST_IFESTV(2)) then
            call err_est_var(EEST_EST,EEST_SIG,VY,NELV,xa,xb)
            do il=1,NELV
                if(EEST_EST(il).gt.EEST_REFT) then
                    ref_mark(il) = 1
                else if(EEST_EST(il).lt.EEST_DREFT.and.
     $                  ref_mark(il).ne.1) then
                    ref_mark(il) = -1
                endif
            enddo
!     stamp the file
            eest_min = glmin(EEST_EST,NELV)
            eest_max = glmax(EEST_EST,NELV)
            if (NID.eq.0) then
                write(*,*) 'EEST, velocity Y'
                write(*,*) 'eest_min = ',eest_min
                write(*,*) 'eest_max = ',eest_max
            endif
        endif

!     Z component
        if (IF3D.and.EEST_IFESTV(3)) then
            call err_est_var(EEST_EST,EEST_SIG,VZ,NELV,xa,xb)
            do il=1,NELV
                if(EEST_EST(il).gt.EEST_REFT) then
                    ref_mark(il) = 1
                else if(EEST_EST(il).lt.EEST_DREFT.and.
     $                  ref_mark(il).ne.1) then
                    ref_mark(il) = -1
                endif
            enddo
!     stamp the file
            eest_min = glmin(EEST_EST,NELV)
            eest_max = glmax(EEST_EST,NELV)
            if (NID.eq.0) then
                write(*,*) 'EEST, velocity Z'
                write(*,*) 'eest_min = ',eest_min
                write(*,*) 'eest_max = ',eest_max
            endif
        endif

      endif

!     temperature and passive scalars
      if (IFHEAT) then
         do il=2,NFIELD
            if (EEST_IFESTV(il+2)) then
!     check T versus V mesh
                call err_est_var(EEST_EST,EEST_SIG,
     $           T(1,1,1,1,il-1),NELFLD(il),xa,xb)
                do jl=1,NELFLD(il)
                    if(EEST_EST(jl).gt.EEST_REFT) then
                        ref_mark(jl) = 1
                    else if(EEST_EST(jl).lt.EEST_DREFT.and.
     $                  ref_mark(jl).ne.1) then
                        ref_mark(jl) = -1
                    endif
                enddo
!     stamp the file
                eest_min = glmin(EEST_EST,NELFLD(il))
                eest_max = glmax(EEST_EST,NELFLD(il))
                if (NID.eq.0) then
                    write(*,*) 'EEST, Temp/PS, ', il-1
                    write(*,*) 'eest_min = ',eest_min
                    write(*,*) 'eest_max = ',eest_max
                endif

!     for tests with advected cone
!     overwrite initial values and mark different positions
!            if (1.eq.0) then
!                call ifill(ref_mark,-1,LELT)
!                kl=nx1*ny1*nz1
!                do jl=1,NELFLD(il)
!                    rmin = 1.0E10
!                    rmax = 0.0
!                    do ll=kl*(jl-1)+1,kl*jl
!                        x = xloc - XM1(ll,1,1,1)
!                        y = yloc - YM1(ll,1,1,1)
!                        if (IF3D) then
!                            z = zloc - ZM1(ll,1,1,1)
!                        else
!                            z = 0.0
!                        endif
!     sphere
!                        x = sqrt(x*x+y*y+z*z)
!     cylinder
!     xy
!                        x = sqrt(x*x+y*y)
!     cylinder
!     xz
!                        x = sqrt(x*x+z*z)
!                        rmin = min(rmin,x)
!                        rmax = max(rmax,x)
!                    enddo
!                    if (rmin.lt.dr) ref_mark(jl) = 1
!                    if (rmin.lt.(0.1+dr).and.
!     &                  rmax.gt.(0.1-dr)) ref_mark(jl) = 1
!                enddo
!            endif
!     test; end

            endif
         enddo                  ! passive scalar loop
      endif

!     Place for MHD


      return
      end
!=======================================================================
!| @brief Initialise error estimator
      subroutine err_est_init
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'
      include 'ERR_EST'

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl
!-----------------------------------------------------------------------
!     stamp logs
      logs = 'Error est start.'
      call nekp4est_log(NP4_LP_PRD,logs)

!     set cutoff parameters
!     used for values
      EEST_SMALL = 1.e-14
!     used for ratios
      EEST_SMALLR = 1.e-10
!     used for gradients
      EEST_SMALLG = 1.e-5
!     used for sigma and rtmp in error calculations
      EEST_SMALLS = 0.2

!     refinement and derefinement thresholds
!     This should be taken from .rea file, but for now I just set it here.
      EEST_REFT  = 1.0e-6
      EEST_DREFT = 5.0e-7

!     number of points in fitting
      EEST_NP = 4
!     last modes skipped
      EEST_ELR = 0

!     correctness check
      if (EEST_NP.gt.EEST_NP_MAX) then
         logs = 'EEST_NP greater than EEST_NP_MAX.'
         call nekp4est_abort(logs)
      endif
      il = EEST_NP+EEST_ELR
      jl = min(LX1,LY1)
      if (IF3D) jl = min(jl,LZ1)
      if (il.gt.jl) then
         logs = 'EEST_NP+EEST_ELR greater than L?1'
         call nekp4est_abort(logs)
      endif

!     initalise coefficient mapping
      call err_est_cff_init

!     logical flags for variables that will be taken into account
!     velx - 1
!     vely - 2
!     velz - 3
!     temp - 4
!     ps   - 5...
      do il =1,LDIMT3
        EEST_IFESTV(il) = .FALSE.
      enddo

!     for this setup only temperature is tested
      EEST_IFESTV(4) = .TRUE.

      return
      end
!=======================================================================
!> @brief Initialise spectral coefficient mapping
      subroutine err_est_cff_init
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'WZ_DEF'
      include 'WZ'
      include 'NEKP4EST'
      include 'ERR_EST'

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl, kl, nn
!     Legendre polynomial
      real plegx(LX1),plegy(LY1),plegz(LZ1)
      real z, rtmp
!-----------------------------------------------------------------------
!     check polynomial order and numer of points for extrapolation
      if (min(NX1,NY1).le.(EEST_NP+EEST_ELR)) then
        logs = 'Error: increase L[XYZ]1'
        call nekp4est_abort(logs)
      endif

      if (IF3D.and.(NZ1.le.(EEST_NP+EEST_ELR))) then
        logs = 'Error: increase L[XYZ]1'
        call nekp4est_abort(logs)
      endif

!     initialise arrays
!     X - direction
      nn = NX1-1
      do jl= 1, NX1
!     Legendre polynomial
        z = ZGM1(jl,1)
        call legendre_poly(plegx,z,nn)
        do il=1, NX1
            EEST_XMAP(il,jl) = plegx(il)*WXM1(jl)
        enddo
      enddo
!     Y - direction; transposed
      nn = NY1-1
      do jl= 1, NY1
!     Legendre polynomial
        z = ZGM1(jl,2)
        call legendre_poly(plegy,z,nn)
        do il=1, NY1
            EEST_YTMAP(jl,il) = plegy(il)*WYM1(jl)
        enddo
      enddo
!     Z - direction; transposed
      if(IF3D) then
        nn = NZ1-1
        do jl= 1, NZ1
!     Legendre polynomial
            z = ZGM1(jl,3)
            call legendre_poly(plegz,z,nn)
            do il=1, NZ1
                EEST_ZTMAP(jl,il) = plegz(il)*WZM1(jl)
            enddo
        enddo
      endif

!     multiplicity factor
      rtmp = 1.0/2.0**NDIM
      do kl=1,NZ1
        do jl=1,NY1
            do il=1,NX1
                EEST_FAC(il,jl,kl) = (2.0*(il-1)+1.0)*(2.0*(jl-1)+1.0)*
     $             (2.0*(kl-1)+1.0)*rtmp
            enddo
        enddo
      enddo

      return
      end
!=======================================================================
!> @brief Get squared polynomial coefficients for error estimator in single element
!! @param[out]  coeff   modal coefficients
!! @param[in]   var     nodal coefficients
!! @param[in]   xa      work array
!! @param[in]   xb      work array
      subroutine err_est_el_lget(coeff,var,xa,xb)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ERR_EST'

!     argument list
      real coeff(LX1,LY1,LZ1), var(LX1,LY1,LZ1)
      real xa(LX1,LY1,LZ1), xb(LX1,LY1,LZ1)

!     local variables
      integer nxy, nyz, iz
!-----------------------------------------------------------------------
      nxy = NX1*NY1
      nyz = NY1*NZ1

!     test
!      call copy(coeff,EEST_YTMAP,nxy)
!      return

      if (IF3D) then
         call mxm(EEST_XMAP,NX1,var,NX1,xa,nyz)
         do iz = 1,NZ1
            call mxm(xa(1,1,iz),NX1,EEST_YTMAP,NY1,
     $           xb(1,1,iz),NY1)
         enddo
         call mxm(xb,nxy,EEST_ZTMAP,NZ1,coeff,NZ1)
      else
         call mxm(EEST_XMAP,NX1,var,NX1,xa,nyz)
         call mxm(xa,NX1,EEST_YTMAP,NY1,coeff,NY1)
      endif

!     multiply by factor
      iz = NX1*NY1*NZ1
      call col2(coeff,EEST_FAC,iz)

!     square coefficients
      call vsq(coeff, iz)

      return
      end
!=======================================================================
!> @brief Error and sigma for single variable
!! @details Get error estimater and sigma for single variable on a whole
!!  mesh for all directions.
!! @param[out]  est   estimated error
!! @param[out]  sig   estimated exponent
!! @param[in]   var   tested variable
!! @param[in]   nell  element number
!! @param[in]   xa    work array
!! @param[in]   xb    work array
      subroutine err_est_var(est,sig,var,nell,xa,xb)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'ERR_EST'

 !    argument list
      real est(nell), sig(nell)
      real var(LX1,LY1,LZ1,nell)
      integer nell
      real xa(LX1,LY1,LZ1), xb(LX1,LY1,LZ1)

!     local variables
      integer il, jl, kl, ll, j_st, j_en
!     polynomial coefficients
      real coeff(LX1,LY1,LZ1)
!     Legendre coefficients; first value coeff(1,1,1)
      real coef11
!     copy of last EEST_NP columns of coefficients
      real coefx(EEST_NP_MAX,LY1,LZ1),coefy(EEST_NP_MAX,LX1,LZ1),
     $     coefz(EEST_NP_MAX,LX1,LY1)
!     estimated error
      real estx, esty, estz
!     estimated decay rate
      real sigx, sigy, sigz
      real third
      parameter (third = 1.0/3.0)

!     for testing only
!      integer ntt
!      real ncn_test(lx1,ly1,lz1,lelt,10)
!      common /test_ncn/ ncn_test
!-----------------------------------------------------------------------
!     loop over elements
      do il = 1,nell
!     get square of polynomial coefficients for given variable
        call err_est_el_lget(coeff(1,1,1),var(1,1,1,il),xa,xb)

!     lower left corner
        coef11 = coeff(1,1,1)

!     small value; nothing to od
        if (coef11.ge.EEST_SMALL) then

!     extrapolate coefficients
!     X - direction
!     copy last EEST_NP collumns (or less if NX1 is smaller)
!     EEST_ELR allows to exclude last row
            j_st = max(1,NX1-EEST_NP+1-EEST_ELR)
            j_en = max(1,NX1-EEST_ELR)
            do ll = 1,NZ1
                do kl = 1,NY1
                    do jl = j_st,j_en
                        coefx(j_en-jl+1,kl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

!     get extrapolated values
            call err_est_extrap(estx,sigx,coef11,coefx,
     $           j_st,j_en,NY1,NZ1)

!     Y - direction
!     copy last EEST_NP collumns (or less if NY1 is smaller)
!     EEST_ELR allows to exclude last row
            j_st = max(1,NY1-EEST_NP+1-EEST_ELR)
            j_en = max(1,NY1-EEST_ELR)
            do ll = 1,NZ1
                do kl = j_st,j_en
                    do jl = 1,NX1
                        coefy(j_en-kl+1,jl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

!     get extrapolated values
            call err_est_extrap(esty,sigy,coef11,coefy,
     $           j_st,j_en,NX1,NZ1)

            if (IF3D) then
!     Z - direction
!     copy last EEST_NP collumns (or less if NZ1 is smaller)
!     EEST_ELR allows to exclude last row
                j_st = max(1,NZ1-EEST_NP+1-EEST_ELR)
                j_en = max(1,NZ1-EEST_ELR)
                do ll = j_st,j_en
                    do kl = 1,NY1
                        do jl = 1,NX1
                            coefz(j_en-ll+1,jl,kl) = coeff(jl,kl,ll)
                        enddo
                    enddo
                enddo

!     get extrapolated values
                call err_est_extrap(estz,sigz,coef11,coefz,
     $              j_st,j_en,NX1,NY1)

!     average
                est(il) =  sqrt(estx + esty + estz)
                sig(il) =  third*(sigx + sigy + sigz)
            else
                est(il) =  sqrt(estx + esty)
                sig(il) =  0.5*(sigx + sigy)
            endif

        else
!     for testing
            estx = 0.0
            esty = 0.0
            estz = 0.0
            sigx = -1.0
            sigy = -1.0
            sigz = -1.0
!     for testing; end

            est(il) =  0.0
            sig(il) = -1.0
        endif

!     for testing only
!        ntt = lx1*ly1*lz1
!        call copy (ncn_test(1,1,1,il,1) ,coeff(1,1,1),ntt)
!        call cfill(ncn_test(1,1,1,il,2) ,coef11,ntt)
!        call cfill(ncn_test(1,1,1,il,3) ,est(i),ntt)
!        call cfill(ncn_test(1,1,1,il,4) ,sig(i),ntt)
!        call cfill(ncn_test(1,1,1,il,5) ,estx,ntt)
!        call cfill(ncn_test(1,1,1,il,6) ,sigx,ntt)
!        call cfill(ncn_test(1,1,1,il,7) ,esty,ntt)
!        call cfill(ncn_test(1,1,1,il,8) ,sigy,ntt)
!        call cfill(ncn_test(1,1,1,il,9) ,estz,ntt)
!        call cfill(ncn_test(1,1,1,il,10),sigz,ntt)
!     for testing only; end

      enddo

      return
      end
!=======================================================================
!> @brief Get extrapolated values of sigma and error estimator.
!! @details We assume coef(n) = c*exp(-sigma*n) and estiamte sigma and
!!  eror eest = sqrt(2*(coef(N)**2/(2*N+1+\int_(N+1)^\infty coef(n)**2/(n+1) dn)))
!! @param[out]    estx   estimated error
!! @param[out]    sigx   estimated exponent
!! @param[in]     coef11 Legendre coefficients; first base-mode amplitude (coeff(1,1,1))
!! @param[in]     coef   Legendre coefficients; last EEST_NP columns
!! @param[in]     ix_st  first mode in coef
!! @param[in]     ix_en  last mode in coef
!! @param[in]     nyl    number of modes in y direction (aray size)
!! @param[in]     nzl    number of modes in z direction (aray size)
      subroutine err_est_extrap(estx,sigx,coef11,coef,
     $           ix_st,ix_en,nyl,nzl)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'ERR_EST'

!     argument list
      integer ix_st,ix_en,nyl,nzl
!     Legendre coefficients; last EEST_NP columns
      real coef(EEST_NP_MAX,nyl,nzl)
!     Legendre coefficients; first value coeff(1,1,1)
      real coef11
!     estimated error and decay rate
      real estx, sigx

!     local variables
      integer il, jl, kl, ll  ! loop index
      integer nsigt, pnr, nzlt
      real sigt, smallr, cmin, cmax, cnm, rtmp, rtmp2, rtmp3
      real sumtmp(4), cffl(EEST_NP_MAX)
      real stmp, estt, clog, ctmp, cave, erlog
      logical cuse(EEST_NP_MAX)
!-----------------------------------------------------------------------
!     initial values
      estx =  0.0
      sigx = -1.0

!     relative cutoff
      smallr = coef11*EEST_SMALLR

!     number of points
      pnr = ix_en - ix_st +1

!     to few points to interpolate
!      if ((ix_en - ix_st).le.1) return

!     for averaging, initial values
      sigt = 0.0
      nsigt = 0

!     loop over all face points
      nzlt = max(1,nzl - EEST_ELR) !  for 2D runs
      do il=1,nzlt
!     weight
        rtmp3 = 1.0/(2.0*(il-1)+1.0)
        do jl=1,nyl - EEST_ELR

!     find min and max coef along single row
            cffl(1) = coef(1,jl,il)
            cmin = cffl(1)
            cmax = cmin
            do kl =2,pnr
                cffl(kl) = coef(kl,jl,il)
                cmin = min(cmin,cffl(kl))
                cmax = max(cmax,cffl(kl))
            enddo

!     are coefficients sufficiently big
            if((cmin.gt.0.0).and.(cmax.gt.smallr)) then
!     mark array position we use in iterpolation
                do kl =1,pnr
                    cuse(kl) = .TRUE.
                enddo
!     max n for polynomial order
                cnm = real(ix_en)

!     check if all the points should be taken into account
!     in original code by Catherine Mavriplis this part is written
!     for 4 points, so I place if statement first
                if (pnr.eq.4) then
!     should we neglect last values
                    if ((cffl(1).lt.smallr).and.
     &                  (cffl(2).lt.smallr)) then
                        if (cffl(3).lt.smallr) then
                            cuse(1) = .FALSE.
                            cuse(2) = .FALSE.
                            cnm = real(ix_en-2)
                        else
                            cuse(1) = .FALSE.
                            cnm = real(ix_en-1)
                        endif
                    else
!     should we take stronger gradient
                        if ((cffl(1)/cffl(2).lt.EEST_SMALLG).and.
     $                      (cffl(3)/cffl(4).lt.EEST_SMALLG)) then
                            cuse(1) = .FALSE.
                            cuse(3) = .FALSE.
                            cnm = real(ix_en-1)
                        elseif ((cffl(2)/cffl(1).lt.EEST_SMALLG).and.
     $                          (cffl(4)/cffl(3).lt.EEST_SMALLG)) then
                            cuse(2) = .FALSE.
                            cuse(4) = .FALSE.
                        endif
                    endif
                endif

!     get sigma for given face point
                do kl =1,4
                    sumtmp(kl) = 0.0
                enddo
!     find new min and count number of points
                cmin = cmax
                cmax = 0.0
                do kl =1,pnr
                    if(cuse(kl)) then
                        rtmp  = real(ix_en-kl)
                        rtmp2 = log(cffl(kl))
                        sumtmp(1) = sumtmp(1) +rtmp2
                        sumtmp(2) = sumtmp(2) +rtmp
                        sumtmp(3) = sumtmp(3) +rtmp*rtmp
                        sumtmp(4) = sumtmp(4) +rtmp2*rtmp
!     find new min and count used points
                        cmin = min(cffl(kl),cmin)
                        cmax = cmax + 1.0
                    endif
                enddo
!     decay rate along single row
                stmp = (sumtmp(1)*sumtmp(2) - sumtmp(4)*cmax)/
     $                 (sumtmp(3)*cmax - sumtmp(2)*sumtmp(2))
!     for averaging
                sigt = sigt + stmp
                nsigt = nsigt + 1

!     get error estimator depending on calculated decay rate
                estt = 0.0
                if (stmp.lt.EEST_SMALLS) then
                    estt = cmin
                else
!     get averaged constant in front of c*exp(-sig*n)
                    clog = (sumtmp(1)+stmp*sumtmp(2))/cmax
                    ctmp = exp(clog)
!     average exponent
                    cave = sumtmp(1)/cmax
!     check quality of approximation comparing is to the constant cave
                    do kl =1,2
                        sumtmp(kl) = 0.0
                    enddo
                    do kl =1,pnr
                        if(cuse(kl)) then
                            erlog = clog - stmp*real(ix_en-kl)
                            sumtmp(1) = sumtmp(1)+
     $                          (erlog-log(cffl(kl)))**2
                            sumtmp(2) = sumtmp(2)+
     $                          (erlog-cave)**2
                        endif
                    enddo
                    rtmp = 1.0 - sumtmp(1)/sumtmp(2)
                    if (rtmp.lt.EEST_SMALLS) then
                        estt = cmin
                    else
!     last coefficient is not included in error estimator
                        estt = ctmp/stmp*exp(-stmp*cnm)
                    endif
                endif
!     add contribution to error estimator; variable weight
                estx = estx + estt/(2.0*(jl-1)+1.0)*rtmp3
            endif  ! if((cmin.gt.0.0).and.(cmax.gt.smallr))
        enddo
      enddo
!     constant weight
!     Multiplication by 4 in 2D / 8 in 3D
!     Normalization of the error by the volume of the reference element
!     which is equal to 4 in 2D / 8 in 3D
!     ==> Both operations cancel each other
      estx = estx/(2.0*(ix_en-1)+1.0)

!     final everaging
!     sigt = 2*sigma so we divide by 2
      if (nsigt.gt.0) then
        sigx = 0.5*sigt/nsigt
      endif

      return
      end
!=======================================================================
