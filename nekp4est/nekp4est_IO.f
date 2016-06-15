!> @file nekp4est_IO.f
!! @ingroup nekp4est
!! @brief set of I/O routines for nekp4est
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Get free unit number
!! @param[out] iunit     file unit
!! @param[out] ierr      error mark
!! @note Identical routine in Nek_framework/io/io_tools/io_tools.f
      subroutine io_file_freeid(iunit, ierr)
      implicit none

!     argument list
      integer iunit
      integer ierr

!     keeep track of max iunit generated
      integer io_iunit_min, io_iunit_max
      common /io_iunit/ io_iunit_min, io_iunit_max

!     local variables
      logical ifcnnd            ! is unit connected
!-----------------------------------------------------------------------
!     initialise variables
      ierr=0
      iunit = io_iunit_min

      do
         inquire(unit=iunit,opened=ifcnnd,iostat=ierr)
         if(ifcnnd) then
            iunit = iunit +1
         else
            exit
         endif
      enddo

      if (iunit.gt.io_iunit_max) io_iunit_max = iunit

      return
      end
!=======================================================================
!> @brief Close opened files
!! @note Identical routine in Nek_framework/io/io_tools/io_tools.f
      subroutine io_file_close()
      implicit none

!     keeep track of max iunit generated
      integer io_iunit_min, io_iunit_max
      common /io_iunit/ io_iunit_min, io_iunit_max

!     local variables
      integer iunit, ierr
      logical ifcnnd            ! is unit connected
!-----------------------------------------------------------------------
      do iunit = io_iunit_min, io_iunit_max
         inquire(unit=iunit,opened=ifcnnd,iostat=ierr)
         if(ifcnnd) close(iunit)
      enddo
      io_iunit_max = io_iunit_min

      return
      end
!=======================================================================
!     following subroutines are modiffications of
!
!     mfo_open_files, mfo_outfld from prepost.f;
!     mfi from ic.f
!=======================================================================
!> @brief Generate file name according to nek rulles
!! @param[out]  fname     file name
!! @param[in]   bname     base name
!! @param[in]   prefix    prefix
!! @param[out]  ierr      error mark
!! @note Identical routine in Nek_framework/io/io_tools/io_tools.f
      subroutine io_mfo_fname(fname,bname,prefix,ierr)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IFREGUO
      include 'RESTART_DEF'
      include 'RESTART'         ! NFILEO

!     argument list
      character*132  fname, bname
      character*3 prefix
      integer ierr

!     local variables
      integer ndigit, itmp
      real rfileo

      character*6  six
      save         six
      data         six / "??????" /
!-----------------------------------------------------------------------
!     initialise variables
      ierr = 0
      fname = ''

!     numbe or IO nodes
#ifdef MPIIO
      rfileo = 1
#else
      rfileo = nfileo
#endif
      ndigit = log10(rfileo) + 1

!     Add directory
      if (ifdiro) fname = 'A'//six(1:ndigit)//'/'

!     Add prefix
      if (prefix(1:1).ne.' '.and.prefix(2:2).ne.' '
     $    .and.prefix(3:3).ne.' ')
     $     fname = trim(fname)//trim(adjustl(prefix))

!     Add SESSION
      fname = trim(fname)//trim(adjustl(bname))

      if (ifreguo) fname = trim(fname)//'_reg'

!     test string length
      itmp = len_trim(fname)
      if (itmp.eq.0) then
         write(*,*) 'ERROR: io_mfo_fname; zero lenght fname.'
         ierr = 1
         return
      elseif ((itmp+ndigit+2+5).gt.132) then
         write(*,*) 'ERROR: io_mfo_fname; fname too long.'
         write(*,*) 'Fname: ',trim(fname)
         ierr = 1
         return
      endif

      !  Add file-id holder and .f appendix
      fname = trim(fname)//six(1:ndigit)//'.f'

      return
      end
!=======================================================================
!> @brief Write current mesh data (GLL points) to the file
!! @param[in]   prefix    prefix
!! @param[in]   fnumber   file number
      subroutine nekp4est_mfo(prefix, fnumber)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'GEOM_DEF'
      include 'GEOM'

!     argument list
      character*3 prefix
      integer fnumber

!     local variables
      integer il, ierr, ioflds, nout

      character*6  str

      integer lwdsizo
      logical lifxyo, lifpo, lifvo, lifto, lifreguo, lifpso(LDIMT1)

!     file strides
      integer*8 offs0,offs,stride,strideB,nxyzo8

      character*132  fname

      real tiostart, tio   ! simple timing
      real dnbyte          ! bandwidth count
!     functions
      real dnekclock_sync, glsum
!-----------------------------------------------------------------------
      tiostart=dnekclock_sync()

!     copy and set output parameters
      lwdsizo= WDSIZO
      WDSIZO = 8

      lifreguo= IFREGUO
      IFREGUO = .false.
      lifxyo= IFXYO
      IFXYO = .true.
      IFXYO_ = .true.
      lifpo= IFPO
      IFPO = .false.
      lifvo= IFVO
      IFVO = .false.
      lifto= IFTO
      IFTO = .false.
      do il=1,LDIMT1
         lifpso(il)= IFPSO(il)
         IFPSO(il) = .false.
      enddo

      nout = NELT
      NXO  = NX1
      NYO  = NY1
      NZO  = NZ1

!     set offset
      offs0 = iHeaderSize + 4 + isize*nelgt

      ierr = 0
      if (NID.eq.PID0) then         ! open files on i/o nodes
!     create file name
         call io_mfo_fname(fname,SESSION,prefix,ierr)
!     file number
         if (ierr.eq.0) then
            write(str,'(i5.5)') fnumber
            fname = trim(fname)//trim(str)
            call mbyte_open(fname,fid0,ierr)
         endif
      endif

      call err_chk(ierr,'Error opening file in nekp4est_mfo. Abort. $')
      call bcast(ifxyo_,lsize)
      ifxyo = ifxyo_
      call mfo_write_hdr                     ! create element mapping +
                                             ! write hdr
      nxyzo8  = NXO*NYO*NZO
      strideB = nelB * nxyzo8*WDSIZO
      stride  = nelgt* nxyzo8*WDSIZO

      ioflds = 0
      ! dump all fields based on the t-mesh to avoid different
      ! topologies in the post-processor
      ! write mesh only
      offs = offs0 + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)
      call mfo_outv(XM1,YM1,ZM1,nout,NXO,NYO,NZO)
      ioflds = ioflds + NDIM

      dnbyte = 1.*ioflds*nout*WDSIZO*NXO*NYO*NZO

!     put output variables back
      WDSIZO = lwdsizo

      IFREGUO = lifreguo
      IFXYO = lifxyo
      IFXYO_ = lifxyo
      IFPO = lifpo
      IFVO = lifvo
      IFTO = lifto
      do il=1,LDIMT1
         IFPSO(il) = lifpso(il)
      enddo

      ierr = 0
      if (NID.eq.PID0)
#ifdef MPIIO
     &   call byte_close_mpi(ifh_mbyte,ierr)
#else
     &   call byte_close(ierr)
#endif
      call err_chk(ierr,'Error closing file in nekp4est_mfo. Abort. $')

      tio = dnekclock_sync()-tiostart
      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4. +ISIZE*NELGT
      dnbyte = dnbyte/1024/1024
      if(NID.eq.0) write(6,7)  ISTEP,TIME,dnbyte,dnbyte/tio,
     &             NFILEO
    7 format(/,i9,1pe12.4,' done :: Write checkpoint',/,
     &       30X,'file size = ',3pG12.2,'MB',/,
     &       30X,'avg data-throughput = ',0pf7.1,'MB/s',/,
     &       30X,'io-nodes = ',i5,/)

      return
      end
!=======================================================================
!> @brief Read current mesh data (GLL points) from the file
!! @param[in]   prefix    prefix
!! @param[in]   fnumber   file number
!! @remarks This routine uses global scratch space SCRNS
      subroutine nekp4est_mfi(prefix, fnumber)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'RESTART_DEF'
      include 'RESTART'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'GEOM_DEF'
      include 'GEOM'

!     argument list
      character*3 prefix
      integer fnumber

!     local variables
      character*132  fname

      character*6  str

      integer iofldsr, ierr

!     file strides
      integer*8 offs0,offs,stride,strideB,nxyzr8

      real tiostart, tio   ! simple timing
      real dnbyte          ! bandwidth count
      integer*8 nbyte

!     scratch arrays
      integer lwk
      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      real wk(lwk)
      common /scrns/ wk

!     functions
      real dnekclock, glsum
!-----------------------------------------------------------------------
      tiostart=dnekclock()

!     create file name
      call io_mfo_fname(fname,SESSION,prefix,ierr)
      call err_chk(ierr,'Error opening file, in nekp4est_mfi.$')
!     file number
      write(str,'(i5.5)') fnumber
      fname = trim(fname)//trim(str)

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping

      offs0   = iHeadersize + 4 + ISIZE*NELGR
      nxyzr8  = NXR*NYR*NZR
      strideB = nelBr* nxyzr8*WDSIZR
      stride  = nelgr* nxyzr8*WDSIZR

!     read arrays
      iofldsr = 0
!     resid array
      offs = offs0 + NDIM*strideB
      call byte_set_view(offs,ifh_mbyte)
      call mfi_getv(XM1,YM1,ZM1,wk,lwk,.false.)

      iofldsr = iofldsr + NDIM

      nbyte = 0
      if(nid.eq.pid0r) nbyte = iofldsr*nelr*wdsizr*nxr*nyr*nzr

      ierr = 0

c     close files
#ifdef MPIIO
      if (nid.eq.pid0r) call byte_close_mpi(ifh_mbyte,ierr)
#else
      if (nid.eq.pid0r) call byte_close(ierr)
#endif

      call err_chk(ierr,'Error closing file, in nekp4est_mfi.$')

      tio = dnekclock()-tiostart

      dnbyte = nbyte
      nbyte = glsum(dnbyte,1)
      nbyte = nbyte + iHeaderSize + 4 + isize*nelgr

      if(nid.eq.0) write(6,7) istep,time,
     &             nbyte/tio/1024/1024/10,
     &             nfiler
    7 format(/,i9,1pe12.4,' done :: Read checkpoint data',/,
     &       30X,'avg data-throughput = ',f7.1,'MBps',/,
     &       30X,'io-nodes = ',i5,/)

      return
      end
!=======================================================================
!> @brief Load tree information from the file
      subroutine nekp4est_mesh_load()
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     local variables
      real t1, tmp    ! simple timing

!     functions
      real dnekclock
!-----------------------------------------------------------------------
!     simple timing
      t1 = dnekclock()

!     load GLL points from MSH file
      if (NP4_IOSTART.le.0)
     $    call nekp4est_abort('Error; nekprest_mfi: wrong file number')
      call nekp4est_mfi('MSH',NP4_IOSTART)
!     fill XC, YC, ZC arrays
!!!!      call nekp4est_fillxyzc

!     simple timing
      tmp= dnekclock() -t1

!     timers
!     total
      NP4_TC = NP4_TC + tmp
!     initialisation
      NP4_TCI = NP4_TCI + tmp

      return
      end
!=======================================================================
