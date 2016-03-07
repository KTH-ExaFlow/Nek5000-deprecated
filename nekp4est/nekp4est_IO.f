!> @brief set of I/O routines for nekp4est
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Initialise files counters
      subroutine IO_file_init()
      implicit none

!     keeep track of max iunit generated
      integer IO_iunit_min, IO_iunit_max
      common /IO_iunit/ IO_iunit_min, IO_iunit_max

!     local variables
      logical ifcall            ! first call
      save ifcall
      data ifcall /.TRUE./
!-----------------------------------------------------------------------
!     set initial IO_iunit_min/max
      if (ifcall) then
        ifcall = .FALSE.
        IO_iunit_min = 200
        IO_iunit_max = IO_iunit_min
      endif

      return
      end
!=======================================================================
!> @brief Get free unit number
!! @parameter[out] iunit     file unit
!! @parameter[out] ierr      error mark
      subroutine IO_file_freeid(iunit, ierr)
      implicit none

!     argument list
      integer iunit
      integer ierr
!     keeep track of max iunit generated
      integer IO_iunit_min, IO_iunit_max
      common /IO_iunit/ IO_iunit_min, IO_iunit_max
!     local variables
      logical ifcnnd            ! is unit connected
!-----------------------------------------------------------------------
      ierr=0
!     to not interact with the whole nek5000 I/O
      iunit = IO_iunit_min

      do
         inquire(unit=iunit,opened=ifcnnd,iostat=ierr)
         if(ifcnnd) then
            iunit = iunit +1
         else
            exit
         endif
      enddo

      if (iunit.gt.IO_iunit_max) IO_iunit_max = iunit

      return
      end
!=======================================================================
!> @brief Close opened files
      subroutine IO_file_close()
      implicit none

!     keeep track of max iunit generated
      integer IO_iunit_min, IO_iunit_max
      common /IO_iunit/ IO_iunit_min, IO_iunit_max

!     local variables
      integer iunit, ierr
      logical ifcnnd            ! is unit connected
!-----------------------------------------------------------------------
      do iunit = IO_iunit_min, IO_iunit_max
         inquire(unit=iunit,opened=ifcnnd,iostat=ierr)
         if(ifcnnd) close(iunit)
      enddo

      return
      end
!=======================================================================
!     following subroutines are modiffications of
!
!     mfo_open_files, mfo_outfld from prepost.f;
!     mfi from ic.f
!=======================================================================
!> @brief Generate name according to nek rulles
!! @parameter[out]  fname     file name
!! @parameter[in]   prefix    prefix
!! @parameter[in]   bname     basename
      subroutine IO_mfo_fname(fname,bname,prefix)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'RESTART_DEF'
      include 'RESTART'

!     argument list
      character*132  fname, bname
      character*3 prefix

!     local variables
      integer ndigit
      real rfileo

      character*6  six
      save         six
      data         six / "??????" /

      character*1 slash,dot
      save        slash,dot
      data        slash,dot  / '/' , '.' /
!-----------------------------------------------------------------------
      fname = ''

!     numbe or IO nodes
#ifdef MPIIO
      rfileo = 1
#else
      rfileo = nfileo
#endif
      ndigit = log10(rfileo) + 1

!     Add directory
      if (ifdiro) fname = 'A'//six(1:ndigit)//slash

!     Add prefix
      if (prefix(1:1).ne.' '.and.prefix(2:2).ne.' '
     $    .and.prefix(3:3).ne.' ') fname = trim(fname)//trim(prefix)

!     Add SESSION
      fname = trim(fname)//trim(bname)

      if (ifreguo) fname = trim(fname)//'_reg'

      !  Add file-id holder and .f appendix
      fname = trim(fname)//six(1:ndigit)//'.f'

      return
      end
!=======================================================================
!> @brief Write current mesh data (GLL points) to the file
!! @parameter[in]   prefix    prefix
!! @parameter[in]   fnumber   file number
      subroutine nekp4est_mfo(prefix, fnumber)  ! muti-file output
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
         call IO_mfo_fname(fname,SESSION,prefix)
!     file number
         write(str,'(i5.5)') fnumber
         fname = trim(fname)//trim(str)
         call mbyte_open(fname,fid0,ierr)
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
!! @parameter[in]   prefix    prefix
!! @parameter[in]   fnumber   file number
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
      call IO_mfo_fname(fname,SESSION,prefix)
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
