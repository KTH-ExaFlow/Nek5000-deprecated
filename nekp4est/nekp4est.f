!> @brief Main interface for nekp4est
!! @author Adam Peplinski
!! @date Feb 26, 2016
!=======================================================================
!> @brief Initialisation of sc and p4est libraries
!! @parameter[in] intracomm mpi communicator
      subroutine nekp4est_init(intracomm)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'

!     input parameters
      integer intracomm

!     local variables
!     simple timing
      real t1, t2, tmp

!     functions
      real dnekclock
!-----------------------------------------------------------------------
!     simple timing
      t1 = dnekclock()

!     reset AMR timers
      NP4_TC = 0.0
      NP4_TCI = 0.0
      NP4_TCF = 0.0
      NP4_TCE = 0.0
      NP4_TCER = 0.0
      NP4_TCP4 = 0.0
      NP4_TCP5 = 0.0
      NP4_TCP6 = 0.0
      NP4_TCP7 = 0.0
      NP4_TCP8 = 0.0
      NP4_TCPN = 0.0
      NP4_TCPS = 0.0
      NP4_TCPD = 0.0
      NP4_TCPM = 0.0
      NP4_TCPP = 0.0
      NP4_TCPPR = 0.0
      NP4_TCPR = 0.0
      NP4_TCL = 0.0
      NP4_TCC = 0.0
      NP4_TCS = 0.0
      NP4_TCSC = 0.0
      NP4_TCN = 0.0

!     reset refinement/coarsening count
      NP4_REF_CNT = 0
      NP4_EER_CNT = 0

!     reset mapping count
      NP4_IMAP = 0

!     set initial default log level
      NP4_LP_DEF=NP4_LP_PRD
      call fsc_set_log(NP4_LP_DEF)

!     sc init
      call fsc_init(intracomm, 1, 0, NP4_LP_DEF)

!     p4est init
      call fp4est_init(NP4_LP_DEF)

      call nekp4est_log(NP4_LP_DEF,'Starting nekp4est logs.')
 
!     check consistency with SIZE
      if (N_DIM.ne.LDIM)
     $ call nekp4est_abort('Error: nekp4est ldim inconsistent')

      if (N_PSCL.ne.LDIMT)
     $call nekp4est_abort('Error: nekp4est ldimt inconsistent')

!     for now as interpolation for hanging faces does not support it yet
      if (IF3D) then
        if(mod(LX1,2).eq.1.or.mod(LY1,2).eq.1.or.mod(LZ1,2).eq.1)
     $call nekp4est_abort('Error: nekp4est LX1 must be even.')
      else
        if(mod(LX1,2).eq.1.or.mod(LY1,2).eq.1)
     $call nekp4est_abort('Error: nekp4est LX1 must be even.')
      endif

!     this is temporary
#ifdef MOAB
      call nekp4est_abort('Error: nekp4est does not support moab')
#endif

!     set reset flag
      NP4_IFRESET = .FALSE.

!     simple timing
      t2 = dnekclock()
      tmp = t2 - t1
!     total
      NP4_TC = NP4_TC + tmp
!     initialisation
      NP4_TCI = NP4_TCI + tmp
      return
      end
!=======================================================================
!> @brief Finalisation of sc and p4est in nekton
      subroutine nekp4est_end()
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'

!     local variables
!     file name
      integer ls
      character*132 filename
      character filename1(132)
      equivalence (filename,filename1)

!     simple timing
      real t1, t2, tmp

!     functions
      integer ltrunc
      real dnekclock
!-----------------------------------------------------------------------
!     simple timing
      t1 = dnekclock()

!     prepare file name
      call blank(filename,132)
      ls = ltrunc(session,132)
      call chcopy(filename1(1),session,ls)
      call chcopy(filename1(ls+1),'_end.tree',9)
      call chcopy(filename1(ls+10),CHAR(0),1)

!     save current tree structure
      call nekp4est_log(NP4_LP_PRD,'Saving forest data.')
!      call fp4est_tree_save(1,filename)

!     save GLL points to MSH file
      if (NP4_IOSTOP.eq.NP4_IOSTART) then
        if (NIO.eq.0) write(6,*) 'Warning; nekp4est_mfo: equal ',
     $      'input/output file numbers, adjusting'
        NP4_IOSTOP = NP4_IOSTOP +1
      endif
      call nekp4est_mfo('MSH',NP4_IOSTOP)

!     release memory
!     delete mesh and ghost cells
      call fp4est_mesh_del()
      call fp4est_ghost_del()
!     delete tree and connectivity
      call fp4est_tree_del()
      call fp4est_cnn_del()

!     sc end
      call fsc_pkg_print(NP4_LP_DEF)
      call fsc_finalize()

!     simple timing
      t2 = dnekclock()
      tmp = t2 - t1
!     total
      NP4_TC = NP4_TC + tmp
!     finalisation
      NP4_TCF = NP4_TCF + tmp

      if (NID.eq.0) then
         write(6,*) 
         write(6,'(A)') 'Finalized AMR:'
         write(6,'(A,i6,A)')    'Cumulative time over ',NP4_REF_CNT,
     $      ' calls'
         write(6,'(A,g13.5,A)') 'Total time in AMR: ',NP4_TC  ,' sec'
         write(6,'(A,g13.5,A)') '      initialisation: ',NP4_TCI ,' sec'
         write(6,'(A,g13.5,A)') '        finalisation: ',NP4_TCF ,' sec'
         write(6,'(A,g13.5,A)') '           evolution: ',NP4_TCE ,' sec'
         write(6,'(A,g13.5,A,i6,A)') '        error estimator: ',
     $      NP4_TCER,' sec; ',NP4_EER_CNT,' calls'
         write(6,'(A,g13.5,A)') '             p4est part: ',NP4_TCP4,
     $      ' sec'
         write(6,'(A,g13.5,A)') '                mesh/ghost: ',NP4_TCP5,
     $      ' sec'
         write(6,'(A,g13.5,A)') '           refinement fill: ',NP4_TCP6,
     $      ' sec'
         write(6,'(A,g13.5,A)') '            refine/balance: ',NP4_TCP7,
     $      ' sec'
         write(6,'(A,g13.5,A)') '                 partition: ',NP4_TCP8,
     $      ' sec'
         write(6,*)
         write(6,'(A,g13.5,A)') '        regenarate mesh: ',NP4_TCPN,
     $      ' sec'
         write(6,'(A,g13.5,A)') '                 mesh size: ',NP4_TCPS,
     $      ' sec'
         write(6,'(A,g13.5,A)') '        p4est-nek transfer: ',NP4_TCPD,
     $      ' sec'
         write(6,'(A,g13.5,A)') '           element mapping: ',NP4_TCPM,
     $      ' sec'
         write(6,'(A,g13.5,A)') '                ParMetis call (F): ',
     $      NP4_TCPP,' sec'
         write(6,'(A,g13.5,A)') '                ParMetis call (R): ',
     $      NP4_TCPPR,' sec'
         write(6,'(A,g13.5,A)') '            redistribution: ',NP4_TCPR,
     $      ' sec'
         write(6,*)
         write(6,'(A,g13.5,A)') '         ref/crs  local: ',NP4_TCL ,
     $      ' sec'
         write(6,'(A,g13.5,A)') '       transfer/sorting: ',NP4_TCC ,
     $      ' sec'
         write(6,'(A,g13.5,A)') '         solver restart: ',NP4_TCS ,
     $      ' sec'
         write(6,'(A,g13.5,A)') '               crs restart: ',NP4_TCSC,
     $      ' sec'
         write(6,'(A)') 'Averaged values:'
         NP4_TCE = NP4_TCE/real(max(1,NP4_EER_CNT))
         write(6,'(A,g13.5,A)') 'Averaged cycle AMR: ',NP4_TCE ,' sec'
         write(6,'(A,g13.5,A)') 'Averaged cycle nek: ',NP4_TCN ,' sec'
         write(6,*)
      endif

      return
      end
