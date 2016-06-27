!> @file nekp4est.f
!! @ingroup nekp4est
!! @brief Main interface for nekp4est
!! @author Adam Peplinski
!! @date Feb 26, 2016
!=======================================================================
!> @brief Initialisation of sc and p4est libraries
!! @param[in] intracomm mpi communicator
      subroutine nekp4est_init(intracomm)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'

!     argument list
      integer intracomm

!     local variables
!     simple timing
      real t1, tmp

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
      if (N_DIM.ne.LDIM) then
         call nekp4est_abort('Error: nekp4est ldim inconsistent')
      endif

      if (N_PSCL.ne.LDIMT) then
         call nekp4est_abort('Error: nekp4est ldimt inconsistent')
      endif

!     for now as interpolation for hanging faces does not support it yet
      if (IF3D) then
        if(mod(LX1,2).eq.1.or.mod(LY1,2).eq.1.or.mod(LZ1,2).eq.1) then
           call nekp4est_abort('Error: nekp4est LX1 must be even.')
        endif
      else
        if(mod(LX1,2).eq.1.or.mod(LY1,2).eq.1) then
           call nekp4est_abort('Error: nekp4est LX1 must be even.')
        endif
      endif

!     this is temporary
#ifdef MOAB
      call nekp4est_abort('Error: nekp4est does not support moab')
#endif

!     set reset flag
      NP4_IFRESET = .FALSE.

!     simple timing
      tmp = dnekclock() - t1
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
      character(len=NP4_LSTL_FNM) filename ! file name

!     simple timing
      real t1, tmp

!     functions
      integer ltrunc
      real dnekclock
!-----------------------------------------------------------------------
!     simple timing
      t1 = dnekclock()

!     save current tree structure
      call nekp4est_log(NP4_LP_PRD,'Saving forest data.')

!     prepare file name
      filename = trim(adjustl(SESSION))
      if (len_trim(filename).gt.(NP4_LSTL_FNM-12)) then
         call nekp4est_abort
     $        ('ERROR: nekp4est_end; too long output file name')
      endif
      filename = trim(filename)//'_end.tree'//CHAR(0)
      call fp4est_tree_save(1,filename)

!     save GLL points to MSH file
      if (NP4_IOSTOP.eq.NP4_IOSTART) then
         call nekp4est_log(NP4_LP_ESS,
     $        'Warning; nekp4est_mfo: equal input/output file numbers')
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
      tmp = dnekclock() - t1
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
!=======================================================================
!> @brief Load tree information from the file
      subroutine nekp4est_tree_load()
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'

!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     common blocks
      integer nidl,npl,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidl,npl,nekcomm,nekgroup,nekreal

!     local variables
      character(len=NP4_LSTL_FNM) filename ! file name

      integer itmp, lcbc, lrbc

!     simple timing
      real t1, tmp

!     functions
      integer iglmax, iglsum
      real dnekclock
      logical nekp4est_ifmesh
!-----------------------------------------------------------------------
!     simple timinng
      t1 = dnekclock()

!     read forest data
      call nekp4est_log(NP4_LP_PRD,'Reading forest data.')

!     prepare file name
      filename = trim(adjustl(SESSION))
      if (len_trim(filename).gt.(NP4_LSTL_FNM-8)) then
         call nekp4est_abort
     $        ('ERROR: nekp4est_tree_load; too long output file name')
      endif
      filename = trim(filename)//'.tree'//CHAR(0)
      call fp4est_tree_load(nekcomm, 1,filename)

!     get b.c. position in bc and cbc array
      IBC = 2
      if (IFFLOW) IBC = 1

      NFLDT = 1
      if (IFHEAT) NFLDT = 2+NPSCAL
      if (IFMHD ) NFLDT = 2+NPSCAL+1

!
!     If p32 = 0.1, there will be no bcs read in
!
      if (PARAM(32).gt.0) NFLDT = IBC + PARAM(32)-1

!     it must done afrer rdparam is executed
!     set max refinement level
      NP4_LMAX = int(abs(PARAM(35)))

!     create ghost zones
      call fp4est_ghost_new()

!     get new mesh
      call fp4est_mesh_new()

!     get mesh size
      call fp4est_msh_get_size(NP4_NELGT,NP4_NELIT,NP4_NELT,
     $    NP4_NELV,NP4_MLEV)

!     get global max quadrant level
      NP4_MLEV = iglmax(NP4_MLEV,1)

!     check if we generate mesh from .rea or .mesh
!!!      if (nekp4est_ifmesh()) then
      if (.TRUE.) then
!     we use nek5000 standard method to generate the mesh

!     check consistency of p4est structure and .rea file
!     T mesh
        if (NELGT.ne.NP4_NELGT) then
           call nekp4est_abort
     $          ('Error: nekp4est_tree_load; NELGT inconsistent')
        endif
!     V mesh
!     get globalnumber of V elements
        itmp = iglsum(NP4_NELV,1)
        if (NELGV.ne.itmp) then
           call nekp4est_abort
     $          ('Error: nekp4est_tree_load; NELGV inconsistent')
        endif

      else
!     we skip mesh generation part in nek5000 and mesh structure data in .rea

!     global element count
        NELGT = NP4_NELGT
!     get globalnumber of V elements
        NELGV = iglsum(NP4_NELV,1)

!     local element count related to p4est element distribution
        NELT = NP4_NELT
        NELV = NP4_NELV

!     stamp log
        if (nid.eq.0) then
            write(*,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
            write(*,12) 'lx1  /lx2  /lx3 :',lx1,lx2,lx3
 12         format(1X,A,4I12,/,/)
        endif

!     check array sizes for p4est element distribution
        call chk_nel

!     initialize arrays
!     some of them are done in initdat but not all. It would be good to do it consistently
!     for now I'm leaving it as it is
        call rzero(CURVE ,72*LELT)
        call blank(CCURVE,12*LELT)
        lcbc=18*LELT*(LDIMT1 + 1)
        lrbc=30*LELT*(LDIMT1 + 1)
        call rzero(BC ,lrbc)
        call blank(CBC,lcbc)

!     reset elemnt counts
        EL_COUNT = 0
        NP4_MAP_NR = 0
        NP4_RFN_NR = 0
        NP4_CRS_NR = 0

!     load mesh data to local arrays
!     in this version we load correct bc data and only mark curved faces
!     this should be enough as mesh generation is skipped and curvature is
!     used only in setrzer and setdef
!!!        call fp4est_msh_get_dat

!     new element distribution
!!!        call mapelpr

!     mesh transfer
!     important topology variables;
!     originally set in stup_topo, but I use them, so I do it here
!!!        call nekp4est_setup_topo

!     redistribute data
!!!        call nekp4est_tree_transfer

!     GLL points (XM1, Ym1, ZM1) and element vertices (XC, YC, ZC)
!     will be filled after io module is initialised
      endif

!$$$!     testing
!$$$      call fp4est_vtk_write('test')

!     simple timing
      tmp = dnekclock() - t1
!     total
      NP4_TC = NP4_TC + tmp
!     initialisation
      NP4_TCI = NP4_TCI + tmp

      return
      end
