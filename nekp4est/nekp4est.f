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
!! @todo Upgrade curvature data; check setdef and setzer.
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

      integer itmp, lcbc, lrbc, el, il

!     tmp array for curvature data
      integer crvl(6,LELT)

!     simple timing
      real t1, tmp

!     functions
      integer iglmax, iglsum, iglmin
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

!     get global max element level
      NP4_MLEV = iglmax(NP4_MLEV,1)

!     check if we generate mesh from .rea or .mesh
!!!      if (nekp4est_ifmesh()) then
      if (.FALSE.) then
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
        call izero(crvl ,6*LELT)
        call rzero(CURVE ,72*LELT)
        call blank(CCURVE,12*LELT)
        lcbc=18*LELT*(LDIMT1 + 1)
        lrbc=30*LELT*(LDIMT1 + 1)
        call rzero(BC ,lrbc)
        call blank(CBC,lcbc)

!     load mesh data to local arrays
!     In this version we load correct bc data and only mark for curved faces.
!     For now this should be enough as mesh generation is skipped and curvature is
!     used only in setrzer and setdef
        call fp4est_msh_get_dat
     $       (IBC,NFLDT,NELGV,LELT,IGROUP,NP4_LEVEL,crvl,BC,CBC)

!     curvature data
!     I mark deformed as 'D' without going into details
!     it should be fine as setrzer and setdef do simple check
!     if(CCURVE(,).ne.' ')
        do el = 1, NELT
           do il=1,6
              if (crvl(il,el).ne.0) CCURVE(il,el) = 'D'
           enddo
        enddo

!     reset refinement history
        call fp4est_msh_get_hst(NP4_MAP_NR,NP4_RFN_NR,NP4_RFN_NR,
     $                        NP4_GLGL_MAP,NP4_GLGL_RFN, NP4_GLGL_CRS)

!     new element distribution
!!!        call mapelpr

!     mesh transfer
!     important topology variables;
!     originally set in stup_topo, but I use them, so I do it here
!!!        call nekp4est_setup_topo

!     redistribute data
!!!        call nekp4est_tree_transfer

!     GLL points (XM1, YM1, ZM1) and element vertices (XC, YC, ZC)
!     will be filled after io module is initialised
      endif

!     get mesh topology
      call nekp4est_topol_get()

#ifdef DEBUG
!     testing
      call fp4est_vtk_write('mesh_test'//char(0))
#endif

!     simple timing
      tmp = dnekclock() - t1
!     total
      NP4_TC = NP4_TC + tmp
!     initialisation
      NP4_TCI = NP4_TCI + tmp

      return
      end
!=======================================================================
!> @brief Get element-processor mapping.
!! @details This routine calculate element-processor mapping storred in
!!  gllnid (PARALLEL). Originally it is provided by ###.map file, but
!!  AMR requires dynamical mesh partitioning. This routine has to be
!!  executed at least once to change p4est distribution to the correct one.
!!  It is used as well for element-processor mapping update after every
!!  mesh refinement/coarseing. We use last partitioning (not p4est one)
!!  assuming all new children are placed on the parent node.
!! @todo There are a few arrays scaling with processor number LP. It would
!!  be good to change it.
!! @todo Add weights to the nodes and faces; Add edges and vertices
!!  to the graph together with wieghts.
!! @todo Add graph partitioning for V and T meshes.
      subroutine nekp4est_get_map()
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKP4EST'

!     MPI communication datails
      integer nidl,npl,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidl,npl,nekcomm,nekgroup,nekreal

!     local variables

      integer eg, il, jl, kl !loop index
      integer itmp1, itmp2

      integer el_countl
      integer graph(4*LELT), graph_offset(LELT)
      integer offset(0:LP),offtmp(0:LP)
      integer order(LELT),sizes(2*LP), order2(LELT),sizes2(0:2*LP)

!     simple timing
      real t1, t2

!     functions
      real dnekclock

!#ifdef DEBUG
!     for testing
      character*2 str
      integer iunit, ierr
!#endif
!-----------------------------------------------------------------------
!     init array
      call izero(GLLNID,NELGT)

# ifdef NEKPMETIS
      call nekp4est_log(NP4_LP_PRD,'Runing ParMETIS.')
!     for now we run it for hydro only
      if (NELGV.eq.NELGT) then
         call fp4est_msh_get_graph(el_countl, graph, graph_offset)

         call izero(offset,NP+1)
         offset(NID+1) = el_countl
!     global sum to exchange data
!     offtmp is used as scratch
         call igop(offset,offtmp,'+  ',NP+1)

         do eg =1,NP
            offset(eg+1) = offset(eg+1) + offset(eg)
         enddo

!     save for future backward communication (fom nek5000 to p4est element distribution)
         call icopy(NP4_NELNID,offset,NP+1)

#ifdef DEBUG
!     testing
      call io_file_freeid(iunit, ierr)
      write(str,'(i2.2)') NID
      open(unit=iunit,file='graph.txt'//str)
      write(iunit,*) (offset(il),il=0,NP)
      do jl=1,offset(nid+1)-offset(NID)
         write(iunit,*) jl+offset(NID)-1,
     $   graph_offset(jl+1)-graph_offset(jl),
     $   (graph(il),il=graph_offset(jl)+1,graph_offset(jl+1))
      enddo
      close(iunit)
#endif

!     We have to check if this is the first partitioning or the remap.
!     If mapping is done for the first time full partitioning based on
!     p4est element distribution (not optimal) must be done but this can
!     be slow
        if (NP4_IMAP.eq.0) then     ! check mapping count
!     first partitioning
!     we use sub-domain processor coupling related to p4est partitioning so no
!     initialisation of GLLNID is necessary
            eg=NP
!     simple timing
            t1 = dnekclock()
            call fpmetis_part(nekcomm, offset,graph_offset, graph,
     $          eg, GLLNID(NP4_NELIT+1))
!     this should increase partitioning quality, but I don't see any improvement
!            call fpmetis_refine(nekcomm, offset,graph_offset, graph,
!     $          eg, GLLNID(NP4_NELIT+1))
            t2 = dnekclock()
            NP4_TCPP = NP4_TCPP + t2 - t1

        else ! remapping
!     use last optimal mapping and update it (should be faster)
!     fill old mapping of new elements
!     reset refinement/coarsening counters
            il=1
            jl=1
!     loop over local p4est elements
            do eg=1,NP4_NELT
!     global element number
                itmp1 = NP4_NELIT+eg
!     check what happend to the element
!     if no refinement/coarsening
                itmp2 = NP4_GLGL_MAP(eg)
!     otherwise
                if (itmp2.eq.0) then
                    if (il.le.NP4_RFN_NR.and.
     $                  NP4_GLGL_RFN(1,il).eq.itmp1) then
!     refinement
                        itmp2 = NP4_GLGL_RFN(2,il)
                        il = il+1
                    elseif (jl.le.NP4_CRS_NR.and.
     $                  NP4_GLGL_CRS(1,1,jl).eq.itmp1) then
!     coarsening
                        itmp2 = NP4_GLGL_CRS(2,1,jl)
                        jl = jl+1
                    else
!     error; some option had to be taken
                        call nekp4est_abort
     $                       ('Error: remap; wrong pointer')
                    endif
                endif
                GLLNID(itmp1) = NP4_GLLNID_O(itmp2)
            enddo

            eg=NP
!     simple timing
            t1 = dnekclock()
!     ParMetis itr parameter describing the ratio of inter-processor
!     communication time compared to data redistribution time.
            t2 = 1000.0
            call fpmetis_rpart(nekcomm, offset,graph_offset, graph,
     $          eg, GLLNID(NP4_NELIT+1),t2)
!     this should increase partitioning quality, but I don't see any improvement
!            call fpmetis_refine(nekcomm, offset,graph_offset, graph,
!     $          eg, GLLNID(NP4_NELIT+1))
            t2 = dnekclock()
            NP4_TCPPR = NP4_TCPPR + t2 - t1

        endif ! mapping count
#ifdef DEBUG
!     testing
      eg = offset(NID+1)-offset(NID)
      call fp4est_vtk_iscalar(GLLNID(NP4_NELIT+1),eg,
     $  "metis_part"//CHAR(0))
      call io_file_freeid(iunit, ierr)
      write(str,'(i2.2)') NID
      open(unit=iunit,file='graph_part.txt'//str)
      write(iunit,*) (offset(il),il=0,NP)
      do il=1,offset(NID+1)-offset(NID)
         write(iunit,*) il+offset(NID)-1,
     $   GLLNID(NP4_NELIT+il)
      enddo
      close(iunit)
#endif

!     global sum to exchange data
!     GLLEL is used as scratch
         call igop(GLLNID,GLLEL,'+  ',NELGT)

!     mapping of nek5000 to p4est element distribution
        call izero(NP4_LGLNID,el_countl)
        il = 0
        do eg=1,NELGT
            if (GLLNID(eg).eq.NID) then
                il = il+1
                loop : do jl=0,NP-1
                    if (offset(jl+1).ge.eg) then
                        NP4_LGLNID(il) = jl
                        exit loop
                    endif
                enddo loop
            endif
        enddo

      else
         call nekp4est_abort('Error: ParMETIS requires NELGT=NELGV')
      endif

#else
      call nekp4est_abort('Error: nekp4est requires ParMETIS')
#endif

!     Count number of elements on this processor
      NELT=0
      NELV=0
      do eg=1,NELGT
         if (GLLNID(eg).eq.NID) then
            if (eg.le.NELGV) NELV=NELV+1
            if (eg.le.NELGT) NELT=NELT+1
         endif
      enddo

!     update mapping count
      NP4_IMAP = NP4_IMAP + 1

      return
      end
!=======================================================================
