!> @file nekp4est_gs.f
!! @ingroup nekp4est
!! @brief Gater-scatter related routines
!! @author Adam Peplinski
!! @date Jun 28, 2016
!=======================================================================
!> @brief Get topology information.
!! @details Get global vertex, face and edge numberring together with
!!  hanging node and orinetation information for faces and edges.
      subroutine nekp4est_topol_get()
      implicit none

!     local variables
      integer lnelt ! local number of elements
!-----------------------------------------------------------------------
!     global node numberring; hangign node mark
      call nekp4est_node_get(lnelt)

!     face, edge orientation
      call nekp4est_prm_get(lnelt)

      return
      end
!=======================================================================
!> @brief Get node information.
!! @details Get global vertex, face and edge numberring together with
!!  hanging node mark.
!! @param[out] lnelt  local number of elements; p4est count
!! @todo Check strange transpose of vertex numberring caused by data
!!  redistribution
      subroutine nekp4est_node_get(lnelt)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKP4EST'
      include 'NEKP4EST_TOPOL'

!     argument list
      integer lnelt ! local number of elements; p4est count

!     local variables
      integer lnoden ! local number of indep. nodes; p4est count
      integer gnoden ! local number of independent nodes not shared with other proc.
      integer vnode ! number of nodes; counted vertices, faces, edges
      parameter (vnode = 3**LDIM-1)
      integer*8 lnodes(vnode*LELT) ! global numberring of nodes

      integer el, il, jl ! loop index
      integer itmp
      integer ntot

!     node splitting and sorting
      integer*8 itmp8, lnvrt8, lnfcs8, lnedg8
      integer*8 glnodes ! global number of independent nodes
      integer prm(vnode*LELT) ! permutation array
      integer*8 unodes(vnode*LELT) ! unique local nodes
      integer ioff(vnode*LELT) ! iffset position

!     communicator
      integer gs_handle ! gather-scatter handle

      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

!     functions
      integer*8 i8glsum, i8gl_running_sum, mod1

#ifdef DEBUG
      character*2 str
      integer iunit, ierr
#endif
!-----------------------------------------------------------------------
!     reset refinement flag
!      NP4_IFRESET_DSS = .FALSE.

      call nekp4est_log
     $      (NP4_LP_PRD,'Get node numberring.')

!     create gll points numbering
      call fp4est_lnodes_new()

!     get lnodes numbering to nekton
      call fp4est_msh_get_lnode(lnelt, lnoden, gnoden, lnodes,
     $      NP4_HNG_ELM, NP4_HNG_FSC, NP4_HNG_EDG)

!     destroy p4est lnode
      call fp4est_lnodes_del()

!     test mesh size with respect to p4est
      if (lnelt.ne.NP4_NELT) call nekp4est_abort
     $     ('Error: nekp4est_glbnr_get; lnelt /= np4_nelt')

!     data transfer
!     this data transfer causes strange transpose in vertex numberring
!     for now I cant find the reason
      call nekp4est_node_transfer(lnodes,vnode,lnelt)

!     p4est provides continuous node numberring without distiguisng between
!     different node type. Split and cout nodes in groups: vertices, faces
!     and edges.

!     get global number of unique nodes
      itmp8 = gnoden
      glnodes = i8glsum(itmp8,1)

!     sort nodes
      ntot = vnode*lnelt
      call i8sort(lnodes,prm,ntot)
!     get unique nodes and offset
      itmp = 1 ! unique node counter
      unodes(1) = lnodes(1)
      ioff(1) = 1
      do il=2,ntot
         if (unodes(itmp).ne.lnodes(il)) then
            itmp = itmp +1
            unodes(itmp) = lnodes(il)
            ioff(itmp) = il
         endif
      enddo
      ioff(itmp+1) = ntot + 1

!     update local node number
      lnoden = itmp

!     set gather-scatter communicator
      call gs_setup(gs_handle,unodes,lnoden,nekcomm,mp)

!     mark nodes with process id
      itmp8 = NID
      do il=1,lnoden
         unodes(il) = itmp8
      enddo

!     find min for given node
      call gs_op(gs_handle,unodes,3,3,0)

!     count local and mark non-local nodes
      NP4_GLBNR_LFCS = 0
      NP4_GLBNR_LEDG = 0
      NP4_GLBNR_LVRT = 0
      do il=1,lnoden
         if (unodes(il).eq.itmp8) then
            itmp = mod1(prm(ioff(il)),vnode)
            if (itmp.le.NP4_NFCS) then ! face
               NP4_GLBNR_LFCS = NP4_GLBNR_LFCS  + 1
               unodes(il) = NP4_GLBNR_LFCS
            else if (itmp.le.(NP4_NFCS + NP4_NEDG)) then ! edge
               NP4_GLBNR_LEDG = NP4_GLBNR_LEDG  + 1
               unodes(il) = NP4_GLBNR_LEDG
            else if (itmp.le.(NP4_NFCS + NP4_NEDG + NP4_NVRT)) then ! vertex
               NP4_GLBNR_LVRT = NP4_GLBNR_LVRT  + 1
               unodes(il) = NP4_GLBNR_LVRT
            else

            endif
         else
            unodes(il) = 0
         endif
      enddo

!     get global counts
      NP4_GLBNR_GFCS = i8glsum(NP4_GLBNR_LFCS,1)
      NP4_GLBNR_GEDG = i8glsum(NP4_GLBNR_LEDG,1)
      NP4_GLBNR_GVRT = i8glsum(NP4_GLBNR_LVRT,1)

!     correctenss test
      if (glnodes.ne.(NP4_GLBNR_GFCS+NP4_GLBNR_GEDG+NP4_GLBNR_GVRT))
     $ call nekp4est_abort
     $             ('ERROR: nekp4est_topol_get; inconsistent glnodes')

!     get number of nodes on processes with lower id
      lnfcs8 = i8gl_running_sum(NP4_GLBNR_LFCS) - NP4_GLBNR_LFCS
      lnedg8 = i8gl_running_sum(NP4_GLBNR_LEDG) - NP4_GLBNR_LEDG
      lnvrt8 = i8gl_running_sum(NP4_GLBNR_LVRT) - NP4_GLBNR_LVRT

!     renumber nodes
      do il=1,lnoden
         if (unodes(il).gt.0) then
            itmp = mod1(prm(ioff(il)),vnode)
            if (itmp.le.NP4_NFCS) then ! face
               unodes(il) = unodes(il) + lnfcs8
            else if (itmp.le.(NP4_NFCS + NP4_NEDG)) then ! edge
            write(*,*) 'tutaj',nid
               unodes(il) = unodes(il) + lnedg8
            else if (itmp.le.(NP4_NFCS + NP4_NEDG + NP4_NVRT)) then ! vertex
               unodes(il) = unodes(il) + lnvrt8
            else

            endif
         else
            unodes(il) = 0
         endif
      enddo

!     redistribute node numbers
      call gs_op(gs_handle,unodes,3,1,0)

!     put node numberring back
      do il=1,lnoden
         do jl = ioff(il),ioff(il + 1) - 1
            lnodes(jl) = unodes(il)
         enddo
      enddo

!     reverse permute lnoden
      call i8swap_rip(lnodes,prm,ntot)

!     extract numberring
      itmp = 0
      do el=1,lnelt
!     face
         do il = 1,NP4_NFCS
            itmp = itmp + 1
            NP4_GLBNR_FCS(il,el) = lnodes(itmp)
         enddo
!     edge
         do il = 1,NP4_NEDG
            itmp = itmp + 1
            NP4_GLBNR_EDG(il,el) = lnodes(itmp)
         enddo
!     vertex
         do il = 1,NP4_NVRT
            itmp = itmp + 1
            NP4_GLBNR_VRT(il,el) = lnodes(itmp)
         enddo
      enddo

!     free communicator
      call gs_free (gs_handle)

#ifdef DEBUG
!     for testing
      call io_file_freeid(iunit, ierr)
      write(str,'(i2.2)') NID
      open(unit=iunit,file='topol.txt'//str)

      write(iunit,*) lnelt, glnodes
      write(iunit,*) NP4_GLBNR_GVRT, NP4_GLBNR_GFCS, NP4_GLBNR_GEDG
      write(iunit,*) NP4_GLBNR_LVRT, NP4_GLBNR_LFCS, NP4_GLBNR_LEDG
      write(iunit,*) 'FACE'
      do el=1,lnelt
         write(iunit,*) el,(NP4_GLBNR_FCS(il,el),il=1,NP4_NFCS)
      enddo
#if N_DIM == 3
      write(iunit,*) 'EDGE'
      do el=1,lnelt
         write(iunit,*) el,(NP4_GLBNR_EDG(il,el),il=1,NP4_NEDG)
      enddo
#endif
      write(iunit,*) 'VERTEX'
      do el=1,lnelt
         write(iunit,*) el,(NP4_GLBNR_VRT(il,el),il=1,NP4_NVRT)
      enddo
      close(iunit)
#endif
      return
      end
!=======================================================================
!> @brief Get face, edge orinetation.
!! @param[in] lnelt  local number of elements; p4est count
      subroutine nekp4est_prm_get(lnelt)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'
      include 'NEKP4EST_TOPOL'

!     argument list
      integer lnelt ! local number of elements; p4est count

!     local variables
      integer el, il ! loop index

!#ifdef DEBUG
      character*2 str
      integer iunit, ierr
!#endif
!-----------------------------------------------------------------------
      call nekp4est_log
     $      (NP4_LP_PRD,'Get element elighment.')

!     reset arrays
      il = LELT*NP4_NFCS
      call izero(NP4_PRM_FCS,il)
      il = LELT*NP4_NEDG
      call izero(NP4_PRM_EDG,il)

!     get face, edge orientation


!     redistribute data

!#ifdef DEBUG
!     for testing
      call io_file_freeid(iunit, ierr)
      write(str,'(i2.2)') NID
      open(unit=iunit,file='aligment.txt'//str)
      write(iunit,*) lnelt
      write(iunit,*) 'FACE'
      do el=1,lnelt
         write(iunit,*) el,(NP4_PRM_FCS(il,el),il=1,NP4_NFCS)
      enddo
#if N_DIM == 3
      write(iunit,*) 'EDGE'
      do el=1,lnelt
         write(iunit,*) el,(NP4_PRM_EDG(il,el),il=1,NP4_NEDG)
      enddo
#endif
      close(iunit)
!#endif
      return
      end
!=======================================================================
!> @brief Generate global GLL points numbering to generate communicator
!! @details I use symmetric vertex, face and edge numberring decribed in
!!  setedge. It is consistent as well with p4est ordering. The numberring
!! is consistent with memory aligment and does not take into account
!! face and edge orientation. I simply assume all elements are aligned.
!! @param[out]   glo_num    global node numberring
!! @param[out]   ngv        number of unique nodes
!! @param[in]    nx         number of points in element along single dimension
!! @param[in]    nel        element number
!! @param[in]    ifcenter   do we include element interior
!! @todo Add face and edge orientation
      subroutine nekp4est_setvert(glo_num,ngv,nx,nel,ifcenter)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKP4EST'
      include 'NEKP4EST_TOPOL'

!     argument list
      integer nx, nel
      logical ifcenter
      integer*8 glo_num(nx**LDIM*nel),ngv

!     local variables
      integer nn, mm ! counter
      integer el, il, jl, kl ! loop index
      integer nxl, els, nstrd(NP4_NMAX), nstrdr(NP4_NMAX) ! data stride
      integer nstrt(NP4_NMAX) ! data start

!     functions
      integer iglsum

!     for testing
#ifdef DEBUG
      character*2 str
      integer iunit, ierr
#endif
!-----------------------------------------------------------------------
!     exclude vertices
      nxl = nx - 2
!     reset global numberring
      nn = nx**LDIM*nel
      call i8zero(glo_num,nn)

!     check if we have to work with vertices only
      if (nx.eq.2) then
!     vertices only
         call i8copy(glo_num,NP4_GLBNR_VRT,nn)
         ngv = NP4_GLBNR_GVRT
      elseif (nx.gt.2) then
!     vertices
         nn = 0
!     data stride
         nstrd(1) = 1
         nstrd(2) = nx-1
         nstrd(3) = nx*(nx-2) + 1
         nstrd(4) = nx -1
#if N_DIM == 3
         nstrd(5) = nx*nx*(nx-2) + 1
         nstrd(6) = nx-1
         nstrd(7) = nx*(nx-2) + 1
         nstrd(8) = nx -1
#endif
         do el=1,nel
            do il=1,NP4_NVRT
!     vertex posiotin in the element
               nn = nn + nstrd(il)
               glo_num(nn) = NP4_GLBNR_VRT(il,el)
            enddo
         enddo
         ngv = NP4_GLBNR_GVRT
!     faces
#if N_DIM == 2
         els = nx*nx
!     data stride
         nstrd(1) = nx
         nstrd(2) = nx
         nstrd(3) = 1
         nstrd(4) = 1
!     data start
         nstrt(1) = 1
         nstrt(2) = nx
         nstrt(3) = 1
         nstrt(4) = nx*(nx-1) + 1
         do el=1,nel
!     element start
            mm = (el-1) * els
            do il=1,NP4_NFCS
!     initial face posiotin in the element
               nn = mm +nstrt(il)
               do jl=1,nxl
                  nn = nn + nstrd(il)
                  glo_num(nn) = ngv + nxl*(NP4_GLBNR_FCS(il,el)-1) + jl
               enddo
            enddo
         enddo
         ngv = ngv + NP4_GLBNR_GFCS*nxl
!     element centre
         if (ifcenter) then
            do el=1,nel
!     element start + initial centre posiotin in the element
               nn = (el-1) * els + nx + 1

                  do jl=1,nxl
                     do il=1,nxl
                        nn = nn + 1
                        glo_num(nn) = ngv + ((LGLEL(el)-1)*nxl +
     $                                             jl - 1)*nxl + il
                     enddo
                     nn = nn + 2
                  enddo
            enddo
!     get global sum of counted vertices
!     I do it because I have no information which mesh (V;T) is used
            mm =  iglsum(nel,1)
            ngv = ngv + mm*nxl*nxl
         endif
#else
         els = nx*nx*nx
!     data stride; within 1D row
         nstrd(1) = nx
         nstrd(2) = nx
         nstrd(3) = 1
         nstrd(4) = 1
         nstrd(5) = 1
         nstrd(6) = 1
!     data stride; between 1D rows
         nstrdr(1) = 2*nx
         nstrdr(2) = 2*nx
         nstrdr(3) = nx*(nx-1) + 2
         nstrdr(4) = nx*(nx-1) + 2
         nstrdr(5) = 2
         nstrdr(6) = 2
!     data start
         nstrt(1) = nx*nx + 1
         nstrt(2) = nx*(nx + 1)
         nstrt(3) = nx*nx + 1
         nstrt(4) = nx*(2*nx-1) + 1
         nstrt(5) = nx + 1
         nstrt(6) = nx*nx*(nx-1) + nx + 1
         do el=1,nel
!     element start
            mm = (el-1) * els
            do il=1,NP4_NFCS
!     initial face posiotin in the element
               nn = mm +nstrt(il)
               do kl=1,nxl
                  do jl=1,nxl
                     nn = nn + nstrd(il)
                     glo_num(nn) = ngv + (nxl*(NP4_GLBNR_FCS(il,el)-1) +
     $                             kl-1)*nxl + jl
                  enddo
                  nn = nn + nstrdr(il)
               enddo
            enddo
         enddo
         ngv = ngv + NP4_GLBNR_GFCS*nxl*nxl
!     edge
!     data stride
         nstrd(1)  = 1
         nstrd(2)  = 1
         nstrd(3)  = 1
         nstrd(4)  = 1
         nstrd(5)  = nx
         nstrd(6)  = nx
         nstrd(7)  = nx
         nstrd(8)  = nx
         nstrd(9)  = nx*nx
         nstrd(10) = nx*nx
         nstrd(11) = nx*nx
         nstrd(12) = nx*nx
!     data start
         nstrt(1)  = 1
         nstrt(2)  = nx*(nx-1) + 1
         nstrt(3)  = nx*nx*(nx-1) + 1
         nstrt(4)  = nx*(nx*nx-1) + 1
         nstrt(5)  = 1
         nstrt(6)  = nx
         nstrt(7)  = nx*nx*(nx-1) + 1
         nstrt(8)  = nx*nx*(nx-1) + nx
         nstrt(9)  = 1
         nstrt(10) = nx
         nstrt(11) = nx*(nx-1) + 1
         nstrt(12) = nx*nx
         do el=1,nel
!     element start
            mm = (el-1) * els
            do il=1,NP4_NEDG
!     initial edge posiotin in the element
               nn = mm +nstrt(il)
               do jl=1,nxl
                  nn = nn + nstrd(il)
                  glo_num(nn) = ngv + nxl*(NP4_GLBNR_EDG(il,el)-1) + jl
               enddo
            enddo
         enddo
         ngv = ngv + NP4_GLBNR_GEDG*nxl
!     element centre
         if (ifcenter) then
            do el=1,nel
!     element start + initial centre posiotin in the element
               nn = (el-1) * els + nx*(nx+1) + 1
               do kl=1,nxl
                  do jl=1,nxl
                     do il=1,nxl
                        nn = nn + 1
                        glo_num(nn) = ngv + (((LGLEL(el)-1)*nxl +
     $                                  kl-1)*nxl + jl - 1)*nxl + il
                     enddo
                     nn = nn + 2
                  enddo
                  nn = nn + 2*(nx+1)
               enddo
            enddo
!     get global sum of conted vertices
!     I do it because I have no information which mesh (V;T) is used
            mm =  iglsum(nel,1)
            ngv = ngv + mm*nxl*nxl*nxl
         endif
#endif
      else ! wrong nx; it cannot be lower thant 2
         call nekp4est_abort
     $   ('Error: nekp4est_setvert; nx must be >= 2')
      endif

#ifdef DEBUG
!     for testing
      call io_file_freeid(iunit, ierr)
      write(str,'(i2.2)') NID
      open(unit=iunit,file='setvert.txt'//str)
      write(iunit,*) nel, nx,  ngv, LDIM
      nn = 0
      do el=1,nel
         write(iunit,*) 'ELEMENT ',el
#if N_DIM == 2
         do jl=1,nx
            write(iunit,*) (glo_num(nn + il),il=1,nx)
            nn = nn + nx
         enddo
#else
         do kl=1,nx
            write(iunit,*) 'K ',kl
            do jl=1,nx
               write(iunit,*) (glo_num(nn + il),il=1,nx)
               nn = nn + nx
            enddo
         enddo
#endif
      enddo
      close(iunit)
#endif
      return
      end

!=======================================================================
!     get global vertex numbering
      subroutine nekp4est_vertex_get(vertex)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'
      include 'NEKP4EST_TOPOL'

!     argument list
      integer vertex (NP4_NVRT*LELT)

!     locall variables
      integer el, il, itmp
!-----------------------------------------------------------------------
      call nekp4est_log(NP4_LP_PRD,'Get vertex numbering')

      itmp = 0
      do el = 1, NELT
         do il = 1, NP4_NVRT
            itmp = itmp + 1
            vertex(itmp) = NP4_GLBNR_VRT(il,el)
         enddo
      enddo

      return
      end
!=======================================================================
!     to define global communicator for gll points based
!     on p4est node numbering
      subroutine nekp4est_setupds
     $     (gs_handle,nx,ny,nz,nel,melg,glo_num)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
!      include 'INPUT'
!      include 'PARALLEL'
!      include 'NONCON'

!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     input
      integer gs_handle,nx,ny,nz,nel,melg
      integer*8 glo_num(nx*ny*nz*nel)
#if 0
      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ntot
      integer*8 ngv

!     simple timing
      real t0, t1
!     functions
      real dnekclock
      character*2 str
!-----------------------------------------------------------------------
      t0 = dnekclock()

      write(str,'(i2)') nx
      call nekp4est_log(6,'nekp4est_setupds nx='//str//CHAR(0))

!     Global-to-local mapping for gs
!     p4est gives grid points numbering, so I just copy those vaues
!     to glo_num and count number of unique lnodes
      call nekp4est_setvert(glo_num,ngv,nx,nel,.FALSE.)

!     Adam Peplinski;
!     I'm not sure this subroutine should be called for p4est,
!     but for now I keep it here
      if(nid.eq.0) write(6,*) 'call usrsetvert'
      call usrsetvert(glo_num,nel,nx,ny,nz)
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrsetvert'

!     Initialize gather-scatter code
      ntot      = nx*ny*nz*nel

!     to awoid testing different communication methods
!      call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)
      call nekp4est_gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

!     call gs_chkr(glo_num)

      t1 = dnekclock() - t0
      if (nid.eq.0) then
         write(6,1) t1,gs_handle,nx,ngv,melg
    1    format('   setupds time',1pe11.4,' seconds ',2i3,2i12)
      endif


!     set interpolation operator for hanging fases
      call set_jmat(JMAT,gs_handle,nx,nz)
#endif
      return
      end
!=======================================================================
