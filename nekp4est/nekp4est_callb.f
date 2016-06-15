!> @file nekp4est_callb.f
!! @ingroup nekp4est
!! @brief Callback subroutines required by p4est wrapper
!! @author Adam Peplinski
!! @date Apr 1, 2016
!=======================================================================
!> @brief Initialise mesh data in p4est mem. structures; e.g. curvature
!! and flows boundary conditions.
!! @details This subroutine is not used in this code version
!! @param[in]  tree    tree root number
!! @param[out] imsh    V- mesh marker
!! @param[out] igrp    element group
!! @param[out] cbcl    bc type
!! @param[out] bcl     bc parameters
!! @param[out] crvl    bc curvature
      subroutine nek_init_msh_dat(tree,imsh,igrp,cbcl,bcl,crvl)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'

!     max element and vertex number
      integer lelm,lpts, lfcs
      integer n_fcs, n_vrts,n_fvrts
      parameter (lelm=LELT)
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM,n_fvrts=n_vrts/2)
      parameter (lpts=n_vrts*lelm)
      parameter (lfcs=n_fcs*lelm)
#if N_DIM == 3
      integer n_edg, ledgs
      parameter (n_edg=12,ledgs=n_edg*lelm)
#endif
!     for number of boundary conditions
      integer nlsize
      parameter (nlsize=N_PSCL+1)

!     argument list
      integer imsh, tree, igrp
      character*3  cbcl (6,0:nlsize)
      real*8       bcl  (5,6,0:nlsize)
      integer       crvl(6)

!     local variables
      integer i, j, k, treel

!-----------------------------------------------------------------------

      return
      end
!=======================================================================
!> @brief Get mesh data form p4est mem. structures; e.g. curvature and
!! flows boundary conditions.
!! @details Required by p4est_msh_get_dat
!! @param[in] lnel    local element number
!! @param[in] gnel    global element number
!! @param[in] imsh    V- mesh marker
!! @param[in] igrp    element group
!! @param[in] level   refinement level
!! @param[in] cbcl    bc type
!! @param[in] bcl     bc parameters
!! @param[in] crvl    bc curvature
      subroutine nek_get_msh_dat(lnel,gnel,imsh,igrp,level,
     $           cbcl,bcl,crvl)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'NEKP4EST'

!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     argument list
      integer lnel, gnel, imsh, igrp, level
      character*3  cbcl (6,0:LDIMT1)
      real*8       bcl  (5,6,0:LDIMT1)
      integer      crvl(6)

!     local variables
      integer i, j, k
!-----------------------------------------------------------------------
!     count local number of elements and check wether tree is sorted
      if (lnel.ne.EL_COUNT) 
     $     call nekp4est_abort('ERROR: unsorted tree.')
      EL_COUNT = EL_COUNT +1

!     check correctness of global element numberring
      if ((gnel-lnel).ne.NP4_NELIT)
     $     call nekp4est_abort('ERROR: global numberring.')

!     is it V or T mesh
      if (imsh.eq.0.and.gnel.gt.NELGV)
     $     call nekp4est_abort('ERROR: V, T mesh missmatch.')

!     take element data
!     element group
      IGROUP(EL_COUNT) = igrp

!     element level
      NP4_LEVEL(EL_COUNT) = level

!     boundary conditions
      do j=IBC,NFLDT
        do i=1,n_fcs
            CBC(i,EL_COUNT,j) = cbcl(i,j)
        enddo
      enddo
      do k=IBC,NFLDT
        do j=1,n_fcs
            do i=1,5
                BC(i,j,EL_COUNT,k) = bcl(i,j,k)
            enddo
        enddo
      enddo

!     curvature data
!     I mark deformed as 'D' without going into details
!     it should be fine as setrzer and setdef do simple check
!     if(CCURVE(,).ne.' ')
      do i=1,6
        if (crvl(i).ne.0) CCURVE(i,EL_COUNT) = 'D'
      enddo

      return
      end
!=======================================================================
!> @brief Get refinement history data form p4est
!! @details Required by p4est_msh_get_dat. It is used to perform mesh
!! refinement on Nek5000 side following p4est data.
!! @param[in] gnel       old global element number (unchanged element) or child position (refined element)
!! @param[in] gnel_prnt  old global parent number (refined element)
!! @param[in] gnel_chld  list of removed children (coarsened element)
      subroutine nek_get_msh_hst(gnel,gnel_prnt,gnel_chld)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'

!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     argument list
      integer gnel, gnel_prnt, gnel_chld(n_vrts)

!     local variables
      integer il

!     functions
      integer ivlmax
!-----------------------------------------------------------------------
!     was there any refinement/coarsening performed
      il = ivlmax(gnel_chld,n_vrts)
      if (gnel_prnt.eq.-1.and.il.eq.-1) then
!     no refinement/coarsening
!     count elements
        NP4_MAP_NR = NP4_MAP_NR +1
!     fill arays
        NP4_GLGL_MAP(EL_COUNT) = gnel+1
      else if(gnel_prnt.ne.-1) then
!     refinement
!     count elements
        NP4_RFN_NR = NP4_RFN_NR +1
!     fill arays
        NP4_GLGL_MAP(EL_COUNT) = 0
        NP4_GLGL_RFN(1,NP4_RFN_NR) = EL_COUNT + NP4_NELIT
        NP4_GLGL_RFN(2,NP4_RFN_NR) = gnel_prnt+1
        NP4_GLGL_RFN(3,NP4_RFN_NR) = gnel
      else
!     coarsening
!     count elements
        NP4_CRS_NR = NP4_CRS_NR +1
!     fill arays
        NP4_GLGL_MAP(EL_COUNT) = 0
!     new global position
        NP4_GLGL_CRS(1,1,NP4_CRS_NR) = EL_COUNT + NP4_NELIT
!     old global position
        NP4_GLGL_CRS(2,1,NP4_CRS_NR) = gnel_chld(1)+1
        do il = 2,n_vrts
!     new global position; unknown
            NP4_GLGL_CRS(1,il,NP4_CRS_NR) = 0
!     old global position
            NP4_GLGL_CRS(2,il,NP4_CRS_NR) = gnel_chld(il)+1
        enddo
      endif

      return
      end
!=======================================================================
!> @brief Transfer refinement mark to p4est
!! @details Required by
!! @param[in]  gnel  global element number
!! @param[out] mark  refinement mark; 0 - no change, 1 - refinement, -1 - coarsering
      subroutine nek_refine_mark(gnel,mark)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'NEKP4EST'
      include 'NEKP4EST_REFINE'

!     argument list
      integer gnel, mark

!     local variables
      integer itmp
!-----------------------------------------------------------------------
!     move mark array to p4est
!     global to local element number
      itmp = gnel - NP4_NELIT+1
!     is the local element number correct
      if (itmp.ge.1.and.itmp.le.NP4_NELT) then
        mark = NP4_MARK(itmp)
!     wrong element number
      else
        call nekp4est_abort
     $     ('Error: nek_refine_mark; wrong element num.')
      endif

      return
      end
!=======================================================================
