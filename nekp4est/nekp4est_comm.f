!> @file nekp4est_comm.f
!! @ingroup nekp4est
!! @brief Set of global communication routines to redistribute elements.
!! @author Adam Peplinski
!! @date Jun 05, 2016
!=======================================================================
!> @brief Redistribute initial tree data between processors
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine nekp4est_tree_transfer
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'NEKP4EST'

!     local variables
!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     integer and real array sizes
      integer isw,rsw
!     integer and real arrays depend on number of passive scalars
!     I hope this number doesnt make the transfer arrays bigger than
!     the assumed size of scratch arrays
!     /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
!     /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
      parameter (isw=(1+3*6*(LDIMT1+1)+12+2+1))
      parameter (rsw=5*6*(LDIMT1+1))

!     transfer arrays
      real vr(rsw,LELT)
      integer vi(isw,LELT)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

      integer lnelt,itmp
      integer eg,il,jl,kl       ! loop index
      integer key(1)            ! required by crystal rauter; for sorting

!     temporary storage for face renumbering
      character*3 cbt(6)
      character*1 curvt(6)
      real        bt(5,6)

#ifdef DEBUG
!     for testing; take into accout that debugging can be run on less than 100 processes
      integer iunit, ierr
      character*2 str
#endif
!-----------------------------------------------------------------------
      call nekp4est_log(NP4_LP_PRD,'Redistributing forest data.')

!     take number of local p4est elements
      lnelt = NP4_NELT

!     single send
!     pack real array
      do eg=1,lnelt
!     B.C.
         itmp = 1
         do il=ibc,nfldt
            do jl=1,n_fcs
               do kl=1,5
                  itmp = itmp+1
                  vr(itmp,eg) = BC(kl,jl,eg,il)
               enddo
            enddo
         enddo
      enddo

!     pack integer array
      do eg=1,lnelt
!     global element number
         itmp = np4_nelit + eg
         vi(1,eg) = itmp
!     processor id
         vi(2,eg) = GLLNID(itmp)

!     B.C.
         itmp = 2
         do il=ibc,nfldt
            do jl=1,n_fcs
               do kl=1,3
                  itmp = itmp+1
                  vi(itmp,eg) = ichar(CBC(jl,eg,il)(kl:kl))
               enddo
            enddo
         enddo

!     curvature
         do il=1,12
            itmp = itmp+1
            vi(itmp,eg) = ichar(CCURVE(il,eg))
         enddo

!     element group
         itmp = itmp+1
         vi(itmp,eg) = IGROUP(eg)

!     element level
         itmp = itmp+1
         vi(itmp,eg) = NP4_LEVEL(eg)
      enddo

!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,rsw,2)

!     test local element number
      if (lnelt.ne.NELT)
     $     call nekp4est_abort('Error: tree_transfer; lnelt /= nelt')

!     sort elements acording to their global number
      key(1) = 1
      call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,rsw,key,1)

!     put data back
!     real arrays
      do eg=1,NELT
!     B.C.
         itmp = 1
         do il=ibc,nfldt
            do jl=1,n_fcs
               do kl=1,5
                  itmp = itmp+1
                  BC(kl,jl,eg,il) = vr(itmp,eg)
               enddo
            enddo
         enddo
      enddo

!     integer arrays
      do eg=1,NELT

#ifdef DEBUG
!     check of back communication arrays; just for testing
!     global element number
         if (vi(1,eg).ne.LGLEL(eg)) call nekp4est_abort
     $     ('Error: tree_transfer; back comm. err.; el')
!     origin processor id
         if (vi(2,eg).ne.NP4_LGLNID(eg)) call nekp4est_abort
     $     ('Error: tree_transfer; back comm. err.; nid')
#endif

!     B.C.
         itmp = 2
         do il=ibc,nfldt
            do jl=1,n_fcs
               do kl=1,3
                  itmp = itmp+1
                  CBC(jl,eg,il)(kl:kl) = char(vi(itmp,eg))
               enddo
            enddo
         enddo

!     curvature
         do il=1,12
            itmp = itmp+1
            CCURVE(il,eg) = char(vi(itmp,eg))
         enddo

!     element group
         itmp = itmp+1
         IGROUP(eg) = vi(itmp,eg)

!     element level
         itmp = itmp+1
         NP4_LEVEL(eg) = vi(itmp,eg)
      enddo

!     change ordering of faces
!     BC, CBC
      do kl=ibc,nfldt
         do eg=1,NELT           !  SWAP TO PREPROCESSOR NOTATION
            call chcopy(cbt,CBC(1,eg,kl)  ,n_fcs*3)
            call copy  ( bt, BC(1,1,eg,kl),n_fcs*5)
            do il=1,n_fcs
               call copy  ( BC(1,il,eg,kl), bt(1,EFACE1(il)),5)
               call chcopy(CBC(  il,eg,kl),cbt(  EFACE1(il)),3)
            enddo
!     renumber faces within 'E  ' and 'P  ' bc's
            do il=1,n_fcs
               if (CBC(il,eg,kl).eq.'P  '
     $            .or.CBC(il,eg,kl).eq.'E  ') then
                  jl = BC(2,il,eg,kl)
                  BC(2,il,eg,kl) = EFACE(jl)
               endif
            enddo
         enddo
      enddo

!     CCURVE
      do eg=1,NELT           !  SWAP TO PREPROCESSOR NOTATION
            call chcopy(curvt,CCURVE(1,eg),n_fcs)
            do il=1,n_fcs
               call chcopy(CCURVE(il,eg),curvt(EFACE1(il)),1)
            enddo
      enddo

#ifdef DEBUG
!     testing
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
      open(unit=iunit,file='tree_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv
      do il=1,lnelt
         write(iunit,*) 'element', il, LGLEL(il), vi(1,il)
         write(iunit,*) 'curvature'
         do jl=1,12
            write(iunit,*) jl,CCURVE(jl,il)
         enddo
         write(iunit,*) 'BC'
         do jl=1,n_fcs
            write(iunit,*) jl,CBC(jl,il,1)
            write(iunit,*) jl,(BC(kl,jl,il,1),kl=1,5)
         enddo
      enddo
      close(iunit)
#endif
      return
      end
!=======================================================================
!> @brief Global communication for error estimator mark
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine nekp4est_mark_transfer
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'NEKP4EST'

!     local variables
!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     integer and real array sizes
      integer isw,rsw
      parameter (isw=(1+2))
      parameter (rsw=1)

!     transfer arrays
      real vr(rsw)
      integer vi(isw,LELT)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

      integer lnelt,eg,itmp
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
!     redistribute integer refinement mark back to p4est element distribution
      call nekp4est_log(NP4_LP_PRD,'Redistributing MARK array.')

!     take number of local nek elements
      lnelt = NELT

!     single send
!     pack integer array
      do eg=1,lnelt
!     global element number
         vi(1,eg) = LGLEL(eg)
!     processor id
         vi(2,eg) = NP4_LGLNID(eg)

!     mark array
         vi(3,eg) = NP4_MARK(eg)
      enddo

!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

!     test local element number
      if (lnelt.ne.NP4_NELT) call nekp4est_abort
     $     ('Error: mark_transfer; lnelt /= nelt')

!     sort elements acording to their global number
      key(1) = 1
      call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

!     put data back
!     integer arrays
      do eg=1,NP4_NELT

#ifdef DEBUG
!     check of back communication arrays; just for testing
!     global element number
        itmp = NP4_NELIT + eg
        if (vi(1,eg).ne.itmp) call nekp4est_abort
     $     ('Error: mark_transfer; back comm. err.; el')
!     origin processor id
        if (vi(2,eg).ne.GLLNID(itmp)) call nekp4est_abort
     $     ('Error: mark_transfer; back comm. err.; nid')
#endif

!     mark array
         NP4_MARK(eg) = vi(3,eg)
      enddo

      return
      end
!=======================================================================
!> @brief Global communication for refinement history; MAP part
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine nekp4est_hst_map_transfer
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'NEKP4EST'

!     local variables
!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     integer and real array sizes
      integer isw,rsw
      parameter (isw=(1+2))
      parameter (rsw=1)

!     transfer arrays
      real vr(rsw)
      integer vi(isw,LELT)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

      integer lnelt,itmp
      integer eg,jl             ! loop index
      integer key(1)            ! required by crystal rauter; for sorting

#ifdef DEBUG
!     for testing
      integer iunit,ierr
      character*2 str
!     call number
      integer icalled, icalled_w
      save icalled
      data icalled /0/
      parameter (icalled_w=6)
#endif
!-----------------------------------------------------------------------
      call nekp4est_log(NP4_LP_PRD,'Redistributing history MAP.')

!     take number of local p4est elements
      lnelt = NP4_NELT

!     count array entries
      jl = 0

!     single send
!     pack integer array
      do eg=1,lnelt
!     if the element was unchanged
        itmp = NP4_GLGL_MAP(eg)
        if (itmp.ne.0) then
!     count array entries
            jl=jl+1

!     old global element number
            vi(1,jl) = itmp
!     old processor id
            vi(2,jl) = NP4_GLLNID_O(itmp)

!     current global element number
            vi(3,jl) = eg + NP4_NELIT
        endif
      enddo

!     test local element number
      if (jl.ne.NP4_MAP_NR) call nekp4est_abort
     $     ('Error: map_transfer;wrong jl')
!     correct number of array entries
      lnelt = NP4_MAP_NR

!     reset position array to 0
      call izero(NP4_GLGL_MAP,NP4_NELT_O)
      call ifill(NP4_GLGL_MAP_NID,-1,NP4_NELT_O)

!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

      if(lnelt.gt.0) then
!     test local element number
        if (lnelt.gt.NP4_NELT_O) call nekp4est_abort
     $         ('Error: map_transfer;wrong lnelt')

!     sort elements acording to their global number
        key(1) = 1
        call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

!     put data back
!     integer arrays
        do eg=1,lnelt
!     map array
            itmp = NP4_GLLEL_O(vi(1,eg))
            NP4_GLGL_MAP(itmp) = vi(3,eg)
            NP4_GLGL_MAP_NID(itmp) = GLLNID(vi(3,eg))
        enddo
      endif

!     save current number of "valid" elements
      NP4_MAP_NR = lnelt

#ifdef DEBUG
!     testing
!     count call number
      icalled = icalled+1
      if (icalled.eq.icalled_w) then
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
      open(unit=iunit,file='MAP_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv
      write(iunit,*) 'Global info_old', NP4_NELGT_O, NP4_NELT_O,
     $                NP4_NELV_O
      do eg=1,NP4_NELT_O
         write(iunit,*) eg, NP4_GLGL_MAP(eg),NP4_GLGL_MAP_NID(eg)
      enddo
      close(iunit)
      endif
#endif
      return
      end
!=======================================================================
!> @brief Global communication for refinement history; RFN part
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine nekp4est_hst_rfn_transfer
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'NEKP4EST'

!     local variables
!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     integer and real array sizes
      integer isw,rsw
      parameter (isw=(2+2))
      parameter (rsw=1)

!     transfer arrays
      real vr(rsw)
      integer vi(isw,LELT)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

      integer lnelt
      integer eg,itmp,itmp2,itmp3,il,jl
      integer key(1)            ! required by crystal rauter; for sorting

#ifdef DEBUG
!     for testing
      integer iunit,ierr
      character*2 str
!     call number
      integer icalled, icalled_w
      save icalled
      data icalled /0/
      parameter (icalled_w=6)
#endif
!-----------------------------------------------------------------------
      call nekp4est_log(NP4_LP_PRD,'Redistributing history RFN.')

!     take number of local refined p4est elements
      lnelt = NP4_RFN_NR

!     single send
!     pack integer array
      do eg=1,lnelt
!     old global parent element number
        itmp = NP4_GLGL_RFN(2,eg)
!     old global parent element number
        vi(1,eg) = itmp
!     old parent element processor id
        vi(2,eg) = NP4_GLLNID_O(itmp)

!     current global element number
        vi(3,eg) = NP4_GLGL_RFN(1,eg)
!     child position
        vi(4,eg) = NP4_GLGL_RFN(3,eg)

      enddo

!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)


      if (lnelt.gt.0) then

!     test local element number; it must be multiplication of n_vrts
        if (mod(lnelt,n_vrts).ne.0) call nekp4est_abort
     $         ('Error: rfn_transfer;wrong lnelt')

!     sort elements acording to their global number
        key(1) = 1
        call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

!     put data back
!     count new allocated elements
        il = NP4_NELT_O

!     integer arrays
        do eg=1,lnelt,n_vrts
            itmp2 = NP4_GLLEL_O(vi(1,eg))
            itmp3 = il
            do jl=0,n_vrts-1
!     which child
                itmp = vi(4,eg+jl)
!     new global element number
                NP4_GLGL_RFN(1,eg+itmp)  = vi(3,eg+jl)
!     local parent position
            if (itmp2.ne.NP4_GLLEL_O(vi(1,eg+jl))) call nekp4est_abort
     $         ('Error: rfn_transfer; wrong parent number')
                NP4_GLGL_RFN(2,eg+itmp)  = itmp2
!     local child position
                if (itmp.eq.0) then
                    NP4_GLGL_RFN(3,eg+itmp)  = itmp2
                else
                    il=il+1
                    NP4_GLGL_RFN(3,eg+itmp)  = itmp3+itmp
                endif
            enddo
        enddo

!     check arrays size
        itmp=il
        if (itmp.gt.LELT) call nekp4est_abort
     $    ('Error: rfn_transfer; too many refined blocks')

!     fill map array
!     first reset values from NP4_NELT_O+1 to itmp
        do eg=NP4_NELT_O+1,itmp
            NP4_GLGL_MAP(eg) = 0
            NP4_GLGL_MAP_NID(eg) = -1
        enddo
!     fill checking if every is correct
        do eg=1,lnelt
            il = NP4_GLGL_RFN(3,eg)
            if (NP4_GLGL_MAP(il).eq.0) then
                NP4_GLGL_MAP(il) = NP4_GLGL_RFN(1,eg)
                NP4_GLGL_MAP_NID(il) = GLLNID(NP4_GLGL_RFN(1,eg))
            else
                call nekp4est_abort
     $         ('Error: rfn_transfer;index allready used')
            endif
        enddo
      else
        itmp=NP4_NELT_O
      endif
!     save current number of "valid" elements
      NP4_RFN_NR = lnelt
      NP4_RFN_NR_S = itmp

#ifdef DEBUG
!     testing
!     count call number
      icalled = icalled+1
      if (icalled.eq.icalled_w) then
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
!     RFN
      open(unit=iunit,file='RFN_dat.txt'//str)
      write(iunit,*) 'Global info', nid, NP4_RFN_NR
      do eg=1,NP4_RFN_NR
         write(iunit,*) eg, (NP4_GLGL_RFN(il,eg),il=1,3)
      enddo
      close(iunit)
!     modified MAP
      open(unit=iunit,file='MAP2_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv,
     $                NP4_RFN_NR_S
      write(iunit,*) 'Global info_old', NP4_NELGT_O, NP4_NELT_O,
     $                NP4_NELV_O
      do eg=1,NP4_RFN_NR_S
         write(iunit,*) eg, NP4_GLGL_MAP(eg), NP4_GLGL_MAP_NID(eg)
      enddo
      close(iunit)
      endif
#endif
      return
      end
!=======================================================================
!> @brief Global communication for refinement history; CRS part
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine nekp4est_hst_crs_transfer
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'NEKP4EST'

!     local variables
!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     integer and real array sizes
      integer isw,rsw
      parameter (isw=(2+2))
      parameter (rsw=1)

!     transfer arrays
      real vr(rsw)
      integer vi(isw,LELT)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

!     tmp array
      integer vi_tmp(LELT)
      common /screv/ vi_tmp

      integer lnelt

      integer eg,itmp,itmp2,itmp3,il,jl, kl
      integer key(1)            ! required by crystal rauter; for sorting

!     functions
      integer igl_running_sum

#ifdef DEBUG
!     for testing
      integer iunit,ierr
      character*2 str
!     call number
      integer icalled, icalled_w
      save icalled
      data icalled /0/
      parameter (icalled_w=6)
#endif
!-----------------------------------------------------------------------
      call nekp4est_log(NP4_LP_PRD,'Redistributing history CRS.')

!     first get tmp global numbering of all elements to be removed
!     count number of elements on lower number processors
      itmp  = igl_running_sum(NP4_CRS_NR)

      if (NP4_CRS_NR.gt.0) then
!     calculate tmp global element number
!       starting point
        itmp = NELGT+(itmp - NP4_CRS_NR)*(n_vrts-1)
        do il = 1, NP4_CRS_NR
            do jl=2,n_vrts
                itmp = itmp + 1
                NP4_GLGL_CRS(1,jl,il) = itmp
            enddo
        enddo

!     fill transfer array,
!     this one will be send to processors owning the elements for coarsening
!     count array entries
        kl=0
        do il = 1, NP4_CRS_NR
            itmp =GLLNID(NP4_GLGL_CRS(1,1,il))
            do jl=1,n_vrts
                kl=kl+1
!     old global child element number
                vi(1,kl) = NP4_GLGL_CRS(2,jl,il)
!     old child element processor id
                vi(2,kl) = NP4_GLLNID_O(vi(1,kl))
!     new element global number
                vi(3,kl) = NP4_GLGL_CRS(1,jl,il)
!     new parent element processor id
                vi(4,kl) = itmp
            enddo
        enddo
        lnelt = kl
      else
        lnelt = 0
      endif

!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

      if (lnelt.gt.0) then
!     sort elements acording to their global number
        key(1) = 1
        call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

!     fill MAP arrays
        do il=1,lnelt
            itmp = NP4_GLLEL_O(vi(1,il))
            if (NP4_GLGL_MAP(itmp).eq.0) then
                NP4_GLGL_MAP(itmp) = vi(3,il)
                if (vi(3,il).le.NELGT) then
!     element mover, but not destroyed
                    NP4_GLGL_MAP_NID(itmp) = GLLNID(vi(3,il))
                else
!     element destroyed
                    NP4_GLGL_MAP_NID(itmp) = vi(4,il)
                endif
            else
                call nekp4est_abort
     $         ('Error: crs_transfer;index allready used')
            endif
        enddo
      endif

!     this one will be send to processors owning the elements resulting from coarsening
      if (NP4_CRS_NR.gt.0) then
!     fill transfer array
!     count array entries
        kl=0
        do il = 1, NP4_CRS_NR
            itmp = NP4_GLGL_CRS(1,1,il)
            itmp2 = GLLNID(itmp)
            do jl=1,n_vrts
                kl=kl+1
!     new parent element global number
                vi(1,kl) = itmp
!     new parent element processor id
                vi(2,kl) = itmp2
!     new child element global number
                vi(3,kl) = NP4_GLGL_CRS(1,jl,il)
!     child position
                vi(4,kl) = jl

            enddo
        enddo
        lnelt = kl
      else
        lnelt = 0
      endif

!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

      if (lnelt.gt.0) then

!     test local element number; it must be multiplication of n_vrts
        if (mod(lnelt,n_vrts).ne.0) call nekp4est_abort
     $         ('Error: crs_transfer;wrong lnelt')

!     sort elements acording to their global number
        key(1) = 1
        call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

!     put data back
!     count new allocated elements
        il = NELT

!     integer arrays
        do eg=1,lnelt/n_vrts
            kl=(eg-1)*n_vrts
            itmp = GLLEL(vi(1,kl+1))
            itmp3 = il
            do jl=1,n_vrts
!     local parent position
                if (itmp.ne.GLLEL(vi(1,kl+jl)))
     $          call nekp4est_abort
     $          ('Error: crs_transfer; wrong parent number')

!     which child
                itmp2 = vi(4,kl+jl)
!     new global element number
                NP4_GLGL_CRS(1,itmp2,eg)  = vi(3,kl+jl)

!     local child position
                if (itmp2.eq.1) then
                    NP4_GLGL_CRS(2,itmp2,eg)  = itmp
                else
                    il=il+1
                    NP4_GLGL_CRS(2,itmp2,eg)  = itmp3 + itmp2 - 1
                endif
            enddo
        enddo

!     check arrays size
        itmp=il
        if (itmp.gt.LELT) call nekp4est_abort
     $    ('Error: crs_transfer; too many refined blocks')

        lnelt = lnelt/n_vrts

      endif
!     save current number of "valid" elements
      NP4_CRS_NR = lnelt

!     check if there are still some zeros in MAP
      do eg=1,NP4_RFN_NR_S
        if (NP4_GLGL_MAP(eg).eq.0) call nekp4est_abort
     $    ('Error: crs_transfer; inconsistent MAP')
      enddo

#ifdef DEBUG
!     testing
!     count call number
      icalled = icalled+1
      if (icalled.eq.icalled_w) then
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
!     CRS
      open(unit=iunit,file='CRS_dat.txt'//str)
      write(iunit,*) 'Global info', nid, NP4_CRS_NR
      do eg=1,NP4_CRS_NR
        do il=1,n_vrts
            write(iunit,*) eg,il,(NP4_GLGL_CRS(jl,il,eg),jl=1,2)
        enddo
      enddo
      close(iunit)
!     modified MAP
      open(unit=iunit,file='MAP3_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv,
     $                NP4_RFN_NR_S
      write(iunit,*) 'Global info_old', NP4_NELGT_O, NP4_NELT_O,
     $                NP4_NELV_O
      do eg=1,NP4_RFN_NR_S
         write(iunit,*) eg, NP4_GLGL_MAP(eg),NP4_GLGL_MAP_NID(eg)
      enddo
      close(iunit)
      endif
#endif
      return
      end
!=======================================================================
!> @brief Redistribute single variable after refinement
!! @param[in] vr      redistributed vector
!! @param[in] rsw     array size
!! @param[in] lbuff   array size
!! @remarks This routine uses global scratch space SCRUZ
      subroutine nekp4est_vec_transfer(vr,rsw, lbuff)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'NEKP4EST'

!     argument list
      integer rsw, lbuff
      real vr(rsw,lbuff)

!     face and vertex number
      integer n_fcs, n_vrts
      parameter (n_fcs=2*LDIM, n_vrts=2**LDIM)

!     integer array size
      integer isw
      parameter (isw=2)

!     transfer arrays
      integer vi(isw,LELT)
      integer*8 vl              ! required by crystal rauter
      common /scruz/ vi

      integer lnelt,eg
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
!     take number of local p4est elements
      lnelt = NP4_RFN_NR_S

!     single send
!     pack integer array
      do eg=1,lnelt
!     global element number
         vi(1,eg) = NP4_GLGL_MAP(eg)
!     processor id
         vi(2,eg) = NP4_GLGL_MAP_NID(eg)
      enddo

!     min aray size
      eg = min(LELT,lbuff)
!     transfer arrays
      call crystal_tuple_transfer
     $     (cr_h,lnelt,eg,vi,isw,vl,0,vr,rsw,2)

!     test local element number
      if (lnelt.ne.(NELT+NP4_CRS_NR*(n_vrts-1))) call nekp4est_abort
     $     ('Error: v1_transfer; lnelt /= nelt'//CHAR(0))

!     sort elements acording to their global number
      key(1) = 1
      call crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,rsw,key,1)

      return
      end
!=======================================================================

!-----------------------------------------------------------------------

!=======================================================================
