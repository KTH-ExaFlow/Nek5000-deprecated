!> @file nekp4est_tool.f
!! @ingroup nekp4est
!! @brief Set of tools for nekp4est
!! @author Adam Peplinski
!! @date Feb 26, 2016
!=======================================================================
!> @brief Write log messages
!! @param[in] priority  log priority
!! @param[in] logs      log body
      subroutine nekp4est_log(priority,logs)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     argument list
      integer priority
      character(len=*) logs

!     local variables
      character(len=NP4_LSTL_LOG) llogs
      integer slen
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./
!-----------------------------------------------------------------------
      if (.not.ifcalled) then
        ifcalled=.TRUE.

!     register nekp4est package
        call fsc_pkg_reg(NP4_PKG_ID,priority,'nekp4est'//CHAR(0),
     $  'SEM solver.'//CHAR(0))
      endif

!     check log length
      slen = len_trim(logs)
      if (slen.ge.NP4_LSTL_LOG) then
         llogs = 'Warning; too long log string, shortening'//CHAR(0)
         call fsc_log (NP4_PKG_ID, 1, priority,trim(llogs))
         llogs = logs(1:NP4_LSTL_LOG-1)//CHAR(0)
      else
         llogs = trim(logs)//CHAR(0)
      endif

      call fsc_log (NP4_PKG_ID, 1, priority, trim(llogs))

      return
      end
!=======================================================================
!> @brief SC based abort function
!! @param[in] logs      log body
      subroutine nekp4est_abort(logs)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     argument list
      character(len=*) logs

!     local variables
      character(len=NP4_LSTL_LOG) llogs
      integer slen
!-----------------------------------------------------------------------
!     check log length
      slen = len_trim(logs)
      if (slen.ge.NP4_LSTL_LOG) then
         llogs = 'Warning; too long log string, shortening'//CHAR(0)
         call fsc_log (NP4_PKG_ID, 1, NP4_LP_ERR,trim(llogs))
         llogs = logs(1:NP4_LSTL_LOG-1)//CHAR(0)
      else
         llogs = trim(logs)//CHAR(0)
      endif

      call fsc_abort(trim(llogs))
      return
      end
!=======================================================================
!> @brief SC based check abort function
!! @param[in] ierr      error indicator
!! @param[in] logs      log body
      subroutine nekp4est_chk_abort(ierr,logs)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     argument list
      integer ierr
      character(len=*) logs

!     local variables
      character(len=NP4_LSTL_LOG) llogs
      integer slen
!-----------------------------------------------------------------------
!     check log length
      slen = len_trim(logs)
      if (slen.ge.NP4_LSTL_LOG) then
         llogs = 'Warning; too long log string, shortening'//CHAR(0)
         call fsc_log (NP4_PKG_ID, 1, NP4_LP_ERR,trim(llogs))
         llogs = logs(1:NP4_LSTL_LOG-1)//CHAR(0)
      else
         llogs = trim(logs)//CHAR(0)
      endif

      call fsc_check_abort(ierr,trim(llogs))
      return
      end
!=======================================================================
!> @brief integer8 version of isort
!! @details Use Heap Sort (p 231 Num. Rec., 1st Ed.)
!! @param[inout] as   sorted array
!! @param[out]   ind  permutation
!! @param[in]    nl   array size
      subroutine i8sort(as,ind,nl)
      implicit none

!     argument list
      integer nl
      integer*8 as(nl)
      integer ind(nl)

!     local variables
      integer im, ie, ii, il, jl
      integer*8 aa
!-----------------------------------------------------------------------
!     initial permuatation
      do im=1,nl
         ind(im)=im
      enddo

      if (nl.eq.1) return
!     initial middle and end position
      im=nl/2+1
      ie=nl
!     infinite loop
      do
         if (im.gt.1) then
            im = im-1
            aa = as(im)
            ii = ind(im)
         else
            aa =  as(ie)
            ii = ind(ie)
            as(ie) = as(1)
            ind(ie)= ind(1)
            ie=ie-1
            if (ie.eq.1) then
                as(1) = aa
               ind(1) = ii
               return
            endif
         endif
         il=im
         jl=im+im
         do while (jl.le.ie)
            if (jl.lt.ie) then
               if (as(jl).lt.as(jl+1) ) jl=jl+1
            endif
            if (aa.lt.as(jl)) then
                as(il) = as(jl)
               ind(il) = ind(jl)
               il=jl
               jl=jl+jl
            else
               jl=ie+1
            endif
         enddo
         as(il) = aa
         ind(il) = ii
      enddo

      end
!=======================================================================
!> @brief modiffied version of iswap_ip
!! @details In-place reverse permutation
!! @param[inout] as   permuted array
!! @param[in]    pi   permutation
!! @param[in]    nl   array size
      subroutine i8swap_rip(as,pi,nl)
      implicit none

!     argument list
      integer nl
      integer*8 as(nl)
      integer pi(nl)

!     local variables
      integer il, jl
      integer*8 as1, as2
      integer j, k, loop_start, last, next
      character*100 str
!-----------------------------------------------------------------------
      do il=1,nl
         if (pi(il).gt.0) then   ! not swapped
            as1 = as(il)
            loop_start = il
            next = pi(il)
            loop : do jl=il,nl
               if (next.lt.0) then
                  write (str,*) 'Hey! i8swap_ip problem.',jl,il,nl,next
                  call fsc_abort(str)
               elseif (next.eq.loop_start) then
                  as(next) = as1
                  pi(next) = -pi(next)
                  exit loop
               else
                  as2 = as(next)
                  as(next) = as1
                  as1 = as2
                  last = next
                  next = pi(last)
                  pi(last) = -pi(last)
               endif
            enddo loop
         endif
      enddo

!     reset permutation mark
      do il=1,nl
         pi(il) = -pi(il)
      enddo
      return
      end
!=======================================================================
!> @brief Important topology variables
!! @details Originally set in setup_topo, but I use them for data transfer,
!!  so to be sure they are initialised I do it here.
      subroutine nekp4est_setup_topo
      implicit none
!-----------------------------------------------------------------------
!     prenek to nek face reordering
      call initds

      return
      end
!=======================================================================
