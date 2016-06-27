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
