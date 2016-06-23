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
      character(len=NP4_LSTL_LOG) logs

!     local variables
      integer pkg_id
      save pkg_id
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./
!-----------------------------------------------------------------------
      if (.not.ifcalled) then
        ifcalled=.TRUE.

!     register nekp4est package
        call fsc_pkg_reg(pkg_id,priority,'nekp4est'//CHAR(0),
     $  'SEM solver.'//CHAR(0))
      endif

      call fsc_log (pkg_id, 1, priority, trim(logs)//CHAR(0))

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
      character(len=NP4_LSTL_LOG) logs
!-----------------------------------------------------------------------
      call fsc_abort(trim(logs)//CHAR(0))
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
      character(len=NP4_LSTL_LOG) logs
!-----------------------------------------------------------------------
      call fsc_check_abort(ierr,trim(logs)//CHAR(0))
      return
      end
!=======================================================================
