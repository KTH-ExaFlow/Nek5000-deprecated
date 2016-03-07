!> @brief set of tools for nekp4est
!! @author Adam Peplinski
!! @date Feb 26, 2016
!=======================================================================
!> @brief Write log messages
!! @parameter[in] priority  log priority
!! @parameter[in] logs      log body
      subroutine nekp4est_log(priority,logs)
      implicit none
!     input variables
      integer priority
      character*200 logs

!     local variables
      integer icalld, pkg_id
      save icalld, pkg_id
      data icalld /0/
!-----------------------------------------------------------------------
      if (icalld.eq.0) then
        icalld=1

!     register nekp4est package
        call fsc_pkg_reg(pkg_id,priority,'nekp4est'//CHAR(0),
     $  'SEM solver.'//CHAR(0))
      endif

      call fsc_log (pkg_id, 1, priority, trim(logs)//CHAR(0))

      return
      end
!=======================================================================
!> @brief SC based abort function
!! @parameter[in] logs      log body
      subroutine nekp4est_abort(logs)
      implicit none

!     input variables
      character*200 logs
!-----------------------------------------------------------------------
      call fsc_abort(trim(logs)//CHAR(0))
      return
      end
!=======================================================================
!> @brief SC based check abort function
!! @parameter[in] ierr      error indicator
!! @parameter[in] logs      log body
      subroutine nekp4est_chk_abort(ierr,logs)
      implicit none

!     input variables
      integer ierr
      character*200 logs
!-----------------------------------------------------------------------
      call fsc_check_abort(ierr,trim(logs)//CHAR(0))
      return
      end
!=======================================================================
