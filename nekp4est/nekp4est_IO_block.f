!> @file nekp4est_IO_block.f
!! @brief block data to initialise common block for I/O routines for nekp4est
!! @details Following Nek5000 standard I keep block data in seaprate file
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
      block data IO_common_init
!     keeep track of max iunit generated
      integer IO_iunit_min, IO_iunit_max
      common /IO_iunit/ IO_iunit_min, IO_iunit_max
      data IO_iunit_min/200/
      data IO_iunit_max/200/
      end
