!> @file nekp4est_refine.f
!! @ingroup nekp4est
!! @brief Set of routines to perform interpolation for refinement/coarsening.
!! @author Adam Peplinski
!! @date Jun 07, 2016
!=======================================================================
!> @brief Get interpolation parameters for octal interpolation
      subroutine nekp4est_genwz
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'WZ_DEF'
      include 'WZ'              ! ZGM(123)
      include 'NEKP4EST_REFINE'

!     local variables
      integer il, jl, kl     !     loop indexes
      integer nt2

!     tmp  array to calculate interpolation positions
      integer nmax
      parameter (nmax=84) ! like in speclib
      real tmpl(nmax)
!-----------------------------------------------------------------------
!     zero arrays
!     Volume interpolation
      nt2 = LX1*LX1*2
      call rzero(IXAMR1CF,nt2)
      call rzero(IXTAMR1CF,nt2)
      call rzero(IXAMR1FC,nt2)
      call rzero(IXTAMR1FC,nt2)
      nt2 = LX2*LX2*2
      call rzero(IXAMR2CF,nt2)
      call rzero(IXTAMR2CF,nt2)
      call rzero(IXAMR2FC,nt2)
      call rzero(IXTAMR2FC,nt2)
      nt2 = LY1*LY1*2
      call rzero(IYAMR1CF,nt2)
      call rzero(IYTAMR1CF,nt2)
      call rzero(IYAMR1FC,nt2)
      call rzero(IYTAMR1FC,nt2)
      nt2 = LY2*LY2*2
      call rzero(IYAMR2CF,nt2)
      call rzero(IYTAMR2CF,nt2)
      call rzero(IYAMR2FC,nt2)
      call rzero(IYTAMR2FC,nt2)
      if (IF3D) then
         nt2 = LZ1*LZ1*2
         call rzero(IZAMR1CF,nt2)
         call rzero(IZTAMR1CF,nt2)
         call rzero(IZAMR1FC,nt2)
         call rzero(IZTAMR1FC,nt2)
         nt2 = LZ2*LZ2*2
         call rzero(IZAMR2CF,nt2)
         call rzero(IZTAMR2CF,nt2)
         call rzero(IZAMR2FC,nt2)
         call rzero(IZTAMR2FC,nt2)
      elseif(IFAXIS) then       ! axisymmetric case
         nt2 = LY1*LY1*2
         call rzero(IAAMR1CF,nt2)
         call rzero(IATAMR1CF,nt2)
         call rzero(IAAMR1FC,nt2)
         call rzero(IATAMR1FC,nt2)
         nt2 = LY2*LY2*2
         call rzero(IAAMR2CF,nt2)
         call rzero(IATAMR2CF,nt2)
         call rzero(IAAMR2FC,nt2)
         call rzero(IATAMR2FC,nt2)
      endif
!     multiplicity arrays
      nt2 = LX1*LY1*LZ1
      call rone(IMAMR1,nt2)
      nt2 = LX2*LY2*LZ2
      call rone(IMAMR2,nt2)

!     get interpolation operators
!     Volume interpolation
!     MESH M1
!     coarse -> fine
!     X
!     negative
      do il=1,NX1
         tmpl(il) = 0.5*(ZGM1(il,1) -1.0)
      enddo
      call IGLLM (IXAMR1CF,IXTAMR1CF,ZGM1(1,1),tmpl,NX1,NX1,NX1,NX1)
!     positive; we use symmetry
      do jl=1,NX1
         do il=1,NX1
            IXAMR1CF(NX1-il+1,NX1-jl+1,2)  = IXAMR1CF(il,jl,1)
            IXTAMR1CF(NX1-il+1,NX1-jl+1,2) = IXTAMR1CF(il,jl,1)
         enddo
      enddo

!     Y
!     negative
      do il=1,NY1
         tmpl(il) = 0.5*(ZGM1(il,2) -1.0)
      enddo
      call IGLLM (IYAMR1CF,IYTAMR1CF,ZGM1(1,2),tmpl,NY1,NY1,NY1,NY1)
!     positive; we use symmetry
      do jl=1,NY1
         do il=1,NY1
            IYAMR1CF(NY1-il+1,NY1-jl+1,2)  = IYAMR1CF(il,jl,1)
            IYTAMR1CF(NY1-il+1,NY1-jl+1,2) = IYTAMR1CF(il,jl,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         do il=1,NZ1
            tmpl(il) = 0.5*(ZGM1(il,3) -1.0)
         enddo
         call IGLLM (IZAMR1CF,IZTAMR1CF,ZGM1(1,3),tmpl,NZ1,NZ1,NZ1,NZ1)
!     positive; we use symmetry
         do jl=1,NZ1
            do il=1,NZ1
               IZAMR1CF(NZ1-il+1,NZ1-jl+1,2)  = IZAMR1CF(il,jl,1)
               IZTAMR1CF(NZ1-il+1,NZ1-jl+1,2) = IZTAMR1CF(il,jl,1)
            enddo
         enddo
      else
         IZAMR1CF(NZ1,NZ1,1) = 1.0
         IZTAMR1CF(NZ1,NZ1,1) = 1.0
         IZAMR1CF(NZ1,NZ1,2) = 1.0
         IZTAMR1CF(NZ1,NZ1,2) = 1.0
      endif

!     fine -> coarse
!     X
!     negative
      nt2 = NX1/2 + mod(NX1,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM1(il,1) + 1.0
      enddo
      call IGLLM (IXAMR1FC,IXTAMR1FC,ZGM1(1,1),tmpl,NX1,nt2,NX1,NX1)
!     positive; we use symmetry
      do jl=1,NX1
         do il=1,nt2
            IXAMR1FC(NX1-il+1,NX1-jl+1,2)  = IXAMR1FC(il,jl,1)
            IXTAMR1FC(NX1-jl+1,NX1-il+1,2) = IXTAMR1FC(jl,il,1)
         enddo
      enddo

!     Y
!     negative
      nt2 = NY1/2 + mod(NY1,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM1(il,2) + 1.0
      enddo
      call IGLLM (IYAMR1FC,IYTAMR1FC,ZGM1(1,2),tmpl,NY1,nt2,NY1,NY1)
!     positive; we use symmetry
      do jl=1,NY1
         do il=1,nt2
            IYAMR1FC(NY1-il+1,NY1-jl+1,2)  = IYAMR1FC(il,jl,1)
            IYTAMR1FC(NY1-jl+1,NY1-il+1,2) = IYTAMR1FC(jl,il,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         nt2 = NZ1/2 + mod(NZ1,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZGM1(il,3) + 1.0
         enddo
         call IGLLM (IZAMR1FC,IZTAMR1FC,ZGM1(1,3),tmpl,NZ1,nt2,NZ1,NZ1)
!     positive; we use symmetry
         do jl=1,NZ1
            do il=1,nt2
               IZAMR1FC(NZ1-il+1,NZ1-jl+1,2)  = IZAMR1FC(il,jl,1)
               IZTAMR1FC(NZ1-jl+1,NZ1-il+1,2) = IZTAMR1FC(jl,il,1)
            enddo
         enddo
      else
         nt2 = NZ1/2 + mod(NZ1,2)
         IZAMR1FC(nt2,NZ1,1) = 1.0
         IZTAMR1FC(NZ1,nt2,1) = 1.0
         IZAMR1FC(nt2,NZ1,2) = 1.0
         IZTAMR1FC(NZ1,nt2,2) = 1.0
      endif

!     MESH M2
!     coarse -> fine
!     X
!     negative
      do il=1,NX2
         tmpl(il) = 0.5*(ZGM2(il,1) -1.0)
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IXAMR2CF,IXTAMR2CF,ZGM2(1,1),tmpl,NX2,NX2,NX2,NX2)
      else                      ! P-N-P_N-2
         call IGLM (IXAMR2CF,IXTAMR2CF,ZGM2(1,1),tmpl,NX2,NX2,NX2,NX2)
      endif

!     positive; we use symmetry
      do jl=1,NX2
         do il=1,NX2
            IXAMR2CF(NX2-il+1,NX2-jl+1,2)  = IXAMR2CF(il,jl,1)
            IXTAMR2CF(NX2-il+1,NX2-jl+1,2) = IXTAMR2CF(il,jl,1)
         enddo
      enddo

!     Y
!     negative
      do il=1,NY2
         tmpl(il) = 0.5*(ZGM2(il,2) -1.0)
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IYAMR2CF,IYTAMR2CF,ZGM2(1,2),tmpl,NY2,NY2,NY2,NY2)
      else                      ! P-N-P_N-2
         call IGLM (IYAMR2CF,IYTAMR2CF,ZGM2(1,2),tmpl,NY2,NY2,NY2,NY2)
      endif
!     positive; we use symmetry
      do jl=1,NY2
         do il=1,NY2
            IYAMR2CF(NY2-il+1,NY2-jl+1,2)  = IYAMR2CF(il,jl,1)
            IYTAMR2CF(NY2-il+1,NY2-jl+1,2) = IYTAMR2CF(il,jl,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         do il=1,NZ2
            tmpl(il) = 0.5*(ZGM2(il,3) -1.0)
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLLM (IZAMR2CF,IZTAMR2CF,ZGM2(1,3),tmpl,
     $           NZ2,NZ2,NZ2,NZ2)
         else                   ! P-N-P_N-2
            call IGLM (IZAMR2CF,IZTAMR2CF,ZGM2(1,3),tmpl,
     $           NZ2,NZ2,NZ2,NZ2)
         endif
!     positive; we use symmetry
         do jl=1,NZ2
            do il=1,NZ2
               IZAMR2CF(NZ2-il+1,NZ2-jl+1,2)  = IZAMR2CF(il,jl,1)
               IZTAMR2CF(NZ2-il+1,NZ2-jl+1,2) = IZTAMR2CF(il,jl,1)
            enddo
         enddo
      else
         IZAMR2CF(NZ2,NZ2,1) = 1.0
         IZTAMR2CF(NZ2,NZ2,1) = 1.0
         IZAMR2CF(NZ2,NZ2,2) = 1.0
         IZTAMR2CF(NZ2,NZ2,2) = 1.0
      endif

!     fine -> coarse
!     X
!     negative
      nt2 = NX2/2 + mod(NX2,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM2(il,1) + 1.0
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IXAMR2FC,IXTAMR2FC,ZGM2(1,1),tmpl,NX2,nt2,NX2,NX2)
      else                      ! P-N-P_N-2
         call IGLM (IXAMR2FC,IXTAMR2FC,ZGM2(1,1),tmpl,NX2,nt2,NX2,NX2)
      endif
!     positive; we use symmetry
      do jl=1,NX2
         do il=1,nt2
            IXAMR2FC(NX2-il+1,NX2-jl+1,2)  = IXAMR2FC(il,jl,1)
            IXTAMR2FC(NX2-jl+1,NX2-il+1,2) = IXTAMR2FC(jl,il,1)
         enddo
      enddo

!     Y
!     negative
      nt2 = NY2/2 + mod(NY2,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM2(il,2) + 1.0
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IYAMR2FC,IYTAMR2FC,ZGM2(1,2),tmpl,NY2,nt2,NY2,NY2)
      else                      ! P-N-P_N-2
         call IGLM (IYAMR2FC,IYTAMR2FC,ZGM2(1,2),tmpl,NY2,nt2,NY2,NY2)
      endif
!     positive; we use symmetry
      do jl=1,NY2
         do il=1,nt2
            IYAMR2FC(NY2-il+1,NY2-jl+1,2)  = IYAMR2FC(il,jl,1)
            IYTAMR2FC(NY2-jl+1,NY2-il+1,2) = IYTAMR2FC(jl,il,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         nt2 = NZ2/2 + mod(NZ2,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZGM2(il,3) + 1.0
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLLM (IZAMR2FC,IZTAMR2FC,ZGM2(1,3),tmpl,
     $           NZ2,nt2,NZ2,NZ2)
         else                   ! P-N-P_N-2
            call IGLM (IZAMR2FC,IZTAMR2FC,ZGM2(1,3),tmpl,
     $           NZ2,nt2,NZ2,NZ2)
         endif
!     positive; we use symmetry
         do jl=1,NZ2
            do il=1,nt2
               IZAMR2FC(NZ2-il+1,NZ2-jl+1,2)  = IZAMR2FC(il,jl,1)
               IZTAMR2FC(NZ2-jl+1,NZ2-il+1,2) = IZTAMR2FC(jl,il,1)
            enddo
         enddo
      else
         nt2 = NZ2/2 + mod(NZ2,2)
         IZAMR2FC(nt2,NZ2,1) = 1.0
         IZTAMR2FC(NZ2,nt2,1) = 1.0
         IZAMR2FC(nt2,NZ2,2) = 1.0
         IZTAMR2FC(NZ2,nt2,2) = 1.0
      endif


!     special treatment of axisymmetric 2D case
      if (.not.IF3D.and.IFAXIS) then
!     copy current interpolation operators
!     mesh 1
         call copy(ICAMR1CF,IYAMR1CF,NY1*NY1*2)
         call copy(ICAMR1FC,IYAMR1FC,NY1*NY1*2)
         call copy(ICTAMR1CF,IYTAMR1CF,NY1*NY1*2)
         call copy(ICTAMR1FC,IYTAMR1FC,NY1*NY1*2)
!     mesh 2
         call copy(ICAMR2CF,IYAMR2CF,NY2*NY2*2)
         call copy(ICAMR2FC,IYAMR2FC,NY2*NY2*2)
         call copy(ICTAMR2CF,IYTAMR2CF,NY2*NY2*2)
         call copy(ICTAMR2FC,IYTAMR2FC,NY2*NY2*2)

!     get intrpolation operators
!     mesh 1
!     coarse -> fine
!     Y
!     negative
         do il=1,NY1
            tmpl(il) = 0.5*(ZAM1(il) -1.0)
         enddo
         call IGLJM (IAAMR1CF,IATAMR1CF,ZAM1,tmpl,NY1,NY1,NY1,NY1)
!     positive; we use symmetry
         do jl=1,NY1
            do il=1,NY1
               IAAMR1CF(NY1-il+1,NY1-jl+1,2)  = IAAMR1CF(il,jl,1)
               IATAMR1CF(NY1-il+1,NY1-jl+1,2) = IATAMR1CF(il,jl,1)
            enddo
         enddo

!     fine -> coarse
!     Y
!     negative
         nt2 = NY1/2 + mod(NY1,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZAM1(il) + 1.0
         enddo
         call IGLJM (IAAMR1FC,IATAMR1FC,ZAM1,tmpl,NY1,nt2,NY1,NY1)
!     positive; we use symmetry
         do jl=1,NY1
            do il=1,nt2
               IAAMR1FC(NY1-il+1,NY1-jl+1,2)  = IAAMR1FC(il,jl,1)
               IATAMR1FC(NY1-jl+1,NY1-il+1,2) = IATAMR1FC(jl,il,1)
            enddo
         enddo

!     MESH M2
!     coarse -> fine
!     Y
!     negative
         do il=1,NY2
            tmpl(il) = 0.5*(ZAM2(il) -1.0)
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLJM (IAAMR2CF,IATAMR2CF,ZAM2,tmpl,
     $           NY2,NY2,NY2,NY2)
         else                   ! P-N-P_N-2
            call IGJM (IAAMR2CF,IATAMR2CF,ZAM2,tmpl,
     $           NY2,NY2,NY2,NY2)
         endif
!     positive; we use symmetry
         do jl=1,NY2
            do il=1,NY2
               IAAMR2CF(NY2-il+1,NY2-jl+1,2)  = IAAMR2CF(il,jl,1)
               IATAMR2CF(NY2-il+1,NY2-jl+1,2) = IATAMR2CF(il,jl,1)
            enddo
         enddo

!     fine -> coarse
!     Y
!     negative
         nt2 = NY2/2 + mod(NY2,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZAM2(il) + 1.0
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLJM (IAAMR2FC,IATAMR2FC,ZAM2,tmpl,
     $           NY2,nt2,NY2,NY2)
         else                   ! P-N-P_N-2
            call IGJM (IAAMR2FC,IATAMR2FC,ZAM2,tmpl,
     $           NY2,nt2,NY2,NY2)
         endif
!     positive; we use symmetry
         do jl=1,NY2
            do il=1,nt2
               IAAMR2FC(NY2-il+1,NY2-jl+1,2)  = IAAMR2FC(il,jl,1)
               IATAMR2FC(NY2-jl+1,NY2-il+1,2) = IATAMR2FC(jl,il,1)
            enddo
         enddo

      endif                     ! axisymmetric

!     for multiplicity
!     mesh 1
!     X
      if (mod(NX1,2).eq.1) then
         il = NX1/2 + 1
         do kl= 1, NZ1
            do jl=1, NY1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            enddo
         enddo
      endif

!     Y
      if (mod(NY1,2).eq.1) then
         jl = NY1/2 + 1
         do kl= 1, NZ1
            do il=1, NX1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            enddo
         enddo
         if (mod(NX1,2).eq.1) then
            il = NX1/2 + 1
            do kl= 1, NZ1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            enddo
         endif
      endif

!     Z
      if (IF3D) then
         if (mod(NZ1,2).eq.1) then
            kl = NZ1/2 + 1
            do jl= 1, NY1
               do il=1, NX1
                  IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
               enddo
            enddo
            if (mod(NX1,2).eq.1) then
               il = NX1/2 + 1
               do jl=1,NY1
                  IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NY1,2).eq.1) then
               jl = NY1/2 + 1
               do il=1,NX1
                  IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NX1,2).eq.1.and.mod(NY1,2).eq.1) then
               il = NX1/2 + 1
               jl = NY1/2 + 1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            endif
         endif
      endif

!     calculate inverse
      nt2 = NX1*NY1*NZ1
      call invcol1(IMAMR1,nt2)

!     mesh 2
!     X
      if (mod(NX2,2).eq.1) then
         il = NX2/2 + 1
         do kl= 1, NZ2
            do jl=1, NY2
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            enddo
         enddo
      endif

!     Y
      if (mod(NY2,2).eq.1) then
         jl = NY2/2 + 1
         do kl= 1, NZ2
            do il=1, NX2
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            enddo
         enddo
         if (mod(NX2,2).eq.1) then
            il = NX2/2 + 1
            do kl= 1, NZ2
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            enddo
         endif
      endif

!     Z
      if (IF3D) then
         if (mod(NZ2,2).eq.1) then
            kl = NZ2/2 + 1
            do jl= 1, NY2
               do il=1, NX2
                  IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
               enddo
            enddo
            if (mod(NX2,2).eq.1) then
               il = NX2/2 + 1
               do jl=1,NY2
                  IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NY2,2).eq.1) then
               jl = NY2/2 + 1
               do il=1,NX2
                  IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NX2,2).eq.1.and.mod(NY2,2).eq.1) then
               il = NX2/2 + 1
               jl = NY2/2 + 1
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            endif
         endif
      endif

!     calculate inverse
      nt2 = NX2*NY2*NZ2
      call invcol1(IMAMR2,nt2)

      return
      end
!=======================================================================
!> @brief Map single variable of coarse element to the fine one
!! @param[out]   vf      fine element vector
!! @param[in]    vc      coarse element vector
!! @param[in]    iel     element number
!! @param[in]    ch_pos  child position
!! @param[in]    limesh   mesh mark (velocity, preasure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,2)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
      subroutine nekp4est_mapcf(vf,vc,iel,ch_pos,limesh,tmp,lnx,lny,lnz)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'GEOM_DEF'
      include 'GEOM'            ! IFRZER
      include 'NEKP4EST'
      include 'NEKP4EST_REFINE'

!     argument list
      integer lnx,lny,lnz
      real vf(lnx,lny,lnz), vc(lnx,lny,lnz)
      integer iel,ch_pos(3) ! ch_pos is child position in 3D
      integer limesh
      real tmp(lnx,lny,lnz,2) ! work array

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer nxy, nyz
      integer iz  ! loop index
!-----------------------------------------------------------------------
      nxy = lnx*lny
      nyz = lny*lnz

!     Use the appropriate derivative- and interpolation operator in
!     the y-direction (= radial direction if axisymmetric).

      if (limesh.eq.1) then ! velocity mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR1CF,IATAMR1CF,lny*lny*2)
            else
               call copy(IYTAMR1CF,ICTAMR1CF,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR1CF(1,1,ch_pos(1)),lnx,vc,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR1CF(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR1CF(1,1,ch_pos(3)),lnz,
     $           vf,lnz)
         else
            call mxm(IXAMR1CF(1,1,ch_pos(1)),lnx,vc,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR1CF(1,1,ch_pos(2)),lny,vf,lny)
         endif
      elseif (limesh.eq.2) then !preasure mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR2CF,IATAMR2CF,lny*lny*2)
            else
               call copy(IYTAMR2CF,ICTAMR2CF,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR2CF(1,1,ch_pos(1)),lnx,vc,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR2CF(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR2CF(1,1,ch_pos(3)),lnz,
     $           vf,lnz)
         else
            call mxm(IXAMR2CF(1,1,ch_pos(1)),lnx,vc,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR2CF(1,1,ch_pos(2)),lny,vf,lny)
         endif
      else  ! there should be as well place for mhd mesh
         logs = 'ERROR: nekp4est_mapcf; wrong limesh'
         call nekp4est_abort(logs)
      endif

      return
      end
!=======================================================================
!> @brief Map single variable of fine element to the coarse one
!! @param[out]   vc      coarse element vector
!! @param[in]    vf      fine element vector
!! @param[in]    iel     element number
!! @param[in]    ch_pos  child position
!! @param[in]    limesh   mesh mark (velocity, preasure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,2)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
      subroutine nekp4est_mapfc(vc,vf,iel,ch_pos,limesh,tmp,lnx,lny,lnz)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'GEOM_DEF'
      include 'GEOM'            ! IFRZER
      include 'NEKP4EST'
      include 'NEKP4EST_REFINE'

!     argument list
      integer lnx,lny,lnz
      real vc(lnx,lny,lnz), vf(lnx,lny,lnz)
      integer iel,ch_pos(3)   !     ch_pos is child position in 3D
      integer limesh
      real tmp(lnx,lny,lnz,2) ! work array

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer nxy, nyz
      integer iz  ! loop index
!-----------------------------------------------------------------------
      nxy = lnx*lny
      nyz = lny*lnz

!     Use the appropriate derivative- and interpolation operator in
!     the y-direction (= radial direction if axisymmetric).

      if (limesh.eq.1) then ! velocity mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR1FC,IATAMR1FC,lny*lny*2)
            else
               call copy(IYTAMR1FC,ICTAMR1FC,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR1FC(1,1,ch_pos(1)),lnx,vf,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR1FC(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR1FC(1,1,ch_pos(3)),lnz,
     $           vc,lnz)
         else
            call mxm(IXAMR1FC(1,1,ch_pos(1)),lnx,vf,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR1FC(1,1,ch_pos(2)),lny,vc,lny)
         endif
      elseif (limesh.eq.2) then !preasure mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR2FC,IATAMR2FC,lny*lny*2)
            else
               call copy(IYTAMR2FC,ICTAMR2FC,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR2FC(1,1,ch_pos(1)),lnx,vf,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR2FC(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR2FC(1,1,ch_pos(3)),lnz,
     $           vc,lnz)
         else
            call mxm(IXAMR2FC(1,1,ch_pos(1)),lnx,vf,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR2FC(1,1,ch_pos(2)),lny,vc,lny)
         endif
      else  ! there should be as well place for mhd mesh
         logs = 'ERROR: nekp4est_mapfc; wrong limesh'
         call nekp4est_abort(logs)
      endif

      return
      end
!=======================================================================
!> @brief Perform single refinement operation of given variable
!! @param[inout] vcf     refined vector
!! @param[in]    el_lst  element list
!! @param[in]    limesh   mesh mark (velocity, preasure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine nekp4est_refine_vs(vcf,el_lst,limesh,tmp,
     $           lnx,lny,lnz,leln)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     argument list
      integer lnx,lny,lnz,leln
      real vcf(lnx,lny,lnz,leln)
      integer el_lst(NP4_NCHILD)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

!     local variables
      integer ch_pos(3)   !     child position in 3D
      integer iel         !     element position
      integer il, nl
!-----------------------------------------------------------------------
!     I assume el_lst(1) gives position of the coarse block
!     and final ch_pos() = 1,1,1
!     copy coarse element
      nl=lnx*lny*lnz
      call copy(tmp,vcf(1,1,1,el_lst(1)),nl)
!     loop over all the children 
      do il= 1,NP4_NCHILD
!     get child position
         iel = el_lst(il)
         ch_pos(3) = (il-1)/4 +1
         ch_pos(2) = mod((il-1)/2,2) +1
         ch_pos(1) = mod(il-1,2) +1
!     refine
         call nekp4est_mapcf(vcf(1,1,1,iel),tmp(1,1,1,1),iel,ch_pos,
     $                       limesh,tmp(1,1,1,2),lnx,lny,lnz)
      enddo
      return
      end
!=======================================================================
!> @brief Perform single coarsening operation of given variable
!! @param[inout] vfc       coarsened vector
!! @param[in]    el_lst    element list
!! @param[in]    limesh   mesh mark (velocity, preasure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine nekp4est_coarse_vs(vfc,el_lst,limesh,tmp,
     $           lnx,lny,lnz,leln)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'
      include 'NEKP4EST_REFINE'

!     argument list
      integer lnx,lny,lnz,leln
      real vfc(lnx,lny,lnz,leln)
      integer el_lst(NP4_NCHILD)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer ch_pos(3)    !     child position in 3D
      integer iel          !     element position
      integer il, nl
!-----------------------------------------------------------------------
!     I assume el_lst(1) gives position of the coarse block
!     and inital ch_pos() = 1,1,1
      nl=lnx*lny*lnz
!     loop over all the children 
      do il= 1,NP4_NCHILD
!     get child position
         iel = el_lst(il)
         ch_pos(3) = (il-1)/4 +1
         ch_pos(2) = mod((il-1)/2,2) +1
         ch_pos(1) = mod(il-1,2) +1
!     coarsen
         call nekp4est_mapfc(tmp(1,1,1,1),vfc(1,1,1,iel),iel,ch_pos,
     $                       limesh,tmp(1,1,1,2),lnx,lny,lnz)
!     sum contributions
         if (il.eq.1) then
            call copy(vfc(1,1,1,iel),tmp,nl)
         else
            call add2(vfc(1,1,1,el_lst(1)),tmp,nl)
         endif
      enddo
!     take into account multiplicity of points
      if (limesh.eq.1) then ! velocity mesh
         call col2(vfc(1,1,1,el_lst(1)),IMAMR1,nl)
      elseif (limesh.eq.2) then !preasure mesh
         call col2(vfc(1,1,1,el_lst(1)),IMAMR2,nl)
      else  ! there should be as well place for mhd mesh
         logs = 'ERROR: nekp4est_coarsen_vs; wrong limesh'
         call nekp4est_abort(logs)
      endif

      return
      end
!=======================================================================
!> @brief Perform all refinement operations for given variable
!! @details Threre are two possible implementations. In following rourtine
!! we performe all block refinemnts for a single variable.
!! @param[inout] vcf     refined vector
!! @param[in]    limesh   mesh mark (velocity, preasure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine nekp4est_refine_vm(vcf,limesh,tmp,lnx,lny,lnz,leln)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     argument list
      integer lnx,lny,lnz,leln
      real vcf(lnx,lny,lnz,leln)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

!     local variables
      integer el_lst(NP4_NCHILD)   ! local element list for refinement
      integer il, jl ! loop index
!-----------------------------------------------------------------------
      do il=0,NP4_RFN_NR-1,NP4_NCHILD
        do jl=1,NP4_NCHILD
            el_lst(jl) = NP4_GLGL_RFN(3,il+jl)
        enddo
        call nekp4est_refine_vs(vcf,el_lst,limesh,tmp,lnx,lny,lnz,leln)
      enddo

      return
      end
!=======================================================================
!> @brief Perform all coarsening operations for given variable
!! @details Threre are two possible implementations. In following rourtine
!! we performe all block coarsenings for a single variable.
!! @param[inout] vfc       coarsened vector
!! @param[in]    limesh   mesh mark (velocity, preasure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine nekp4est_coarse_vm(vfc,limesh,tmp,
     $           lnx,lny,lnz,leln)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'
      include 'NEKP4EST_REFINE'

!     argument list
      integer lnx,lny,lnz,leln
      real vfc(lnx,lny,lnz,leln)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

!     local variables
      integer el_lst(NP4_NCHILD)   ! local element list for refinement
      integer il, jl ! loop index
!-----------------------------------------------------------------------
      do il=1,NP4_CRS_NR
        do jl=1,NP4_NCHILD
            el_lst(jl) = NP4_GLGL_CRS(2,jl,il)
        enddo
        call nekp4est_coarse_vs(vfc,el_lst,limesh,tmp,lnx,lny,lnz,leln)
      enddo

      return
      end
!=======================================================================
!> @brief Perform single refinement operation on set of variables
!! @details Threre are two possible implementations. In following rourtine
!! we performe single block refinement for all variables.
!! @param[in]    el_lst    element list
      subroutine nekp4est_refine_el(el_lst)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NBDINP
      include 'GEOM_DEF'
      include 'GEOM'            ! IFGEOM
      include 'NEKP4EST'

!     argument list
      integer el_lst(NP4_NCHILD)

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl ! loop index
      integer fs, fe
      integer limesh  ! mesh mark (velocity, preasure, mhd)
!     work arrays
      real tmpv(LX1,LY1,LZ1,3), tmpvp(LX2,LY2,LZ2,3),
     $     tmpb(LBX1,LBY1,LBZ1,3), tmpbp(LBX2,LBY2,LBZ2,3)
!-----------------------------------------------------------------------
!     is the firs element a correct one
!     I assume el_lst(1) is an existing parent
!     compare with old global element number
      if (el_lst(1).gt.NP4_NELT_O) then
         logs = 'Error: nekp4est_refine_el el_lst>NP4_NELT_O'
         call nekp4est_abort(logs)
      endif

!     refine variables
!     coordinates; mesh 1
      limesh = 1
      call nekp4est_refine_vs(XM1,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELT)
      call nekp4est_refine_vs(YM1,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELT)
      if(IF3D) call nekp4est_refine_vs(ZM1,el_lst,limesh,tmpv,
     $              LX1,LY1,LZ1,LELT)

!     set every step by setprop in nek_advance
!     interpolate density (VTRANS) and viscosity (VDIFF)
!      fs = 2
!      fe = NFIELD
!      if (IFFLOW) fs = 1
!      if (IFMHD) fe = fe + 1
!      do il=fs, fe
!         call nekp4est_refine_vs(VDIFF (1,1,1,1,il),el_lst,limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!         call nekp4est_refine_vs(VTRANS(1,1,1,1,il),el_lst,limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!      enddo

!     check if the first block is in velocity mesh
!     be carefull; for now the old block number is valid
      if (el_lst(1).le.NP4_NELT_O) then
!     velocity
         call nekp4est_refine_vs(VX,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELV)
         call nekp4est_refine_vs(VY,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELV)
         if(IF3D) call nekp4est_refine_vs(VZ,el_lst,limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)

!     lag arrays velocity
         if (IFFLOW.or.IFMHD) then
            do il =1, NBDINP-1
               call nekp4est_refine_vs(VXLAG(1,1,1,1,il),el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
               call nekp4est_refine_vs(VYLAG(1,1,1,1,il),el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
               if (IF3D) call nekp4est_refine_vs(VZLAG(1,1,1,1,il),
     $              el_lst,limesh,tmpv,LX1,LY1,LZ1,LELV)
            enddo
         endif

!     arrays for time integration
         if (IFTRAN) then
            call nekp4est_refine_vs(ABX1,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            call nekp4est_refine_vs(ABY1,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            if(IF3D) call nekp4est_refine_vs(ABZ1,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)

            call nekp4est_refine_vs(ABX2,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            call nekp4est_refine_vs(ABY2,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            if(IF3D) call nekp4est_refine_vs(ABZ2,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
         endif

!     mesh 2
         limesh = 2
!     preasure
         call nekp4est_refine_vs(PR,el_lst,limesh,
     $              tmpvp,LX2,LY2,LZ2,LELV)


!     lag arrays preasure PRLAG
         if (NBDINP.eq.3) call nekp4est_refine_vs(PRLAG,el_lst,limesh,
     $              tmpvp,LX2,LY2,LZ2,LELV)
      endif                     ! block in velocity mesh

!     mesh 2
      limesh = 2
!     user defined divergence USRDIV
      call nekp4est_refine_vs(USRDIV,el_lst,limesh,
     $              tmpvp,LX2,LY2,LZ2,LELT)

!     temperature and passive scalars
!     mesh 1
      limesh = 1
      if (IFHEAT) then
!     varialbe T; temperature nad passive scalars
         do il=2,NFIELD
!     check T versus V mesh 
            if (IFTMSH(il).or.el_lst(1).le.NELV) then
               call nekp4est_refine_vs(T(1,1,1,1,il-1),el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)

!     lag arrays TLAG
               do jl =1, NBDINP-1
                  call nekp4est_refine_vs(TLAG(1,1,1,1,jl,il-1),el_lst,
     $              limesh,tmpv,LX1,LY1,LZ1,LELT)
               enddo

!     arrays time integration VGRADT[12]
               call nekp4est_refine_vs(VGRADT1(1,1,1,1,il-1),el_lst,
     $              limesh,tmpv,LX1,LY1,LZ1,LELT)
               call nekp4est_refine_vs(VGRADT2(1,1,1,1,il-1),el_lst,
     $              limesh,tmpv,LX1,LY1,LZ1,LELT)
            endif               ! T versus V mesh
         enddo                  ! passive scalar loop
      endif                     ! IFHEAT

!     mhd
!     NOT TESTED AND NOT SURE IT IS OK!!!!!!!!!!!!!!!!
      if (IFMHD) then
!     mesh 3
         limesh = 3
!     mhd blocks reside in velocity mesh
!     be carefull; for now the old block number is valid
         if (el_lst(1).le.NP4_NELT_O) then
!     magnetic field
            call nekp4est_refine_vs(BX,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_refine_vs(BY,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            if(IF3D) call nekp4est_refine_vs(BZ,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)

!     arrays for time integration
            call nekp4est_refine_vs(BBX1,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_refine_vs(BBY1,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            if(IF3D) call nekp4est_refine_vs(BBZ1,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)

            call nekp4est_refine_vs(BBX2,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_refine_vs(BBY2,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            if(IF3D) call nekp4est_refine_vs(BBZ2,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)

!     lag arrays 
            do il =1, NBDINP-1
               call nekp4est_refine_vs(BXLAG(1,il),el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
               call nekp4est_refine_vs(BYLAG(1,il),el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
               if (IF3D) call nekp4est_refine_vs(BZLAG(1,il),el_lst,
     $              limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
            enddo

!     mesh 4
            limesh = 4
            call nekp4est_refine_vs(PM,el_lst,limesh,
     $              tmpbp,LBX2,LBY2,LBZ2,LBELV)

!     lag arrays
            if (NBDINP.eq.3) call nekp4est_refine_vs(PMLAG,el_lst,
     $              limesh,tmpbp,LBX2,LBY2,LBZ2,LBELV)

         endif                  ! block in velocity mesh
      endif

!     moving mesh; BM1LAG
      if (IFGEOM) then
         logs = 'Error: nekp4est_refine_el no moving mesh yet.'
         call nekp4est_abort(logs)
      endif

!     it could be place for perturbation, but interpolation would not
!     give correct base flow structure, so I skip it for now

      return
      end
!=======================================================================
!> @brief Perform single coarsening operation on set of variables
!! @details Threre are two possible implementations. In following rourtine
!! we performe single block coarsening for all variables.
!! @param[in]    el_lst    element list
      subroutine nekp4est_coarse_el(el_lst)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NBDINP
      include 'GEOM_DEF'
      include 'GEOM'            ! IFGEOM
      include 'NEKP4EST'

!     argument list
      integer el_lst(NP4_NCHILD)

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl ! loop index
      integer fs, fe
      integer limesh  ! mesh mark (velocity, preasure, mhd)
!     work arrays
      real tmpv(LX1,LY1,LZ1,3), tmpvp(LX2,LY2,LZ2,3),
     $     tmpb(LBX1,LBY1,LBZ1,3), tmpbp(LBX2,LBY2,LBZ2,3)
!-----------------------------------------------------------------------
!     is the firs element a correct one
!     I assume el_lst(1) is an existing parent
      if (el_lst(1).gt.NELT) then
         logs = 'Error: nekp4est_coarse_el el_lst>NELT'
         call nekp4est_abort(logs)
      endif

!     coarsen variables
!     coordinates; mesh 1
      limesh = 1
      call nekp4est_coarse_vs(XM1,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELT)
      call nekp4est_coarse_vs(YM1,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELT)
      if(IF3D) call nekp4est_coarse_vs(ZM1,el_lst,limesh,tmpv,
     $              LX1,LY1,LZ1,LELT)

!     set every step by setprop in nek_advance
!     interpolate density (VTRANS) and viscosity (VDIFF)
!      fs = 2
!      fe = NFIELD
!      if (IFFLOW) fs = 1
!      if (IFMHD) fe = fe + 1
!      do il=fs, fe
!         call nekp4est_coarse_vs(VDIFF (1,1,1,1,il),el_lst,limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!         call nekp4est_coarse_vs(VTRANS(1,1,1,1,il),el_lst,limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!      enddo

!     check if the first block is in velocity mesh
      if (el_lst(1).le.NELV) then
!     velocity
         call nekp4est_coarse_vs(VX,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELV)
         call nekp4est_coarse_vs(VY,el_lst,limesh,tmpv,LX1,LY1,LZ1,LELV)
         if(IF3D) call nekp4est_coarse_vs(VZ,el_lst,limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)

!     lag arrays velocity
         if (IFFLOW.or.IFMHD) then
            do il =1, NBDINP-1
               call nekp4est_coarse_vs(VXLAG(1,1,1,1,il),el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
               call nekp4est_coarse_vs(VYLAG(1,1,1,1,il),el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
               if (IF3D) call nekp4est_coarse_vs(VZLAG(1,1,1,1,il),
     $              el_lst,limesh,tmpv,LX1,LY1,LZ1,LELV)
            enddo
         endif

!     arrays for time integration
         if (IFTRAN) then
            call nekp4est_coarse_vs(ABX1,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            call nekp4est_coarse_vs(ABY1,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            if(IF3D) call nekp4est_coarse_vs(ABZ1,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)

            call nekp4est_coarse_vs(ABX2,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            call nekp4est_coarse_vs(ABY2,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            if(IF3D) call nekp4est_coarse_vs(ABZ2,el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
         endif

!     mesh 2
         limesh = 2
!     preasure
         call nekp4est_coarse_vs(PR,el_lst,limesh,
     $              tmpvp,LX2,LY2,LZ2,LELV)

!     lag arrays preasure PRLAG
         if (NBDINP.eq.3) call nekp4est_coarse_vs(PRLAG,el_lst,limesh,
     $              tmpvp,LX2,LY2,LZ2,LELV)
      endif                     ! block in velocity mesh

!     mesh 2
      limesh = 2
!     user defined divergence USRDIV
      call nekp4est_coarse_vs(USRDIV,el_lst,limesh,
     $              tmpvp,LX2,LY2,LZ2,LELT)

!     temperature and passive scalars
!     mesh 1
      limesh = 1
      if (IFHEAT) then
!     varialbe T; temperature nad passive scalars
         do il=2,NFIELD
!     check T versus V mesh 
            if (IFTMSH(il).or.el_lst(1).le.NELV) then
               call nekp4est_coarse_vs(T(1,1,1,1,il-1),el_lst,limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)

!     lag arrays TLAG
               do jl =1, NBDINP-1
                  call nekp4est_coarse_vs(TLAG(1,1,1,1,jl,il-1),el_lst,
     $              limesh,tmpv,LX1,LY1,LZ1,LELT)
               enddo

!     arrays time integration VGRADT[12]
               call nekp4est_coarse_vs(VGRADT1(1,1,1,1,il-1),el_lst,
     $              limesh,tmpv,LX1,LY1,LZ1,LELT)
               call nekp4est_coarse_vs(VGRADT2(1,1,1,1,il-1),el_lst,
     $              limesh,tmpv,LX1,LY1,LZ1,LELT)
            endif               ! T versus V mesh
         enddo                  ! passive scalar loop
      endif                     ! IFHEAT

!     mhd
!     NOT TESTED AND NOT SURE IT IS OK!!!!!!!!!!!!!!!!
      if (IFMHD) then
!     mesh 3
         limesh = 3
!     mhd blocks reside in velocity mesh
         if (el_lst(1).le.NELV) then
!     magnetic field
            call nekp4est_coarse_vs(BX,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_coarse_vs(BY,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            if(IF3D) call nekp4est_coarse_vs(BZ,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)

!     arrays for time integration
            call nekp4est_coarse_vs(BBX1,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_coarse_vs(BBY1,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            if(IF3D) call nekp4est_coarse_vs(BBZ1,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)

            call nekp4est_coarse_vs(BBX2,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_coarse_vs(BBY2,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
            if(IF3D) call nekp4est_coarse_vs(BBZ2,el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)

!     lag arrays 
            do il =1, NBDINP-1
               call nekp4est_coarse_vs(BXLAG(1,il),el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
               call nekp4est_coarse_vs(BYLAG(1,il),el_lst,limesh,
     $              tmpb,LBX1,LBY1,LBZ1,LBELV)
               if (IF3D) call nekp4est_coarse_vs(BZLAG(1,il),el_lst,
     $              limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
            enddo

!     mesh 4
            limesh = 4
            call nekp4est_coarse_vs(PM,el_lst,limesh,
     $              tmpbp,LBX2,LBY2,LBZ2,LBELV)

!     lag arrays
            if (NBDINP.eq.3) call nekp4est_coarse_vs(PMLAG,el_lst,
     $              limesh,tmpbp,LBX2,LBY2,LBZ2,LBELV)

         endif                  ! block in velocity mesh
      endif

!     moving mesh; BM1LAG
      if (IFGEOM) then
         logs = 'Error: nekp4est_coarse_el no moving mesh yet.'
         call nekp4est_abort(logs)
      endif

!     it could be place for perturbation, but interpolation would not
!     give correct base flow structure, so I skip it for now

      return
      end
!=======================================================================
!> @brief Local refinement
!! @details Threre are two possible implementations. In following rourtine
!! we performe single block refinement for all variables implemented in
!! nekp4est_refine_el and next turn to the next refinement block. This
!! version gives more possibility for error check.
      subroutine nekp4est_refine_local_el
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer el_lst(NP4_NCHILD)  ! local element list for refinement
      integer il, jl ! loop index
!-----------------------------------------------------------------------
      logs = 'Local data refinement.'
      call nekp4est_log(NP4_LP_PRD,logs)
      do il=0,NP4_RFN_NR-1,NP4_NCHILD
        do jl=1,NP4_NCHILD
            el_lst(jl) = NP4_GLGL_RFN(3,il+jl)
        enddo
        call nekp4est_refine_el(el_lst)
      enddo

      return
      end
!=======================================================================
!> @brief Local corsening
!! @details Threre are two possible implementations. In following rourtine
!! we performe single block coarsening for all variables implemented in
!! nekp4est_coarse_el and next turn to the next coarsening block. This
!! version gives more possibility for error check.
      subroutine nekp4est_coarsen_local_el
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKP4EST'

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer el_lst(NP4_NCHILD)   ! local element list for refinement
      integer il, jl ! loop index
!-----------------------------------------------------------------------
      logs = 'Local data coarsening.'
      call nekp4est_log(NP4_LP_PRD,logs)
      do il=1,NP4_CRS_NR
        do jl=1,NP4_NCHILD
            el_lst(jl) = NP4_GLGL_CRS(2,jl,il)
        enddo
        call nekp4est_coarse_el(el_lst)
      enddo

      return
      end
!=======================================================================
!> @brief Local refinement
!! @details Threre are two possible implementations. In following rourtine
!! we performe all block refinemnts for a single variable implemented in
!! nekp4est_refine_vm and next turn to the variable. This version could
!! be slightly faster.
      subroutine nekp4est_refine_local_v
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NBDINP
      include 'GEOM_DEF'
      include 'GEOM'            ! IFGEOM
      include 'NEKP4EST'

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl ! loop index
      integer fs, fe
      integer limesh  ! mesh mark (velocity, preasure, mhd)
!     work arrays
      real tmpv(LX1,LY1,LZ1,3), tmpvp(LX2,LY2,LZ2,3),
     $     tmpb(LBX1,LBY1,LBZ1,3), tmpbp(LBX2,LBY2,LBZ2,3)
!-----------------------------------------------------------------------
      logs = 'Local data refinement.'
      call nekp4est_log(NP4_LP_PRD,logs)

!     refine variables
!     coordinates; mesh 1
      limesh = 1
      call nekp4est_refine_vm(XM1,limesh,tmpv,LX1,LY1,LZ1,LELT)
      call nekp4est_refine_vm(YM1,limesh,tmpv,LX1,LY1,LZ1,LELT)
      if(IF3D) call nekp4est_refine_vm(ZM1,limesh,tmpv,LX1,LY1,LZ1,LELT)

!     set every step by setprop in nek_advance
!     interpolate density (VTRANS) and viscosity (VDIFF)
!      fs = 2
!      fe = NFIELD
!      if (IFFLOW) fs = 1
!      if (IFMHD) fe = fe + 1
!      do il=fs, fe
!         call nekp4est_refine_vm(VDIFF (1,1,1,1,il),limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!         call nekp4est_refine_vm(VTRANS(1,1,1,1,il),limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!      enddo


!     velocity
      call nekp4est_refine_vm(VX,limesh,tmpv,LX1,LY1,LZ1,LELV)
      call nekp4est_refine_vm(VY,limesh,tmpv,LX1,LY1,LZ1,LELV)
      if(IF3D) call nekp4est_refine_vm(VZ,limesh,tmpv,LX1,LY1,LZ1,LELV)

!     lag arrays velocity
      if (IFFLOW.or.IFMHD) then
         do il =1, NBDINP-1
            call nekp4est_refine_vm(VXLAG(1,1,1,1,il),limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            call nekp4est_refine_vm(VYLAG(1,1,1,1,il),limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
            if (IF3D) call nekp4est_refine_vm(VZLAG(1,1,1,1,il),limesh,
     $              tmpv,LX1,LY1,LZ1,LELV)
         enddo
      endif

!     arrays for time integration
      if (IFTRAN) then
         call nekp4est_refine_vm(ABX1,limesh,tmpv,LX1,LY1,LZ1,LELV)
         call nekp4est_refine_vm(ABY1,limesh,tmpv,LX1,LY1,LZ1,LELV)
         if(IF3D) call nekp4est_refine_vm(ABZ1,limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)

         call nekp4est_refine_vm(ABX2,limesh,tmpv,LX1,LY1,LZ1,LELV)
         call nekp4est_refine_vm(ABY2,limesh,tmpv,LX1,LY1,LZ1,LELV)
         if(IF3D) call nekp4est_refine_vm(ABZ2,limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)
      endif

!     mesh 2
      limesh = 2
!     preasure
      call nekp4est_refine_vm(PR,limesh,tmpvp,LX2,LY2,LZ2,LELV)

!     lag arrays preasure PRLAG
      if (NBDINP.eq.3) call nekp4est_refine_vm(PRLAG,limesh,tmpvp,
     $              LX2,LY2,LZ2,LELV)

!     mesh 2
      limesh = 2
!     user defined divergence USRDIV
      call nekp4est_refine_vm(USRDIV,limesh,tmpvp,LX2,LY2,LZ2,LELT)

!     temperature and passive scalars
!     mesh 1
      limesh = 1
      if (IFHEAT) then
!     varialbe T; temperature nad passive scalars
         do il=2,NFIELD
            call nekp4est_refine_vm(T(1,1,1,1,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)

!     lag arrays TLAG
            do jl =1, NBDINP-1
               call nekp4est_refine_vm(TLAG(1,1,1,1,jl,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)
            enddo

!     arrays time integration VGRADT[12]
            call nekp4est_refine_vm(VGRADT1(1,1,1,1,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)
            call nekp4est_refine_vm(VGRADT2(1,1,1,1,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)
         enddo                  ! passive scalar loop
      endif                     ! IFHEAT

!     mhd
!     NOT TESTED AND NOT SURE IT IS OK!!!!!!!!!!!!!!!!
      if (IFMHD) then
!     mesh 3
         limesh = 3
!     magnetic field
         call nekp4est_refine_vm(BX,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         call nekp4est_refine_vm(BY,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         if(IF3D) call nekp4est_refine_vm(BZ,limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)

!     arrays for time integration
         call nekp4est_refine_vm(BBX1,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         call nekp4est_refine_vm(BBY1,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         if(IF3D) call nekp4est_refine_vm(BBZ1,limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)

         call nekp4est_refine_vm(BBX2,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         call nekp4est_refine_vm(BBY2,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         if(IF3D) call nekp4est_refine_vm(BBZ2,limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)

!     lag arrays
         do il =1, NBDINP-1
            call nekp4est_refine_vm(BXLAG(1,il),limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_refine_vm(BYLAG(1,il),limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)
            if (IF3D) call nekp4est_refine_vm(BZLAG(1,il),limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)
         enddo

!     mesh 4
         limesh = 4
         call nekp4est_refine_vm(PM,limesh,tmpbp,LBX2,LBY2,LBZ2,LBELV)

!     lag arrays
         if (NBDINP.eq.3) call nekp4est_refine_vm(PMLAG,limesh,tmpbp,
     $              LBX2,LBY2,LBZ2,LBELV)

      endif

!     moving mesh; BM1LAG
      if (IFGEOM) then
         logs = 'Error: nekp4est_refine_local_v no moving mesh yet.'
         call nekp4est_abort(logs)
      endif

!     it could be place for perturbation, but interpolation would not
!     give correct base flow structure, so I skip it for now


      return
      end
!=======================================================================
!> @brief Local corsening
!! @details Threre are two possible implementations. In following rourtine
!! we performe all block coarsening for a single variable implemented in
!! nekp4est_coarse_vm and next turn to the variable. This version could
!! be slightly faster.
      subroutine nekp4est_coarsen_local_v
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'           ! IF3D
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'           ! NBDINP
      include 'GEOM_DEF'
      include 'GEOM'            ! IFGEOM
      include 'NEKP4EST'

!     local variables
      character(len=NP4_LSTL_LOG) logs  ! log string
      integer il, jl ! loop index
      integer fs, fe
      integer limesh  ! mesh mark (velocity, preasure, mhd)
!     work arrays
      real tmpv(LX1,LY1,LZ1,3), tmpvp(LX2,LY2,LZ2,3),
     $     tmpb(LBX1,LBY1,LBZ1,3), tmpbp(LBX2,LBY2,LBZ2,3)
!-----------------------------------------------------------------------
      logs = 'Local data coarsening.'
      call nekp4est_log(NP4_LP_PRD,logs)

!     coarsen variables
!     coordinates; mesh 1
      limesh = 1
      call nekp4est_coarse_vm(XM1,limesh,tmpv,LX1,LY1,LZ1,LELT)
      call nekp4est_coarse_vm(YM1,limesh,tmpv,LX1,LY1,LZ1,LELT)
      if(IF3D) call nekp4est_coarse_vm(ZM1,limesh,tmpv,
     $              LX1,LY1,LZ1,LELT)

!     set every step by setprop in nek_advance
!     interpolate density (VTRANS) and viscosity (VDIFF)
!      fs = 2
!      fe = NFIELD
!      if (IFFLOW) fs = 1
!      if (IFMHD) fe = fe + 1
!      do il=fs, fe
!         call nekp4est_coarse_vm(VDIFF (1,1,1,1,il),limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!         call nekp4est_coarse_vm(VTRANS(1,1,1,1,il),limesh,
!     $                           tmpv,LX1,LY1,LZ1,LELT)
!      enddo

!     velocity
      call nekp4est_coarse_vm(VX,limesh,tmpv,LX1,LY1,LZ1,LELV)
      call nekp4est_coarse_vm(VY,limesh,tmpv,LX1,LY1,LZ1,LELV)
      if(IF3D) call nekp4est_coarse_vm(VZ,limesh,tmpv,LX1,LY1,LZ1,LELV)

!     lag arrays velocity
      if (IFFLOW.or.IFMHD) then
         do il =1, NBDINP-1
            call nekp4est_coarse_vm(VXLAG(1,1,1,1,il),limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)
            call nekp4est_coarse_vm(VYLAG(1,1,1,1,il),limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)
            if (IF3D) call nekp4est_coarse_vm(VZLAG(1,1,1,1,il),
     $              limesh,tmpv,LX1,LY1,LZ1,LELV)
         enddo
      endif

!     arrays for time integration
      if (IFTRAN) then
         call nekp4est_coarse_vm(ABX1,limesh,tmpv,LX1,LY1,LZ1,LELV)
         call nekp4est_coarse_vm(ABY1,limesh,tmpv,LX1,LY1,LZ1,LELV)
         if(IF3D) call nekp4est_coarse_vm(ABZ1,limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)

         call nekp4est_coarse_vm(ABX2,limesh,tmpv,LX1,LY1,LZ1,LELV)
         call nekp4est_coarse_vm(ABY2,limesh,tmpv,LX1,LY1,LZ1,LELV)
         if(IF3D) call nekp4est_coarse_vm(ABZ2,limesh,tmpv,
     $              LX1,LY1,LZ1,LELV)
      endif

!     mesh 2
      limesh = 2
!     preasure
      call nekp4est_coarse_vm(PR,limesh,tmpvp,LX2,LY2,LZ2,LELV)

!     lag arrays preasure PRLAG
      if (NBDINP.eq.3) call nekp4est_coarse_vm(PRLAG,limesh,tmpvp,
     $              LX2,LY2,LZ2,LELV)


!     mesh 2
      limesh = 2
!     user defined divergence USRDIV
      call nekp4est_coarse_vm(USRDIV,limesh,tmpvp,LX2,LY2,LZ2,LELT)

!     temperature and passive scalars
!     mesh 1
      limesh = 1
      if (IFHEAT) then
!     varialbe T; temperature nad passive scalars
         do il=2,NFIELD
            call nekp4est_coarse_vm(T(1,1,1,1,il-1),limesh,tmpv,
     $              LX1,LY1,LZ1,LELT)

!     lag arrays TLAG
            do jl =1, NBDINP-1
               call nekp4est_coarse_vm(TLAG(1,1,1,1,jl,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)
            enddo

!     arrays time integration VGRADT[12]
            call nekp4est_coarse_vm(VGRADT1(1,1,1,1,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)
            call nekp4est_coarse_vm(VGRADT2(1,1,1,1,il-1),limesh,
     $              tmpv,LX1,LY1,LZ1,LELT)
         enddo                  ! passive scalar loop
      endif                     ! IFHEAT

!     mhd
!     NOT TESTED AND NOT SURE IT IS OK!!!!!!!!!!!!!!!!
      if (IFMHD) then
!     mesh 3
         limesh = 3
!     magnetic field
         call nekp4est_coarse_vm(BX,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         call nekp4est_coarse_vm(BY,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         if(IF3D) call nekp4est_coarse_vm(BZ,limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)

!     arrays for time integration
         call nekp4est_coarse_vm(BBX1,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         call nekp4est_coarse_vm(BBY1,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         if(IF3D) call nekp4est_coarse_vm(BBZ1,limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)

         call nekp4est_coarse_vm(BBX2,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         call nekp4est_coarse_vm(BBY2,limesh,tmpb,LBX1,LBY1,LBZ1,LBELV)
         if(IF3D) call nekp4est_coarse_vm(BBZ2,limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)

!     lag arrays
         do il =1, NBDINP-1
            call nekp4est_coarse_vm(BXLAG(1,il),limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)
            call nekp4est_coarse_vm(BYLAG(1,il),limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)
            if (IF3D) call nekp4est_coarse_vm(BZLAG(1,il),limesh,tmpb,
     $              LBX1,LBY1,LBZ1,LBELV)
         enddo

!     mesh 4
         limesh = 4
         call nekp4est_coarse_vm(PM,limesh,tmpbp,LBX2,LBY2,LBZ2,LBELV)

!     lag arrays
         if (NBDINP.eq.3) call nekp4est_coarse_vm(PMLAG,limesh,tmpbp,
     $              LBX2,LBY2,LBZ2,LBELV)
      endif

!     moving mesh; BM1LAG
      if (IFGEOM) then
         logs = 'Error: nekp4est_coarse_loca_v no moving mesh yet.'
         call nekp4est_abort(logs)
      endif

!     it could be place for perturbation, but interpolation would not
!     give correct base flow structure, so I skip it for now


      return
      end
!=======================================================================
!> @brief Divide by BM1 before refining/coarsening
      subroutine nekp4est_remove_bm
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'MASS_DEF'
      include 'MASS'

!     local variables
      integer ntotv,ntott,il
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
!     arrays for time integration
      if (IFTRAN) then
        if (IFFLOW) then
            call invcol2(ABX1,BM1,ntotv)
            call invcol2(ABY1,BM1,ntotv)
            if(IF3D) call invcol2(ABZ1,BM1,ntotv)

            call invcol2(ABX2,BM1,ntotv)
            call invcol2(ABY2,BM1,ntotv)
            if(IF3D) call invcol2(ABZ2,BM1,ntotv)
        endif


!     temperature and passive scalars
        if (IFHEAT) then
!     varialbe T; temperature nad passive scalars
            do il=2,NFIELD
!     arrays time integration VGRADT[12]
                call invcol2(VGRADT1(1,1,1,1,il-1),BM1,ntott)
                call invcol2(VGRADT2(1,1,1,1,il-1),BM1,ntott)
            enddo                 ! passive scalar loop
        endif                     ! IFHEAT
      endif                       ! IFTRAN

      return
      end
!=======================================================================
!> @brief Multiply by BM1 after refining/coarsening
      subroutine nekp4est_mult_bm
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'MASS_DEF'
      include 'MASS'

!     local variables
      integer ntotv,ntott,il
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
!     arrays for time integration
      if (IFTRAN) then
        if (IFFLOW) then
            call col2(ABX1,BM1,ntotv)
            call col2(ABY1,BM1,ntotv)
            if(IF3D) call invcol2(ABZ1,BM1,ntotv)

            call col2(ABX2,BM1,ntotv)
            call col2(ABY2,BM1,ntotv)
            if(IF3D) call invcol2(ABZ2,BM1,ntotv)
        endif

!     temperature and passive scalars
        if (IFHEAT) then
!     varialbe T; temperature nad passive scalars
            do il=2,NFIELD
!     arrays time integration VGRADT[12]
                call col2(VGRADT1(1,1,1,1,il-1),BM1,ntott)
                call col2(VGRADT2(1,1,1,1,il-1),BM1,ntott)
            enddo                 ! passive scalar loop
        endif                     ! IFHEAT
      endif                       ! IFTRAN

      return
      end
!=======================================================================
