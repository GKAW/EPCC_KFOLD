! ==============================================================================
! Module: CLASS_RNAFOLD
! 
! Purpose: A FORTRAN class structure containing subroutines and data
!          elements required for computing transition probabilites between
!          different RNA secondary structures.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules - RNAVar
! Functions -
! Subroutines - 
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      MODULE CLASS_RNAFOLD

        USE RNAVar, ONLY : rateh,ratem,rated,beta,pnuc,&
                         & iwc,eaup,mxnt

        IMPLICIT NONE

        PRIVATE

        PUBLIC :: LOOP_INIT, LOOP_FIRE

        TYPE, PUBLIC :: RNA_STRUC

          CHARACTER, DIMENSION(:), ALLOCATABLE  :: seq

          INTEGER, DIMENSION(:), ALLOCATABLE  :: iseq
          INTEGER, DIMENSION(:), ALLOCATABLE  :: ibsp
          INTEGER, DIMENSION(:), ALLOCATABLE  :: link

          INTEGER, DIMENSION(:), ALLOCATABLE  :: loop
          INTEGER, DIMENSION(:), ALLOCATABLE  :: nhlx
          INTEGER, DIMENSION(:), ALLOCATABLE  :: nsgl

          INTEGER :: n
          INTEGER :: nl
          INTEGER :: nsum

          DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: psum
          DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ptot

          DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: wrk1
          DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: wrk2

        END TYPE RNA_STRUC

        CONTAINS

        INCLUDE 'loop_init.f90'
        INCLUDE 'loop_resum.f90'

        INCLUDE 'helx_reac.f90'
        INCLUDE 'reac_core.f90'
        INCLUDE 'loop_reac.f90'
        INCLUDE 'loop_fire.f90'

      END MODULE CLASS_RNAFOLD
