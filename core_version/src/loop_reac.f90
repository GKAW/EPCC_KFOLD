! ==============================================================================
! Subroutine: LOOP_REAC (R,INDX)
! 
! Purpose: Recomputes the possible reactions within a loop element and
!          their corrisponding rates. 
!
! Method:  There are 6 possible reactions around the loop:
!
!          (1) Nucleation                   A
!                                          A A
!              A A A A A A U A A  -->  A A A-U A A
!
!          (2) Helix Extension
!                                        x-x
!              A A A x-x U A A  -->  A A A-U A A
!
!          (3) Helix Retraction
!                    x-x
!                A A A-U A A   -->  A A A x-x U A A
!
!          (4) Helix Morphing               x-x
!                  x-x  x-x                 U-A
!                A A-U  U-A A  -->  A A x-x U-A
!
!          (5) Defect Diffusion         x-x
!                    x-x                  U
!                A A A-U U A A -->  A A A-U A A
!
!          (6) Helix Open
!                         U-A       U-A
!                         x-x      x   x
!                         A-U  -->  A-U
!
! Arguments:
!
!             R - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!          INDX - The indx number of the loop element that reactions
!                 will be calculated for.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================


! Safe variables:
! write-only:
! wrk1
! wrk2

! read-only:
! r% ibsp  : for subroutines too, as there is no changing of structure before reaction fires.
! r% loop

! issues:
! atot : iteratively updated, looking at making allocatable and summing at the end 


      SUBROUTINE LOOP_REAC (R,INDX)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        INTEGER, INTENT(IN) :: indx

        !=== VARIABLES ===!
        INTEGER :: i,j,k,n,ip,jp,kp,is,js, alloci
        INTEGER :: ks,ke,l,lmx,icnt,iloop
        INTEGER :: nt,nh,ns,mt,mh,ms,icase

        INTEGER, DIMENSION(:), ALLOCATABLE :: klist, cumcnt

        DOUBLE PRECISION :: x,dg,atot,rate
        ! apartial holds threadwise components of atot        
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: apartial

        !write(*,*) "loop_reac allocatables declared"
        i = r% loop(indx)
        j = r% ibsp(i)
        n = r% n

        IF ( i == n ) j = 1

        nh = r% nhlx(indx)
        ns = r% nsgl(indx)

        nt = ns + 2 * nh

        !=== Internal Loop iloop = 1 ===!
        !=== External Loop iloop = 0 ===!

        IF ( i < j ) iloop = 1
        IF ( i > j ) iloop = 0

        IF ( iloop == 1 ) THEN
          ks = i + 1
          ke = j
        ELSE
          ks = j
          ke = i
        ENDIF


        !=== COMPUTE REACTIONS ===!

        icnt = 0

        k = ks
        !=== Trace loop once to identify allocatable size ===!
        DO WHILE (k <= ke)
          IF ( r% ibsp(k) > 0 ) THEN
            IF ( k /= ke ) THEN
              k = r% ibsp(k)
            ENDIF
          ENDIF

          k = k + 1
          icnt = icnt + 1
        ENDDO

        !=== Allocate arrays ===!
        ALLOCATE( apartial(icnt))
        ALLOCATE( klist(icnt))
        ALLOCATE( cumcnt(icnt))
        !write(*,*) "loop_reac allocatables allocated: ", icnt
        !=== Initialize to zero ===!
        apartial(:) = 0.0d0
        klist(:) = 0
        cumcnt(:) = 0

        !=== Trace loop a second time to populate shared allocatables ===!
        k = ks
        DO alloci=1,icnt
          klist(alloci) = k
          IF ( r% ibsp(k) > 0 ) THEN
            IF ( k /= ke ) THEN
              k = r% ibsp(k)
              cumcnt(alloci+1:icnt) = cumcnt(alloci+1:icnt) + 1
            ENDIF
          ENDIF
          k = k + 1
          cumcnt(alloci+1:icnt) = cumcnt(alloci+1:icnt) + 1
        ENDDO


        ! ======
        !!$omp parallel do default(none) &
        !!$omp& private(klist, cumcnt, k, x, l, nt, iloop, &
        !!$omp& kp, is, js, pnuc, iwc, ip, jp, icase, ke, nh, &
        !!$omp& n, dg, beta, rate, mh, ms, mt, icnt, lmx, ns) &
        !!$omp& shared(r, apartial)

        !write(*,*) "Calling reac_core"
        DO alloci=1,icnt
          CALL REAC_CORE (r, nt, nh, ns, ke, iloop, icnt, alloci, klist, cumcnt, apartial)
          !write(*,*) "reac_core completed: ", alloci, " / ", icnt
        ENDDO

        atot = 0.0d0

        !=== Combine partial sums ===!
        DO alloci=1,icnt
          atot = atot + apartial(alloci)
        ENDDO

        !=== Save Loop Reaction Rate ===!
        !write(*,*) "Saving loop reaction rate"
        IF ( iloop == 1 ) r% wrk1(i) = atot

        !=== Open Internal Helix BP ===!

        k = klist(icnt)
        ! Replicate helix event logic from loop above for safety and readability
        IF ( r% ibsp(k) > 0 ) THEN
          ip = k
          jp = r% ibsp(k)
          IF ( r% link(ip) == 0 ) THEN
          IF ( iloop == 1 .and. k == ke ) THEN

            is = ip + 1
            js = jp - 1

            DO WHILE ( r% link(is) == 0 )

              CALL DELTAG_HI (r,is,js,dg)

              dg = dg / 2.0d0

              x = beta * dg
              x = DEXP(-x) * rateh

              r% wrk1(is) = x
              r% wrk1(js) = 0.0d0
              !write(*,*) "writing wrk2: ", is
              r% wrk2(is) = 0.0d0
              !write(*,*) "writing wrk2: ", js
              r% wrk2(js) = 0.0d0

              atot = atot + x

              is = is + 1
              js = js - 1

            ENDDO

          ENDIF
          ENDIF
        ENDIF
        r% ptot(indx) = atot

        !write(*,*) "atot", atot
        !write(*,*) "Calling loop_resum, atot", atot
        CALL LOOP_RESUM (r,indx)
        !write(*,*) "loop_resum complete"
        DEALLOCATE( cumcnt)
        !write(*,*) "cumcnt DEALLOCATE"
        DEALLOCATE( klist)
        !write(*,*) "klist DEALLOCATE"
        DEALLOCATE( apartial)
        !write(*,*) "apartial DEALLOCATE"

        RETURN

      END SUBROUTINE LOOP_REAC
