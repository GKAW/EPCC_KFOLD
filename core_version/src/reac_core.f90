SUBROUTINE REAC_CORE (R, NT, NH, NS, KE, ILOOP, ICNT, ALLOCI, KLIST, CUMCNT, APARTIAL)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        INTEGER, INTENT(IN) :: nt
        INTEGER, INTENT(IN) :: nh
        INTEGER, INTENT(IN) :: ns
        INTEGER, INTENT(IN) :: ke
        INTEGER, INTENT(IN) :: iloop
        INTEGER, INTENT(IN) :: icnt
        INTEGER, INTENT(IN) :: alloci

        INTEGER, INTENT(IN) :: klist(icnt)
        INTEGER, INTENT(IN) :: cumcnt(icnt)

        DOUBLE PRECISION, INTENT(INOUT) :: apartial(icnt)


        !=== VARIABLES ===!
        INTEGER :: k,n,ip,jp,kp,is,js
        INTEGER :: l,lmx
        INTEGER :: mt,mh,ms,icase

        DOUBLE PRECISION :: x,dg, rate

        n = r% n


          k = klist(alloci)
          !write(*,*) "loop_reac: alloci, k", alloci, k 
          !=== Nucleation Events ===!x

          IF ( r% ibsp(k) == 0 ) THEN

            r% wrk1(k) = 0.0d0
            r% wrk2(k) = 0.0d0

            x = 0.0d0

            l = 2
            lmx = nt / 2 + 1

            IF ( MOD(nt,2) == 0 ) THEN
            IF ( cumcnt(alloci)+1 > lmx-1 ) lmx = lmx - 1
            ENDIF

            IF ( iloop == 0 ) lmx = nt - cumcnt(alloci)

            kp = k + 1
            is = r% iseq(k)
            !write(*,*) "lmx", lmx
            DO WHILE ( l <= lmx )

              IF ( r% ibsp(kp) == 0 ) THEN

                js = r% iseq(kp)

                IF ( l > 4 .and. iwc(is,js) == 1 ) THEN
                  x = x + pnuc(l)
                ENDIF

              ELSE

                l = l + 1
                kp = r% ibsp(kp)

              ENDIF

              l = l + 1
              kp = kp + 1

            ENDDO
 
            r% wrk1(k) = x
            !write(*,*) "nucleation: apartial, alloci", x, alloci
            apartial(alloci) = apartial(alloci) + x

          ENDIF
          !write(*,*) "nucleation done"

          !=== Helix Events ===!

          IF ( r% ibsp(k) > 0 ) THEN

            ip = k
            jp = r% ibsp(k)
            r% wrk1(jp) = 0.0d0
            !write(*,*) "writing wrk2: ", ip
            r% wrk2(ip) = 0.0d0

            IF ( r% link(ip) == 0 ) THEN
            r% wrk1(ip) = 0.0d0
            !write(*,*) "writing wrk2: ", jp
            r% wrk2(jp) = 0.0d0
            ENDIF

            icase = 0

            !=== Helix Extension ===!

            IF ( ip > 1 .and. jp < n ) THEN
            IF ( nh > 1  .or. ns > 4 ) THEN
            IF ( r%ibsp(ip-1) == 0 .and. r%ibsp(jp+1) == 0 ) THEN

              is = r% iseq(ip-1)
              js = r% iseq(jp+1)

              IF ( iwc(is,js) == 1 ) THEN

                icase = 1

                IF ( iloop == 1 ) THEN
                IF ( nh == 2 .and. ns == 2 ) THEN
                  icase = 2
                  IF ( k == ke ) icase = 0
                ELSEIF ( k == ke ) THEN
                  icase = 3
                ENDIF
                ENDIF

              ENDIF

            ENDIF
            ENDIF
            ENDIF

            IF ( icase > 0 ) THEN

              CALL DELTAG_HE (r,ip,jp,dg)

              dg = dg / 2.0d0

              x = beta * dg
              x = DEXP(-x) * rateh
              !!write(*,*) "writing wrk2: ", ip
              r% wrk2(ip) = x
              apartial(alloci) = apartial(alloci) + x

            ENDIF

            icase = 0

            !=== Helix Retraction ===!

            IF ( ip /= n .and. jp /= 1 ) THEN

              IF ( r% ibsp(ip+1) == jp-1 ) THEN

                icase = 1

                IF ( iloop == 1 ) THEN
                IF ( k == ke ) icase = 2
                ENDIF

                rate = rateh

              ELSEIF ( iloop == 0 .or. k /= ke ) THEN

                icase = 3

                l  = r% link(ip)
                mh = r% nhlx(l)
                ms = r% nsgl(l)

                mt = ms + 2 * mh

                IF ( iloop == 1 ) THEN
                  l = MIN(nt,mt)
                ELSE
                  l = mt
                ENDIF

                rate = pnuc(l)

              ENDIF

            ENDIF

            IF ( icase > 0 ) THEN

              CALL DELTAG_HR (r,ip,jp,dg)

              IF ( icase /= 3 ) THEN
                dg = dg / 2.0d0
              ENDIF

              x = beta * dg
              x = DEXP(-x) * rate

              IF ( icase == 2 ) THEN
                r% wrk1(ip) = x
              ELSE
                r% wrk1(jp) = x
              ENDIF

              apartial(alloci) = apartial(alloci) + x

            ENDIF

            icase = 0

            !=== Helix Morphing ===!

            IF ( iloop == 0 .or. nh > 2 ) THEN
            IF ( ip > 1 .and. jp < n ) THEN

              is = r% iseq(ip-1)
              js = r% iseq(jp+1)

              IF ( iwc(is,js) == 1 ) icase = 1

              is = r% ibsp(ip-1)
              js = r% ibsp(jp+1)

              IF ( is /= 0 ) THEN
              IF ( r% link(is) /= 0 ) icase = 0
              ENDIF

              IF ( js /= 0 ) THEN
              IF ( r% link(jp+1) /= 0 ) icase = 0
              ENDIF

              IF ( is == 0 .and. js == 0 ) icase = 0

            ENDIF
            ENDIF

            IF ( icase > 0 ) THEN

              CALL DELTAG_HM (r,ip,jp,dg)

              dg = dg / 2.0d0

              x = beta * dg
              x = DEXP(-x) * ratem

              apartial(alloci) = apartial(alloci) + x

            ENDIF

            !=== Defect Diffusion ===!

            !=== PUSH ===!

            icase = 0

            IF ( r% link(ip) == 0 ) THEN

              icase = 2

              IF ( iloop == 1 ) THEN
              IF ( nh == 2 .and. ns == 1 ) icase = 3
              IF ( nh == 1 .and. ns == 3 ) icase = 0
              ENDIF

            ELSEIF ( iloop == 0 .or. k /= ke ) THEN

              icase = 1

              IF ( iloop == 1 ) THEN
              IF ( nh == 2 .and. ns == 1 ) icase = 4
              ENDIF

            ENDIF

            IF ( icase > 0 ) THEN

              !=== Push 5' End ===!

              kp = ip - 1

              IF ( kp >= 1 ) THEN
              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(kp)
                js = r% iseq(jp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  apartial(alloci) = apartial(alloci) + x

                ENDIF

              ENDIF
              ENDIF

              !=== Push 3' End ===!

              kp = jp + 1

              IF ( kp <= n ) THEN
              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(ip)
                js = r% iseq(kp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  apartial(alloci) = apartial(alloci) + x

                ENDIF

              ENDIF
              ENDIF

            ENDIF

            !=== PULL === !

            icase = 0

            IF ( r% link(ip) /= 0 ) THEN
            IF ( iloop == 0 .or. k /= ke ) THEN

               l = r% link(ip)
              mh = r% nhlx(l)
              ms = r% nsgl(l)

              icase = 1

              IF ( mh == 1 .and. ms == 3 ) icase = 0
              IF ( mh == 2 .and. ms == 1 ) icase = 4

            ENDIF
            ENDIF

            IF ( icase > 0 ) THEN

              !=== Pull 5' End ===!

              kp = ip + 1

              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(kp)
                js = r% iseq(jp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  apartial(alloci) = apartial(alloci) + x

                ENDIF

              ENDIF

              !=== Pull 3' End ===!

              kp = jp - 1

              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(ip)
                js = r% iseq(kp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  apartial(alloci) = apartial(alloci) + x

                ENDIF

              ENDIF

            ENDIF
            !write(*,*) "helix tings: apartial, alloci", x, alloci

            ! NOTE: Moved saving the loop propensity to the end, and compensated for Open Internal Helix BP section below.

          ENDIF

    RETURN

END SUBROUTINE REAC_CORE