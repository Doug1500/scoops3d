        SUBROUTINE SeedSearch(iseed,nseed,nres,nzres,zzsrchres,isrchnum,jsrchnum)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine determines new search centers in the Search lattice
!       based on centers that produced minimum FOS values on the ground
!       in the previous iteration.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called by Scoops3D
!
!        VARIABLES
!
!        delxy -- DEM grid resolution
!        fostol -- maximum percent difference in FOS between seed iterations to
!           not generate a new seed
!        fsmin(i,j) -- minimum factor of safety at each DEM column
!        fsminold(i,j) -- last minimum factor of safety at each DEM column used to
!           compare change in FOS during seed iterations to see if criterion is met
!        fsx(nx,ny),fsy(nx,ny),fsz(nx,ny) -- search sphere center of minimum FOS 
!           scoop at each DEM column
!        idem,jdem -- DEM grid loop counters
!        ismin,ismax -- minimum and maximum boundaries of search grid on x-axis,
!           specified relative to DEM cells i=1 through nx
!        isrch,jsrch,ksrch -- search lattice loop counters
!        isrchnum,jsrchnum -- the number of search lattice nodes in x and y directions
!        iseed -- grid refinement iteration count also used to fill srchlatt
!           lattice. iseed=1 represents initial coarse grid locations.
!        jsmin,jsmax -- minimum and maximum boundaries of search grid on y-axis,
!           specified relative to DEM cells j=1 through ny
!        ksmax -- maximum boundary of search grid on z-axis calculated from 
!           user-specified zsmax and z resolution
!        multres -- resolution multiplier for coarse search lattice
!        nkseed(nxsrch,nysrch) -- number of k values for current seed at each i,j
!           search node
!        nsrchres -- resolution of finest search grid relative to DEM resolution.
!         Acts as a multiplier of DEM resolution delxy
!        nx -- number of DEM cells in x direction
!        ny -- number of DEM cells in y direction
!        rnull -- positive real null value used for initializations and set 
!           in CommonData. 
!        srchlatt(isrchnum,jsrchnum,ksmax) -- search lattice array with seed
!           iteration number at each searched node. 
!        xdem,ydem -- x and y locations of column lower left node
!        zdem(nx,ny) -- DEM elevations
!        zsmin,zsmax -- minimum and maximum elevations of search grid
!        zsrch -- search elevation at new search center
!        zsrchlatt(isrchnum,jsrchnum,ksmax) -- search lattice array with 
!           elevations of new search centers generated from seed.
!        zsrchres -- resolution of finest search grid on z-axis
!        zzsrchres -- resolution of current search grid on z-axis
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: nx,ny,xcen,ycen,zcen,zdem,delxy
        USE SearchData, ONLY: zsmin,zsmax,zsrchres,zsrchlatt,srchlatt,ksmax,&
                          nkseed,fostol,ismin,ismax,jsmin,jsmax,nsrchres,multres
        USE FosData, ONLY: fsmin,fsx,fsy,fsz,fsminold
        
        IMPLICIT NONE

        INTEGER(int2), INTENT(inout) :: iseed
        INTEGER, INTENT(out) :: nseed
        INTEGER, INTENT(in) :: nres,nzres,isrchnum,jsrchnum
        INTEGER :: idem,jdem,isrch,jsrch,ksrch,iisrch,jjsrch,iin,jin
        INTEGER :: iii,jjj,kkk,icen,jcen,k,kminarray(1),kmin
        REAL(pr) :: fsdiff,zsrch
        REAL(pr), INTENT (in) :: zzsrchres
        REAL, ALLOCATABLE :: sorter(:)

        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
        
        errmessage = ' '
        solution = ' '
        problemtype = 'identifying nodes for search refinement'         
        ALLOCATE (sorter(ksmax))
        sorter = REAL(zsmax + zsrchres)
        zsrchlatt = REAL(zsmax + zsrchres)
        nkseed = 0
        nseed = 0
        iseed = iseed + 1_int2
        
!     Loop through DEM cells to find min FOS search center at each grid cell.
!     After initial coarse search (iseed=2), assign 1 to all search centers 
!     associated with minimum FOS on the ground.        
        IF (iseed.eq.2_int2) THEN
          fsdiff = nullhi
          DO jdem = 1,ny
            DO idem = 1,nx
              IF (fsmin(idem,jdem).lt.nullhi) THEN
                IF (fsx(idem,jdem).gt.0) THEN
                 iisrch = INT(fsx(idem,jdem)/delxy) + 1
                ELSE
                 iisrch = INT(fsx(idem,jdem)/delxy)
                END IF 
                IF (fsy(idem,jdem).gt.0) THEN 
                 jjsrch = INT(fsy(idem,jdem)/delxy) + 1
                ELSE
                 jjsrch = INT(fsy(idem,jdem)/delxy)
                END IF 
                isrch = NINT(float(iisrch-ismin)/float(nsrchres))+1
                jsrch = NINT(float(jjsrch-jsmin)/float(nsrchres))+1
                ksrch = NINT((fsz(idem,jdem)-zsmin)/zsrchres) + 1
                IF (isrch.lt.1.or.isrch.gt.isrchnum.or.&
                  jsrch.lt.1.or.jsrch.gt.jsrchnum.or.&
                  ksrch.lt.1.or.ksrch.gt.ksmax) THEN
                   errmessage = 'search lattice array out of bounds'
                   Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                     
                END IF
                srchlatt(isrch,jsrch,ksrch) = 1_int2
                fsminold(idem,jdem) = fsmin(idem,jdem)
              END IF
            END DO
          END DO
        END IF  !  iseed = 2
         
!    Check all search nodes adjacent to current seed. If they are valid nodes
!    and have not been searched before, determine whether to place new seed  
!    for next search iteration.         
        DO jdem = 1,ny
          DO idem = 1,nx
            IF (fsmin(idem,jdem).lt.nullhi) THEN
!    Compare DEM cell factor of safety to value from last iteration.
              IF (iseed.gt.2_int2) THEN
                fsdiff = (fsminold(idem,jdem) - fsmin(idem,jdem))/fsminold(idem,jdem)
                fsminold(idem,jdem) = fsmin(idem,jdem)
              END IF
!     If % change in FOS is larger than tolerance (or this is first iteration 
!     after initial coarse iteration) assign next seed number for next iteration.
              IF (fsdiff.ge.fostol) THEN
                IF (fsx(idem,jdem).gt.0) THEN               
                  iisrch = INT(fsx(idem,jdem)/delxy) + 1
                ELSE 
                  iisrch = INT(fsx(idem,jdem)/delxy)
                END IF
                IF (fsy(idem,jdem).gt.0) THEN 
                  jjsrch = INT(fsy(idem,jdem)/delxy) + 1
                ELSE
                  jjsrch = INT(fsy(idem,jdem)/delxy)
                END IF 
                isrch = NINT(float(iisrch-ismin)/float(nsrchres))+1
                jsrch = NINT(float(jjsrch-jsmin)/float(nsrchres))+1
                ksrch = NINT((fsz(idem,jdem)-zsmin)/zsrchres) + 1 
                IF (isrch.lt.1.or.isrch.gt.isrchnum.or.&
                  jsrch.lt.1.or.jsrch.gt.jsrchnum.or.&
                  ksrch.lt.1.or.ksrch.gt.ksmax) THEN
                   errmessage = 'search lattice array out of bounds'
                   Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                    
                END IF
!     Determine if nodes around last seed have been searched and if not,
!     assign new seed to those nodes.
                If (srchlatt(isrch,jsrch,ksrch).eq.iseed-1_int2) THEN  
                  DO kkk = ksrch-nzres,ksrch+nzres,nzres
                    DO jjj = jsrch-multres,jsrch+multres,multres
                      DO iii = isrch-multres,isrch+multres,multres
                        IF ((iii.ge.1.and.iii.le.isrchnum).and. &
                            (jjj.ge.1.and.jjj.le.jsrchnum).and. &
                            (kkk.ge.1.and.kkk.le.ksmax)) THEN
                          IF (srchlatt(iii,jjj,kkk).eq.0_int2) THEN
                            IF (iii.eq.isrch) THEN
                              iin = iisrch
                            ELSE
                              IF (iii.eq.isrch-multres) THEN
                                iin = iisrch - nres
                              ELSE
                                iin = iisrch + nres
                              END IF
                            END IF
                            IF (jjj.eq.jsrch) THEN
                              jin = jjsrch
                            ELSE
                              IF (jjj.eq.jsrch-multres) THEN
                                jin = jjsrch - nres
                              ELSE
                                jin = jjsrch + nres
                              END IF
                            END IF
                            IF (kkk.eq.ksrch) THEN
                              zsrch = fsz(idem,jdem)
                            ELSE
                              IF (kkk.eq.ksrch-nzres) THEN
                                zsrch = fsz(idem,jdem) - zzsrchres
                              ELSE
                                zsrch = fsz(idem,jdem) + zzsrchres
                              END IF
                            END IF
                            IF (zsrch.le.zsmax) THEN
                              IF (iii.ge.1.and.jjj.ge.1.and.iii.le.nx.and.jjj.le.ny) THEN
                                IF (zdem(iin,jin).eq.rnull.or.zsrch.gt.zdem(iin,jin)) THEN
                                  srchlatt(iii,jjj,kkk) = iseed
                                  nkseed(iii,jjj) = nkseed(iii,jjj) + 1
                                  IF (nkseed(iii,jjj).gt.ksmax)  THEN
                                    errmessage = 'nkseed array exceeded bounds'
                                    Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                                    
                                  END IF
                                  zsrchlatt(iii,jjj,nkseed(iii,jjj)) = REAL(zsrch)
                                  nseed = nseed + 1                                   
                                END IF
                              ELSE
                                srchlatt(iii,jjj,kkk) = iseed
                                nkseed(iii,jjj) = nkseed(iii,jjj) + 1
                                IF (nkseed(iii,jjj).gt.ksmax)  THEN
                                    errmessage = 'nkseed array exceeded bounds'
                                    Call WriteError(1,errmessage,problemtype,'no','no ',0,' ')                                    
                                END IF
                                zsrchlatt(iii,jjj,nkseed(iii,jjj)) = REAL(zsrch)
                                nseed = nseed + 1
                              END IF
                            END IF
                          END IF
                        END IF
                      END DO  ! iii
                    END DO  ! jjj  
                  END DO  ! kkk
                END IF  !  If (srchlatt(isrch,jsrch,ksrch).eq.iseed-1)
              END IF  !  IF (fsdiff.ge.fostol)
            END IF  !  IF (fsmin(idem,jdem).lt.nullhi)
          END DO  ! idem
        END DO  !  jdem
         
        DO jcen = 1,jsrchnum
          DO icen = 1,isrchnum
            IF (nkseed(icen,jcen).lt.1) CYCLE
            DO k = 1,nkseed(icen,jcen)           
              kminarray = MINLOC(zsrchlatt(icen,jcen,1:nkseed(icen,jcen)))
              kmin = kminarray(1)
              sorter(k) = zsrchlatt(icen,jcen,kmin)
              zsrchlatt(icen,jcen,kmin) = REAL(zsmax + zsrchres)
            END DO
            DO k = 1,nkseed(icen,jcen)           
              zsrchlatt(icen,jcen,k) = sorter(k)
            END DO
          END DO
        END DO
         
        DEALLOCATE (sorter)
         
        END SUBROUTINE SeedSearch
