        SUBROUTINE iwater4wt
             
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        This subroutine calculates unit weights at 3D pressure file
!        locations for columns with an unsaturated portion (iwater=4).
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       Called in Readin
!
!       VARIABLES
!

!       ddepth -- depth between last and current unit weight calculation elevations
!       depth -- depth from DEM to pressure data elevations
!       diff -- distance from pressure node to layer boundary
!       dsedz(nx,ny,pmaxk) -- array of effective saturation gradients
!       dthetazdz(nx,ny,pmaxk) -- array of volumetric water content gradients
!       duwtdz(nx,ny,strmaxk) -- array of vertical unit weight gradients
!       gamr(nmat,3) -- array of total, partially saturated, and saturated
!           unit weights for each layer 
!       i,j -- DEM grid array location
!       ktop -- max k value at particular i,j
!       l -- Material layer counter
!       layer(nmat,nx,ny) -- array of bottom elevations for material layers
!       maxpk(nx,ny) -- array of highest k value of pressure data at each DEM cell
!       pdepth -- depth between pressure data elevations
!       pmaxk -- number of 3-d pressure values at each cell
!       presshpz(nx,ny,pmaxk) -- array of z elevations of pressure head data
!       rnull -- real null value used throughout program 
!       se(nx,ny,pmaxk) -- array of effective saturation 
!       thetares(nmat) -- residual water content for each layer
!       thetasat(nmat) -- saturated water content for each layer
!       thetaz(nx,ny,pmaxk) -- array of volumetric water content
!       tres -- residual water content for current layer
!       tsat -- saturated water content for current layer
!       uwt3d(nx,ny,strmaxk) -- unit weight at each depth (depth weighted average)
!       uwtsat -- saturated unit weight for current layer
!       uwtsum -- sum of unit weights weighted for depth
!       zdem(nx,ny) -- DEM elevations
!       zlast -- elevation of previous data point or layer boundary
!
!       OUTPUT FILES
!           20 'inputfilename_out.txt' -- echo of input and results containing
!                overall minimum slip surface data, written in subroutines
!                Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!                Readstrength, and Writeout.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: nx,ny,zdem
        USE MaterialData, ONLY: nmat,gamw,gamr,thetasat,thetares,layer
        USE WaterData, ONLY: pmaxk,maxpk,thetaz,dthetazdz,se,dsedz,presshpz
        USE StrengthData, ONLY: uwt3d,duwtdz
        
        IMPLICIT NONE
        
        INTEGER :: ktop,i,j,k,l,error,last
        REAL(pr) :: zlast,ddepth,pdepth,uwtsum,depth,diff
        REAL(pr) :: zlay(nmat),tsat,tres,uwtz,uwtzlast,uwtsat

        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'calculating unsaturated unit weight'          

        ALLOCATE (uwt3d(nx,ny,pmaxk),STAT=error) 
        ALLOCATE (duwtdz(nx,ny,pmaxk),STAT=error)
        ALLOCATE (se(nx,ny,pmaxk),STAT=error) 
        ALLOCATE (dsedz(nx,ny,pmaxk),STAT=error)
        IF (error.ne.0) THEN
          errmessage = 'unit weight arrays not allocated successfully' 
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF
        uwt3d = rnull
        duwtdz = 0.0_pr
        se = rnull
        duwtdz = 0.0_pr

           
!     Calculate depth-weighted unit weights and gradients at each pressure head location.
        DO j=ny,1,-1
          DO i=nx,1,-1
            IF (zdem(i,j).eq.rnull) CYCLE      
            IF (nmat.gt.1) THEN
              DO l=1,nmat-1
                zlay(l) = layer(l,i,j)
              END DO
            END IF 
!    Set elevation of lowest layer to be lower than 3D pressure current node.             
            zlay(nmat) = presshpz(i,j,1)-1.0_pr          
            uwtsum = 0.0_pr
            ktop = maxpk(i,j) 
            last = 1
            pdepth = 0.0_pr
            depth = 0.0_pr
            uwtsat = gamr(1,3)
            tsat = thetasat(1)
            tres = thetares(1) 
            zlast = zdem(i,j)
            uwtzlast = uwtsat - (tsat-thetaz(i,j,1))*gamw  
                         
            DO k = ktop,1,-1
              IF (k.lt.ktop) THEN
                depth = zdem(i,j) - presshpz(i,j,k)
                pdepth = presshpz(i,j,k+1)-presshpz(i,j,k) 
              END IF
!    Find layer at current depth to determine saturated and residual water contents              
              DO l = last,nmat
                IF (zlay(l).ne.rnull) THEN
                  uwtsat = gamr(l,3)
                  tsat = thetasat(l)
                  tres = thetares(l)
                  IF (presshpz(i,j,k).lt.zlay(l)) THEN  ! If 3D node below this layer
                    ddepth = zlast - zlay(l)
                    diff = zlay(l) - presshpz(i,j,k)
!    Find water content at layer boundary using interpolation of 3D values and add to
!    running sum of column wt.                    
                    uwtz = uwtsat - (tsat-(thetaz(i,j,k)+ dthetazdz(i,j,k)*diff))*gamw
                    uwtsum = uwtsum + ((uwtz+uwtzlast)/2.0_pr) * ddepth
                    zlast = zlay(l)
                    uwtzlast = uwtz
                  ELSE  !   Node is in this layer
                    IF (thetaz(i,j,k).gt.tsat) THEN
                      thetaz(i,j,k) = tsat  
!  Note - don't use WriteError for this warning because we do not want the message written to terminal output                      
!                      PRINT *,'*** WATER CONTENT GREATER THAN SATURATED WATER CONTENT AT NODE: ',i,j,k
!                      PRINT *,'*** WILL SET EQUAL TO RESIDUAL WATER CONTENT AND CONTINUE'
                      WRITE (39,*) 'WARNING -- water content greater than saturated water content at node: ',i,j,k 
                      WRITE (39,*) '   Scoops3D will set water content equal to saturated water content'  
                    END IF
                    IF (thetaz(i,j,k).lt.tres) THEN
                      thetaz(i,j,k) = tres
!                      PRINT *,'*** WATER CONTENT LESS THAN RESIDUAL WATER CONTENT AT NODE: ',i,j,k
!                      PRINT *,'*** WILL SET EQUAL TO RESIDUAL WATER CONTENT AND CONTINUE'
                      WRITE (39,*) ' WARNING -- water content less than residual water content at node: ',i,j,k 
                      WRITE (39,*) '  Scoops3D will set water content equal to residual water content'  
                    END IF                     
                    ddepth = zlast - presshpz(i,j,k)
                    se(i,j,k) = (thetaz(i,j,k) - tres)/(tsat - tres)
                    uwtz = uwtsat - (tsat-thetaz(i,j,k))*gamw
                    uwtsum = uwtsum + ((uwtz+uwtzlast)/2.0_pr) * ddepth
                    IF (depth.gt.0.0_pr) THEN
                      uwt3d(i,j,k) = uwtsum/depth
                    ELSE
                      IF (k.lt.ktop) THEN
                        uwt3d(i,j,k) = uwt3d(i,j,k+1)
                      ELSE
                        uwt3d(i,j,k) = uwtz
                      END IF
                    END IF
                    IF (pdepth.eq.0.0_pr) THEN  !  Set gradient to zero for top node.
                      duwtdz(i,j,k) = 0.0_pr
                      dsedz(i,j,k) = 0.0_pr
                    ELSE
                      duwtdz(i,j,k) = (uwt3d(i,j,k+1) - uwt3d(i,j,k))/pdepth
                      dsedz(i,j,k) = (se(i,j,k+1) - se(i,j,k))/pdepth
                    END IF
                    uwtzlast = uwtz
                    last = l
                    EXIT
                  END IF
                END IF
              END DO                                   
            END DO
          END DO
        END DO

        IF (ALLOCATED(dthetazdz)) DEALLOCATE (dthetazdz)

        RETURN
        
1000    FORMAT(A)

        END
        
