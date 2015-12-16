        SUBROUTINE Readpressh(prfile,prcoords,pcount,pnum)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        This subroutine reads an ASCII data file that contains
!        3-D pressure head information.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called in Readin
!
!        VARIABLES
!
!        coords -- flag for coordinates of pressure-head file
!                  1 = ijk, 2 = ijz, 3 = xyz
!        thetaz(nx,ny,pmaxk) -- array of volumetric water content
!        dthetazdz(nx,ny,pmaxk) -- array of volumetric water content gradients
!        delxy -- DEM grid resolution (delta x, delta y)
!        delzp -- vertical resolution of pressure head data
!        dpdz(nx,ny,pmaxk) -- array of pressure head gradients
!        dzp -- elevation difference between pressure nodes, used for gradient
!        error -- file allocation flag
!        fxa(nmat) -- curve-fitting parameter, a from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!        fxn(nmat) -- curve-fitting parameter, n from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!        fxm(nmat) -- curve-fitting parameter, m from the equation for the soil-water characteristic 
!            curve defined by Fredlund and Xing (1994) 
!        fxr(nmat) -- curve-fitting parameter, psi_r (soil suction associated with residual moisture 
!            content) from the equation for the soil-water characteristic curve defined by 
!            Fredlund and Xing (1994)
!        heading -- character variable for reading input
!        headlong -- length of heading string
!        i,j -- DEM grid array location
!        ilast,jlast -- used to determine z order of data, last x and y values
!        ijorder -- flag for whether pressure data is given as all z data for 
!           each i,j point (=1), or all i,j data given for each z location (=0)
!        invert -- flag indicating whether to invert z data so lowest z is 
!           associated with lowest k.  1=invert
!        ios -- status of read calls
!        ios2 -- status of file opening calls
!        iwater -- flag indicating method for modeling water pressure.
!            0 = no water pressures, 1 = ru approximation, 2 = piezometric surface, 
!            3=input 3-d pressure file, 4=3D variably saturated file containing
!            pressure head and moisture content, 5=3D variably saturated file,
!            moisture content from vanGenuchten SWCC, 6=3D variably saturated file,
!            moisture content from Fredlund and Xing SWCC
!        ktop -- k location at or below DEM 
!        maxi - maximum i value (number of columns) in pressure-head file
!        maxj - maximum j value (number of rows) in pressure-head file
!        maxpk(nx,ny) -- array of highest k value of pressure data at each DEM cell
!        nx -- number of DEM columns
!        ny -- number of DEM rows
!        pcount -- number of DEM cells with piezometric surface found
!        pflag -- equals 1 if piezometric surface found at current DEM cell
!        pgrid -- flag for pressure head grid, 1=regular, 2=irregular grid
!        piezo(i,j) -- array of piezometric elevation at each DEM cell
!        phead -- pressure head value from input file
!        pmaxk -- specified maximum number of 3-d pressure values at each cell
!        pnum -- number of horizontal cells with pressure values
!        prcoords -- coordinate system = 'ijk','ijz' or 'xyz'
!        pressh(nx,ny,pmaxk) -- array of pressure head data
!        presshpz(nx,ny,pmaxk) -- array of z elevations for irregular grid 
!           pressure head data
!        prfile -- pressure file name
!        pzsurf -- flag for whether piezometric surface was found
!        rnull -- real null value used throughout program
!        tempp,tempz -- temporary pressure and elevation arrays used to invert 
!        vga(nmat) -- curve-fitting parameter, alpha from the equation for the 
!           soil-water characteristic curve defined by vanGenuchten (1980) 
!        vgn(nmat) -- curve-fitting parameter, n, from the equation for the soil-water 
!            characteristic curve defined by vanGenuchten (1980)
!        vmc -- volumetric moisture content
!        x,y -- input location variables
!        xll,yll -- x and y origin of DEM grid, read in header lines
!        z -- elevation of pressure data
!        zdem(i,j) -- DEM elevations
!        zlast -- elevation of previous data point
!        zpmin -- minimum elevation of pressure head values
!
!
!        INPUT FILES
!         unit filename
!           13 'inputfilename.#' -- file containing 3D pressure input
!
!        OUTPUT FILES
!         unit filename
!           20 'inputfilename_out.txt' -- echo of input and results containing
!                overall minimum slip surface data, written in subroutines
!                Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!                Readstrength, and Writeout.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY:  nx,ny,zdem,delxy,xll,yll
        USE MaterialData, ONLY: nmat,thetasat,thetares,layer,gamw,&
                             fxa,fxn,fxm,fxr,vga,vgn
        USE WaterData, ONLY: pzsurf,pmaxk,pgrid,piezo,pressh,presshpz,dpdz,&
                             zpmin,delzp,maxpk,thetaz,dthetazdz,iwater
        
        IMPLICIT NONE
        
        INTEGER, INTENT(out) :: pcount,pnum
        INTEGER :: i,j,ilast,jlast,coords,ios,invert,pflag
        INTEGER :: maxi,maxj,k,error,ijorder,ktop,headlong
        INTEGER :: kcount,headlinecount,last,l
        
        REAL(pr) :: phead,vmc,delzpold
        REAL(pr) :: x,y,z,xlast,ylast,zlast,dzp,zlay(nmat)
        REAL(pr) :: tsat,tres,a,n,mswcc,rswcc,cfx
        REAL(pr), ALLOCATABLE :: tempz(:,:,:)
        REAL(pr), ALLOCATABLE :: tempp(:,:,:),tempthetaz(:,:,:)   
        
        CHARACTER*60 :: heading 
        CHARACTER*220, INTENT(in) :: prfile
        CHARACTER*3, INTENT(out) :: prcoords
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'reading 3D pressure-head file'          
        zpmin = 999999.0_pr
        pgrid = 1
        coords = 0
        pzsurf = 0
        pmaxk = 0
        ilast = 0
        jlast = 0
        headlinecount = 0
        pcount = 0  
        zlast = rnull
        x = rnull
        y = rnull
        z = rnull
        maxi = 0
        maxj = 0  
        k = 1
        ios = 1
        delzp=0.0_pr
        ijorder = -1
        invert = -1
        pnum = 0
        error = 0
        headlong = 0

        OPEN (13,STATUS = 'old',FILE = prfile,IOSTAT=ios)
        IF (ios.ne.0) THEN
          errmessage = 'error opening pressure-head file:'
          solution = 'check for existence of file with specified name and location'
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'ch ',0,prfile)                   
        END IF

! Read comment lines or blank lines at beginning of file
        DO 
          headlinecount = headlinecount + 1
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT
          headlong = LEN_TRIM(heading)
          IF (heading(1:1).ne.'#'.and.headlong.gt.0) EXIT 
        END DO
        
        IF (ios.eq.0) THEN 
          DO
! If coordinate system header line is not after comments, then exit           
            IF (heading(1:4).ne.'coor'.and.heading(1:4).ne.'COOR'&
                .and.heading(1:4).ne.'Coor') EXIT
! Read in coordinate system, z spacing format and maximum number of z values   
        
            READ (13,*,IOSTAT=ios) prcoords 
            IF (ios.ne.0) EXIT        
            headlinecount = headlinecount + 1                   
            SELECT CASE (prcoords)
              CASE ('ijk','IJK')
                coords = 1
                prcoords = 'ijk'
              CASE ('ijz','IJZ')              
                coords = 2
                prcoords = 'ijz'                
              CASE ('xyz','XYZ')               
                coords = 3
                prcoords = 'xyz'  
              CASE DEFAULT
                EXIT                  
            END SELECT
 ! If coordinate system is ijk, read in vertical spacing
            IF (coords.eq.1) THEN
              READ (13,1000,IOSTAT=ios) heading
              headlinecount = headlinecount+1
              IF (ios.ne.0.or.(heading(1:4).ne.'delz'.and.heading(1:4).ne.'DELZ')&
                 .and.heading(1:4).ne.'Delz') EXIT 
              READ (13,*,IOSTAT=ios) delzp,zpmin
              headlinecount= headlinecount+1                                                 
            END IF                                        
            EXIT
          END DO    
        END IF   

        solution = '3D pressure-head file must be in 3D file format, as described in Scoops3D manual'
        IF (coords.eq.0) THEN
          errmessage = 'error reading 3D pressure-head file -- header lines must specify coordinates: ijk,ijz or xyz'              
          Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')         
          error = 1
        END IF  

        IF (coords.eq.1.and.delzp.eq.0) THEN
          errmessage = 'error reading 3D pressure-head file -- "delz" is required with ijk coordinates'              
          Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                                      
          error = 1
        END IF                    

        IF (coords.eq.1.and.zpmin.eq.999999.0_pr) THEN   
          errmessage = 'error reading 3D pressure-head file -- "zmin" is required with ijk coordinates'              
          Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                                            
          error = 1
        END IF 
         
        IF (error.eq.1.or.ios.ne.0) THEN
          errmessage = 'error reading 3D pressure-head file'              
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')                 
        END IF        
          
!    Read pressure head data from file first time to determine spacing, and other parameters.
        DO
! read xyz or ijk - the distinction is not important here, so the variables x,y,z are used
          READ (13,*,IOSTAT=ios)x,y,z
          IF (ios.ne.0) EXIT
! If coordinate system is not ijk, determine vertical spacing, minimum z, etc.        
          IF (coords.ne.1) THEN
            delzp = z - zlast
! if this is not the first line, check to see if spacing is irregular    
            IF (zlast.ne.rnull) THEN      
              IF ((xlast.eq.x.and.ylast.eq.y).and.delzp.ne.delzpold) pgrid = 2
            END IF
            zlast = z
            delzpold = delzp
            IF (z.lt.zpmin) zpmin = z
          END IF  
          IF (x.eq.xlast.and.y.eq.ylast) THEN
            kcount = kcount+1
            IF (pmaxk.lt.kcount) pmaxk = kcount
          ELSE 
            kcount = 1
            xlast = x
            ylast = y
          END IF   
        END DO     

        IF (ios.ne.0.and.(x.eq.rnull.or.y.eq.rnull.or.z.eq.rnull)) THEN
          errmessage = 'error reading 3D pressure-head file'              
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
        END IF
! rewind file to read a second time        
        REWIND (13)

! skip through header lines when reading the file the second time      
        DO i=1,headlinecount
          READ (13,1000,IOSTAT=ios) heading
          IF (ios.ne.0) EXIT
        END DO       

        IF (ios.ne.0) THEN
          errmessage = 'error reading 3D pressure-head file'              
          CLOSE (13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
        END IF            

        solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
!    Allocate pressure arrays to be used throughout program.          
        ALLOCATE (pressh(nx,ny,pmaxk),dpdz(nx,ny,pmaxk),maxpk(nx,ny),STAT=error)
        ALLOCATE (presshpz(nx,ny,pmaxk),STAT=error)
!    If variably saturated conditions, allocate thetaz arrays.
        IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) &
           ALLOCATE (thetaz(nx,ny,pmaxk),dthetazdz(nx,ny,pmaxk),STAT=error)
        IF (error.ne.0) THEN
          errmessage = 'pressure head arrays not allocated successfully' 
          CLOSE(13)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
        END IF          
        pressh = 0.0_pr
        dpdz = 0.0_pr
        maxpk = 0
        presshpz = rnull 
        IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) THEN
          thetaz = 0.0_pr 
          dthetazdz = 0.0_pr
        END IF       
 
!    Read pressure head data from file.
        DO          
          IF (coords.eq.1) THEN
            IF (iwater.eq.4) THEN
              READ (13,*,IOSTAT=ios) i,j,k,phead,vmc
            ELSE
              READ (13,*,IOSTAT=ios) i,j,k,phead
            END IF           
            IF (ios.ne.0) EXIT
            z = REAL(k-1,pr)*delzp + zpmin
          ELSE 
            IF (coords.eq.2) THEN
              IF (iwater.eq.4) THEN
                READ (13,*,IOSTAT=ios) i,j,z,phead,vmc
              ELSE
                READ (13,*,IOSTAT=ios) i,j,z,phead
              END IF        
              IF (ios.ne.0) EXIT             
            ELSE
              IF (iwater.eq.4) THEN
                READ (13,*,IOSTAT=ios) x,y,z,phead,vmc 
              ELSE
                READ (13,*,IOSTAT=ios) x,y,z,phead
              END IF              
!     Assumes x and y given for DEM cell center
              IF (ios.ne.0) EXIT
              i = INT((x-xll)/delxy) + 1
              j = INT((y-yll)/delxy) + 1
            END IF
            IF (pgrid.eq.1)   k = NINT((z-zpmin)/delzp) + 1
          END IF 
  
          IF (i.lt.1.or.j.lt.1.or.k.lt.1) THEN
            solution = 'check coordinates of 3D pressure-head file'          
            errmessage = 'Scoops3D found negative values for i,j,or k in 3D pressure-head file' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')          
            PRINT *, '         ijk = ',i,j,k
            WRITE (39,*) '        ijk = ',i,j,k              
            CLOSE(13)
            CLOSE (20)
            CLOSE(39)
            STOP
          END IF
   
!     If irregular grid, determine ordering of data, i.e. whether all z are given for
!     each i and j, or all i,j for each z.      
          IF (pgrid.eq.2) THEN
            IF (ilast.ne.0.and.jlast.ne.0) THEN ! If not first line of input
              IF (ijorder.eq.-1) THEN 
!    Compare first and second line of input to determine order of data.               
                IF (i.eq.ilast.and.j.eq.jlast) THEN
!    All z values are entered at once for each i and j.
                  ijorder = 1
                  IF (z.lt.zlast) THEN
!    Data must be inverted so min z is associated with min k.
                    invert = 1 
                  ELSE
                    invert = 0
                  END IF
                ELSE
!    All i and j are entered for each z value.
                  ijorder = 0
                END IF
              END IF
!    Determine when to increment k or restart at k=1 at new i,j cell.              
              IF (ijorder.eq.1) THEN            
                IF (i.eq.ilast.and.j.eq.jlast) THEN
                  k = k + 1
                ELSE
                  k = 1
                END IF
              ELSE
                IF (z.ne.zlast) THEN
                  k = k + 1
                  IF (invert.eq.-1) THEN
                    IF (z.lt.zlast) THEN
                      invert = 1
                    ELSE
                      invert = 0
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF

          solution = 'check coordinates of 3D pressure-head file -- horizontal dimensions should agree with DEM'
          IF (i.gt.nx) THEN          
            errmessage = 'number of cells in x-direction for 3D pressure-head file does not agree with DEM' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')           
            PRINT *,'maximum value for i, in 3D pressure-head file = ',i,'ncols in DEM  = ',nx
            WRITE (39,*) 'maximum value for i, in 3D pressure-head file = ',i,'ncols in DEM  = ',nx          
            CLOSE (20)
            CLOSE(13)
            CLOSE(39)
            STOP
          END IF
          IF (j.gt.ny) THEN
            errmessage = 'number of cells in y-direction for 3D pressure-head file does not agree with DEM' 
            Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')           
            PRINT *,'maximum value for j, in 3D pressure-head file = ',j,'ncols in DEM  = ',ny
            WRITE (39,*) 'maximum value for j, in 3D pressure-head file = ',j,'ncols in DEM  = ',ny          
            CLOSE (20)
            CLOSE(13)
            CLOSE(39)
            STOP
          END IF    
                 
          pressh(i,j,k) = phead
          IF (iwater.eq.4) thetaz(i,j,k) = vmc
          presshpz(i,j,k) = z 

!    Keep track of number of k values at each i,j.          
          maxpk(i,j) = maxpk(i,j) + 1

          IF (maxi.lt.i) maxi = i
          IF (maxj.lt.j) maxj = j
               
          ilast = i
          jlast = j
          zlast = z     

        END DO
           
        solution = 'reduce memory requirements. See "Practical Considerations" chapter of Scoops3D manual'
!    If grid has low z associated with high k, invert k element of arrays.       
        IF (invert.eq.1) THEN
          ALLOCATE (tempp(nx,ny,pmaxk),tempz(nx,ny,pmaxk),STAT=error)
          IF (iwater.eq.4)  ALLOCATE (tempthetaz(nx,ny,pmaxk),STAT=error)          
          IF (error.ne.0) THEN
            errmessage = 'temporary pressure head arrays not allocated successfully' 
            CLOSE(13)
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')            
          END IF
          DO j=1,ny
            DO i=1,nx
              IF (zdem(i,j).eq.rnull) CYCLE              
              DO k=1,maxpk(i,j)             
                tempp(i,j,k) = pressh(i,j,maxpk(i,j)-(k-1))
                tempz(i,j,k) = presshpz(i,j,maxpk(i,j)-(k-1))
                IF (iwater.eq.4) tempthetaz(i,j,k) = thetaz(i,j,maxpk(i,j)-(k-1))
              END DO
            END DO
          END DO
        
          pressh = tempp
          presshpz = tempz
          DEALLOCATE (tempp) 
          DEALLOCATE (tempz)
          IF (iwater.eq.4) THEN
            thetaz = tempthetaz
            DEALLOCATE (tempthetaz)
          END IF  
        END IF     

!     Calculate moisture content each pressure head location, if using a SWCC
        IF (iwater.eq.5.or.iwater.eq.6) THEN
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
            ktop = maxpk(i,j) 
            last = 1
            tsat = thetasat(1)
            tres = thetares(1)                         
            IF (iwater.eq.5) THEN  ! vanGunuchten curve fit  
                a = vga(1)
                n = vgn(1)                      
            END IF
            IF (iwater.eq.6) THEN ! Fredlund and Xing curve fit 
                a = fxa(1)
                n = fxn(1)
                mswcc = fxm(1)
                rswcc = fxr(1)                  
            END IF                                     
            DO k = ktop,1,-1
!    Find layer at current depth to determine saturated and residual moisture contents              
              DO l = last,nmat
                IF (zlay(l).ne.rnull) THEN
                  tsat = thetasat(l)
                  tres = thetares(l)                  
                  IF (iwater.eq.5) THEN  ! vanGunuchten curve fit  
                      a = vga(l)
                      n = vgn(l)                      
                  END IF
                  IF (iwater.eq.6) THEN ! Fredlund and Xing curve fit 
                      a = fxa(l)
                      n = fxn(l)
                      mswcc = fxm(l)
                      rswcc = fxr(l)                  
                  END IF    
!   If the correct layer has been found, calculate theta                  
                  IF (presshpz(i,j,k).ge.zlay(l)) THEN  ! If 3D node below this layer 
                    IF (iwater.eq.5) THEN  ! vanGunuchten curve fit  
                      thetaz(i,j,k) = tres + (tsat-tres)/(1+(a*abs(pressh(i,j,k))*gamw)**n)**(1-1/n)
                    END IF
                    IF (iwater.eq.6) THEN ! Fredlund and Xing curve fit 
                        cfx = 1 - log(1+abs(pressh(i,j,k))*gamw/rswcc)/log(1+(1000000/rswcc))
        		thetaz(i,j,k) = cfx * tsat / (log(exp(1.) + (abs(pressh(i,j,k))*gamw/a)**n)**mswcc ) 
                    END IF                               
                    last = l
                    EXIT
                  END IF
                END IF
              END DO                                   
            END DO
          END DO
         END DO
        END IF   ! if SWCC (iwater = 5 or 6) 
               
!    Loop through DEM grid to find pressure gradients and piezometric surface.
        DO j = ny,1,-1
          DO i = nx,1,-1
            IF (zdem(i,j).eq.rnull) CYCLE
            IF (maxpk(i,j).eq.0) THEN
              errmessage = 'missing pressure-head data for column ' 
              Call WriteError(2,errmessage,problemtype,'no','no ',0,' ')                
              PRINT *,'       i = ',i,' j = ',j
              WRITE (39,*) '       i = ',i,' j = ',j 
              CLOSE(13)
              CLOSE (39)
              CLOSE (20)
              STOP
            END IF  
            pflag = 0 
            ktop = maxpk(i,j)                                              
            DO k = maxpk(i,j)-1,1,-1
!     If pressure is given for z value above DEM, set equal to DEM.
!     Change the maximum k value to refer to the first one to fall above the DEM.
              IF (presshpz(i,j,k+1).ge.zdem(i,j)) THEN
                IF (presshpz(i,j,k).ge.zdem(i,j)) CYCLE                
                presshpz(i,j,k+1) = zdem(i,j) 
                ktop = k + 1  
              END IF          
              dzp = presshpz(i,j,k+1) - presshpz(i,j,k)
  
!     Find pressure gradient at each depth. Note: 0 gradient for highest depth.
              IF (dzp.eq.0.0_pr) THEN
                dpdz(i,j,k) = 0.0_pr
                IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) dthetazdz(i,j,k) = 0.0_pr
              ELSE
                dpdz(i,j,k) = (pressh(i,j,k+1) - pressh(i,j,k))/dzp
                IF (iwater.eq.4.or.iwater.eq.5.or.iwater.eq.6) dthetazdz(i,j,k) = (thetaz(i,j,k+1) - thetaz(i,j,k))/dzp
              END IF

!     Search for change in sign of pressure head from negative to positive
!     with increasing depth to determine location of piezometric surface. 
!     The highest piezometric surface found at each cell will be kept.  
              IF (pressh(i,j,k).ge.0.0_pr.and.pflag.eq.0) THEN
                pflag = 1
!     If piezometric surface found allocate piezo array and set flag pzsurf.
                IF (pzsurf.eq.0) THEN
                  IF (.NOT.ALLOCATED(piezo)) ALLOCATE (piezo(nx,ny))
                  piezo=rnull
                  pzsurf=1
                END IF  
!     If at k=ktop-1, check whether highest data node (ktop) also has
!     positive pressure and assign piezo surface value to elevation of
!     ktop if so.
                IF (k.eq.ktop-1.and.pressh(i,j,ktop).ge.0.0_pr) THEN
                  piezo(i,j) = presshpz(i,j,ktop) 
                ELSE         
!     Use linear interpolation for piezometric surface location.
                  IF (dpdz(i,j,k).ne.0.0_pr) THEN
                    dzp = -pressh(i,j,k)/dpdz(i,j,k)  
                  ELSE
                    dzp = 0.0_pr
                  END IF  
                  piezo(i,j) = presshpz(i,j,k) + dzp   
                END IF
!    Count number of grid cells where piezometric surface found
                pcount = pcount + 1
              END IF                 
            END DO  !  Loop on k
!    Set maxpk array equal to highest k given if all z fall below DEM, 
!    or k associated with first node to fall above DEM.
            maxpk(i,j) = ktop
          END DO
        END DO        

!   Calculate number of cells with pressure data read.        
        pnum = maxi*maxj      
        
        CLOSE(13)

        RETURN        
        
1000    FORMAT(A)

     END SUBROUTINE readpressh
   
     

