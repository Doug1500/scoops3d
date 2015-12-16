        SUBROUTINE readstrat(sfilin)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!      This subroutine reads and processes the ASCII array data files that 
!      contain stratigraphic information.  There must be a file for each 
!      layer above the bottom layer, containing the elevation of the lower 
!      boundary of the layer at each DEM cell. 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called in Readin
!
!       VARIABLES
!
!        activelay(nmat-1) -- number of active (nonnull) columns in each layer
!        delxy -- DEM grid spacing
!        delxy2 -- resolution read for strat file
!        dlay -- layer thickness in column
!        dlayd -- thickness of dry layer in each column
!        dlayw -- thickness of wet layer in each column
!        duwtlaydz(nx,ny,nmat-1) -- gradient of depth-weighted unit weights between layers
!        ext -- layer number, extension for strat files
!        gamr(nmat,3) -- array of total, partially saturated, and saturated
!           unit weights for each layer
!        i,j -- DEM grid array location
!        ios -- status of read calls
!        iwater -- flag indicating method for modeling water pressure.
!           0=no water pressures, 1=ru approximation, 2= piezometric surface, 
!           3=input 3-d pressure file
!        layer(nmat,nx,ny) -- array of bottom elevations for material layers
!        ll -- layer number
!        llast -- keeps track of layer number of last layer found in  
!           each column
!        maxlayer(nmat-1) --maximum value of each stratigraphic layer file
!        minlayer(nmat-1) -- min value of each stratigraphic layer file
!        moist -- determines which element of soil unit weight array to use: 1-dry weight
!           2-partially saturated, 3-fully saturated.
!        mmoist -- local value for moist
!        nmat -- number of material layers
!        nodata2 -- no data value for strat file
!        nullsame -- logical used to determine if same null value as DEM
!        numdat -- number of data points in each layer file
!        rnull -- real null value used throughout program 
!        nx -- number of DEM cells in x direction
!        nx2 -- number of x cells in strat file
!        ny -- number of DEM cells in y direction
!        ny2 -- number of y cells in strat file
!        piezo(nx,ny) -- array of piezometric surface at each DEM cell
!        pzsurf -- flag for whether piezometric surface was found,1=yes,0=no
!          readline -- counter for number of header lines read
!        sfilin -- root name for strat files
!        strtfile -- strat file name with extension
!        uwtlay(nx,ny,nmat-1) -- depth-weighted average unit weight at each layer
!           bottom with correct moist for water conditions
!        uwtsum -- sum of unit weights weighted for depth
!        vold -- dry volume of layer in column
!        volw -- wet volume of layer in column
!        stuff -- character variable for reading input
!        xll2,yll2 -- x and y origin of strat grid, read in header lines
!        xll,yll -- x and y origin of DEM grid, read in header lines
!        zdem(nx,ny) -- DEM elevations
!        zlay(nmat) -- elevation of layer in each column 
!        zllast -- elevation of last layer found in column
!
!      INPUT FILES
!        unit filename
!         14 'inputfilename_#.asc' -- file of layer boundary elevations.
!               Each layer should have own file with suffix giving the 
!               layer number, i.e. 'layer_1.asc'
!
!      OUTPUT FILES
!        unit filename
!         20 'inputfilename_out.txt' -- echo of input and results containing
!              overall minimum slip surface data, written in subroutines
!              Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!              Readstrength, and Writeout.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: nx,ny,delxy,zdem,xll,yll
        USE WaterData, ONLY: iwater,pzsurf,piezo,moist
        USE MaterialData, ONLY: nmat,layer,gamr,minlayer,maxlayer,activelay,&
                          gsameeach,uwtlay,duwtlaydz
        
        IMPLICIT NONE
        
        INTEGER :: i,ll,j,nx2,ny2,ios,readline,numdat,laywarn,mmoist,llast      
   
        REAL(pr) :: dlay,zlay(nmat),zllast,dlayd,dlayw
        REAL(pr) :: uwtsum,depth
        REAL(pr) :: delxy2,xll2,yll2,nodata2

        LOGICAL :: nullsame

        CHARACTER*220, INTENT(in) :: sfilin
        CHARACTER*7 :: ext
        CHARACTER*12 :: stuff
        CHARACTER*60 :: heading
        CHARACTER*220 :: strtfile
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'reading layer files'  
        nullsame = .TRUE.
        ios = 1 
        
!    Append file name number extension for each layer file.
        DO ll = 1,nmat-1 
          laywarn = 0
          IF (ll.le.9) THEN
            ext = '_'//CHAR(48+ll)//'.asc'
          ELSE
          IF (ll.le.19) THEN
            ext = '_1'//CHAR(48+(ll-10))//'.asc'
          ELSE
            IF (ll.le.29) ext = '_2'//CHAR(48+(ll-20))//'.asc'
            END IF
          END IF    
          strtfile = TRIM(sfilin)//ext
          OPEN (14,STATUS = 'old',FILE = strtfile,IOSTAT=ios)
          IF (ios.ne.0) THEN
            errmessage = 'error opening material layer file:'
            solution = 'check for existence of file with specified name and location'
            Call WriteError(1,errmessage,problemtype,solution,'ch ',0,strtfile)           
          END IF

          nx2 = 0
          ny2 = 0
          delxy2 = 0.0_pr
          nodata2 = 0.0_pr
          xll2 = 0.0_pr
          yll2 = 0.0_pr
          readline = 0

!    Read standard 6 line header lines for all array files.            
          DO
            READ (14,1000,IOSTAT=ios) heading
            IF (ios.ne.0.or.readline.eq.6) EXIT
            
            solution = 'material layer files must be in Esri ASCII raster format, as described in Scoops3D manual'
            SELECT CASE (heading(1:5))
              CASE ('ncols','NCOLS')   
                BACKSPACE 14     
                READ (14,*,IOSTAT=ios) stuff,nx2
                readline = readline + 1
                IF (ios.ne.0) THEN
                  errmessage = 'error reading "ncols" header line of material layer file '              
                  CLOSE(14)
                  Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')                   
                END IF

              CASE ('nrows','NROWS')   
                BACKSPACE 14     
                READ (14,*,IOSTAT=ios) stuff,ny2
                readline = readline + 1
                IF (ios.ne.0) THEN
                  errmessage = 'error reading "nrows" header line of material layer file '              
                  CLOSE(14)
                  Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')                   
                END IF

              CASE ('xllco','XLLCO')   
                BACKSPACE 14     
                READ (14,*,IOSTAT=ios) stuff,xll2
                readline = readline + 1
                IF (ios.ne.0) THEN
                  errmessage = 'error reading "xllcorner" header line of material layer file '              
                  CLOSE(14)
                  Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')                  
                END IF
                           
              CASE ('yllco','YLLCO')   
                BACKSPACE 14     
                READ (14,*,IOSTAT=ios) stuff,yll2
                readline = readline + 1
                IF (ios.ne.0) THEN
                  errmessage = 'error reading "yllcorner" header line of material layer file '              
                  CLOSE(14)
                  Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')                  
                END IF       

              CASE ('cells','CELLS')   
                BACKSPACE 14     
                READ (14,*,IOSTAT=ios) stuff,delxy2
                readline = readline + 1
                IF (ios.ne.0) THEN
                  errmessage = 'error reading "cellsize" header line of material layer file '              
                  CLOSE(14)
                  Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')                  

                END IF
            
              CASE ('NODAT','nodat')   
                BACKSPACE 14     
                READ (14,*,IOSTAT=ios) stuff,nodata2
                readline = readline + 1
                IF (ios.ne.0) THEN
                  errmessage = 'error reading "nodata" header line of material layer file '              
                  CLOSE(14)
                  Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')
                END IF
                EXIT

            END SELECT
          END DO   
        
          IF (readline.ne.6) THEN
            errmessage = 'error reading material layer file -- missing or incorrect header data for layer #'              
            CLOSE(14)
            Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')            
          END IF

          solution = 'check extent and cellsize of input grid files'
!    Check that array size, cell size, and x and y corner locations are equivalent to DEM array.            
          IF ((nx.ne.nx2).or.(ny.ne.ny2)) THEN
            errmessage = 'DEM file and material layer file must have equivalent ncols and nrows, layer #'              
            CLOSE(14)
            Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')           
          END IF
        
          IF (delxy2.ne.delxy) THEN
            errmessage = 'DEM file and material layer file must have equivalent cellsize, layer #'              
            CLOSE(14)
            Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')             
          END IF
        
          IF (xll2.ne.xll.or.yll2.ne.yll) THEN
            errmessage = 'DEM file and material layer file must have equivalent xllcorner and yllcorner, layer #'              
            CLOSE(14)
            Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')             
          END IF
          
          IF (nodata2.ne.rnull) nullsame = .false.
          
          numdat = 0
 
          solution = 'check format of input grid data; see "Program Input" chapter of Scoops3D manual'        
!     Read the bottom elevations of the layer.          
          DO j = ny,1,-1
            READ (14,*,IOSTAT=ios) (layer(ll,i,j),i=1,nx)
            IF (ios.ne.0.and.i.le.nx) THEN
              errmessage = 'error reading material layer file # '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ')            
            END IF
!     Find the highest and lowest elevation of each layer and the number of DEM cells
!     intersected by the layer.        
            DO i=1,nx
              IF (layer(ll,i,j).ne.nodata2) THEN
                IF (layer(ll,i,j).lt.minlayer(ll)) minlayer(ll) = layer(ll,i,j)
                IF (layer(ll,i,j).gt.maxlayer(ll)) maxlayer(ll) = layer(ll,i,j)
                activelay(ll) = activelay(ll) + 1
              END IF
            END DO
            numdat = numdat + nx   
          END DO          

         solution = 'check dimensions of input data, number of data records should equal ncols*nrows'
!    Check that layer data are given for every DEM cell.          
          IF (numdat.ne.(nx*ny)) THEN
            errmessage = 'layer file dimensions are incorrect'              
            CLOSE(14)
            Call WriteError(1,errmessage,problemtype,solution,'int',ll,' ') 
          END IF  
            
!    Reset nodata2 values to standard program null value.          
          IF (.not.nullsame) WHERE (layer.eq.nodata2) layer = rnull          
        
          CLOSE (14)
        
        END DO  ! loop on number of layers
        

!     Find thicknesses of each layer in each column with a valid DEM elevation
        DO j = 1,ny
          DO i = 1,nx 
            IF (zdem(i,j).ne.rnull) THEN
              zllast = zdem(i,j)
              llast = 0
              uwtsum = 0.0_pr
              mmoist = moist
              DO ll = 1,nmat-1                     
                zlay(ll) = layer(ll,i,j)
                IF (zlay(ll).ne.rnull) THEN
!     Find lowest layer that falls above DEM and set higher layer to rnull
                  IF (llast.gt.1.and.ll.gt.1.and.zlay(ll).ge.zdem(i,j).and.zllast.eq.zdem(i,j)) THEN
                    layer(llast,i,j) = rnull
                    layer(ll,i,j) = rnull
                    llast = ll
                    CYCLE 
                  END IF              
!     If layer is not below prior layer, change layer elevation to rnull.
                  IF (zllast.le.zlay(ll)) THEN
                    IF (laywarn.eq.0) THEN
                      laywarn = 1
                      solution = 'check elevations specified in material layer file'
                      errmessage = 'previous material layer or DEM elevation found at or below layer #'              
                      Call WriteError(0,errmessage,problemtype,solution,'int',ll,' ')                     
                      PRINT *,'Scoops3D will set layer thickness to 0 in cells where elevation is below previous layer' 
                      WRITE (39,*)'Scoops3D will set layer thickness to 0 in cells where elevation is below previous layer'
                    END IF
                    layer(ll,i,j) = rnull
                    CYCLE
                  END IF
                  dlay = zllast - zlay(ll)
                  depth = zdem(i,j) - zlay(ll)
!     If no piezometric surface or  partial sat and sat unit weights are same,
!     find depth-weighted unit weights in column.
                  IF ((pzsurf.eq.0).or.(gsameeach.eq.1)) THEN
                    uwtsum = uwtsum + dlay*gamr(ll,mmoist)
!    If piezometric surface determine when to use partially saturated or 
!    saturated unit weights based on location of water table. 
                  ELSE
                    IF (pzsurf.eq.1) THEN    
!    If piezometric surface is above last layer or below current layer                
                      IF ((mmoist.eq.3).or.(piezo(i,j).lt.zlay(ll))) THEN
                        uwtsum = uwtsum + dlay*gamr(ll,mmoist)
                      ELSE  !  piezometric surface falls within layer
                        mmoist = 3
                        dlayd = zllast - piezo(i,j)
                        dlayw = piezo(i,j) - zlay(ll)
                        uwtsum = uwtsum + dlayd*gamr(ll,2)+dlayw*gamr(ll,3)
                      END IF
                    END IF                               
                  END IF
                  uwtlay(i,j,ll) = uwtsum/depth
                  IF (ll.gt.1.and.llast.ne.0)&
                    duwtlaydz(i,j,ll) = (uwtlay(i,j,llast)-uwtlay(i,j,ll))/dlay
                  zllast = zlay(ll)
                  llast = ll
                END IF
              END DO  ! loop on layers
            END IF  
          END DO  ! loop on i
        END DO  ! loop on j

        RETURN

1000    FORMAT(A)

        END
