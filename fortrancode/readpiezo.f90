        SUBROUTINE Readpiezo(pzfile)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine reads and processes the ASCII data file that 
!       contains piezometric information.  The data file should contain
!       the elevation of the piezometric surface at each DEM cell.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called by Readin
!
!      VARIABLES
!
!      delxy -- DEM grid resolution (delta x, delta y)
!      delxy2 -- resolution read for piezometric file
!      heading -- character variable for reading input
!      i,j -- DEM grid array location
!      ios -- status of read calls
!      nodata2 -- no data value for Piezometric file
!      nullsame -- logical used to determine if same null value as DEM
!      nump -- number of piezometric data
!      nx -- number of DEM columns
!      nx2 -- number of columns in piezometric file
!      ny -- number of DEM rows
!      ny2 -- number of rows in piezometric file
!      pzfile -- piezometric file name
!      piezo(i,j) -- array of piezometric surface elevations for
!              each DEM column
!      pzsurf -- flag for whether piezometric surface was found,1=yes,0=no
!      readline -- counter for number of header lines read
!      rnull -- real null value used throughout program 
!      stuff -- character variable for reading input
!      xll2,yll2 -- x and y origin of piezometric grid, read in header lines
!
!      INPUT FILES READ
!       unit filename
!         14 'inputfilename.#' -- file of piezometric elevations at each
!                  DEM cell
!
!      OUTPUT FILES
!       unit filename
!         20 'inputfilename_out.txt' -- echo of input and results containing
!              overall minimum slip surface data, written in subroutines
!              Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!              Readstrength, and Writeout.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
        USE CommonData
        USE GridData, ONLY: nx,ny,delxy,xll,yll
        USE WaterData, ONLY: piezo,pzsurf
            
        IMPLICIT NONE 
      
        INTEGER :: i,j,nx2,ny2,ios,nump,readline
        
        REAL(pr) :: delxy2,xll2,yll2,nodata2

        CHARACTER*220, INTENT(in) :: pzfile
        CHARACTER*60 :: heading
        CHARACTER*12 :: stuff
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
        
        LOGICAL :: nullsame  

        errmessage = ' '
        solution = ' '
        problemtype = 'reading piezometric surface file'                
        ios = 0
        pzsurf = 1
        nump = 0
        readline = 0
        
        OPEN (14,STATUS = 'old',FILE = pzfile,IOSTAT=ios)
        IF (ios.ne.0) THEN
          errmessage = 'error opening piezometric surface file:'
          solution = 'check for existence of file with specified name and location'
          Call WriteError(1,errmessage,problemtype,solution,'ch ',0,pzfile)         
        END IF
	
        DO 
          READ (14,1000,IOSTAT=ios) heading
          IF (ios.ne.0.or.readline.eq.6) EXIT

          solution = 'piezometric surface file must be in Esri ASCII raster format, as described in Scoops3D manual'
          SELECT CASE (heading(1:5))
          CASE ('ncols','NCOLS')   
            BACKSPACE 14     
            READ (14,*,IOSTAT=ios) stuff,nx2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "ncols" header line of piezometric surface file '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
            END IF

          CASE ('nrows','NROWS')   
            BACKSPACE 14     
            READ (14,*,IOSTAT=ios) stuff,ny2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "nrows" header line of piezometric surface file '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF

          CASE ('xllco','XLLCO')   
            BACKSPACE 14     
            READ (14,*,IOSTAT=ios) stuff,xll2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "xllcorner" header line of piezometric surface file '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
                           
          CASE ('yllco','YLLCO')   
            BACKSPACE 14     
            READ (14,*,IOSTAT=ios) stuff,yll2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "yllcorner" header line of piezometric surface file '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF   

          CASE ('cells','CELLS')   
            BACKSPACE 14     
            READ (14,*,IOSTAT=ios) stuff,delxy2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "cellsize" header line of piezometric surface file '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
            
          CASE ('NODAT','nodat')   
            BACKSPACE 14     
            READ (14,*,IOSTAT=ios) stuff,nodata2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "nodata" header line of piezometric surface file '              
              CLOSE(14)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
            EXIT
                                 
          END SELECT
        END DO   
        
        IF (readline.ne.6) THEN
          errmessage = 'error reading piezometric surface file -- missing or incorrect header data'              
          CLOSE(14)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF

      solution = 'check extent and cellsize of input grid files'                
!     Compare DEM and piezometric grid dimensions as given in header info.
        IF ((nx.ne.nx2).or.(ny.ne.ny2)) THEN
          errmessage = 'DEM file and piezometric surface file must have equivalent ncols and nrows'              
          CLOSE(14)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF
        
        IF (delxy2.ne.delxy) THEN
          errmessage = 'DEM file and piezometric surface file must have equivalent cellsize'              
          CLOSE(14)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
        END IF
        
        IF (xll2.ne.xll.or.yll2.ne.yll) THEN
          errmessage = 'DEM file and piezometric surface file must have equivalent xllcorner and yllcorner'              
          CLOSE(14)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
        END IF
                           
        IF (nodata2.ne.rnull) nullsame = .false.

        solution = 'check format of input grid data; see "Program Input" chapter of Scoops3D manual'        
        DO j = ny,1,-1
          READ (14,*,IOSTAT=ios) (piezo(i,j),i=1,nx)
          IF (ios.ne.0.and.i.le.nx) THEN
            errmessage = 'error reading piezometric surface file on line number '              
            CLOSE(14)
            Call WriteError(1,errmessage,problemtype,solution,'int',ny-j+1,' ')            
          END IF
          nump = nump + nx   
        END DO
                                
        solution = 'check dimensions of input data, number of data records should equal ncols*nrows'               
        IF (nump.ne.(nx*ny)) THEN
          errmessage = 'piezometric file dimensions are incorrect'              
          CLOSE(14)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF  
            
        IF (.not.nullsame) WHERE (piezo.eq.nodata2) piezo = rnull               

        CLOSE (14)

        RETURN
        
1000    FORMAT (A)   

        END SUBROUTINE readpiezo
