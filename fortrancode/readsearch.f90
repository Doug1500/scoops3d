        SUBROUTINE Readsearch(sfile,nxsrch,nysrch,nsrchpt)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        This subroutine reads an ASCII data file that contains
!        3-d search grid information
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!        Called by Readin
!
!        VARIABLES
!
!        delxy -- DEM grid resolution (delta x, delta y)
!        delxy2 -- resolution read for search file
!        heading -- character variable for reading input
!        i,j -- DEM grid array location
!        ios -- status of read calls
!        ismin -- minimum i value of search grid relative to DEM grid
!        ismax -- maximim i value of search grid relative to DEM grid
!        jsmin -- minimum j value of search grid relative to DEM grid
!        jsmax -- maximim j value of search grid relative to DEM grid
!        nodata2 -- no data value for search file
!        nsrchpt -- total number of search grid points
!        numsrch -- number of values in file
!        nx -- number of DEM cells in x direction
!        nx2 -- number of x nodes in search file
!        ny -- number of DEM cells in y direction
!        ny2 -- number of y nodes in search file
!        readline -- counter for number of header lines read
!        rnull -- real null value used throughout program 
!        searchgrid - array to specify search grid points
!        sfile -- search file name
!        stuff -- character variable for reading input
!        xll2,yll2 -- x and y origin of search grid, read in header lines
!        xll,yll -- x and y origin of DEM grid, read in header lines
!        xlls,ylls -- calculated x and y origin based on ismin and jsmin
!
!        INPUT FILES
!         unit filename
!           10 'inputfilename.#' -- file that defines search grid
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
        USE SearchData, ONLY:ismin,ismax,jsmin,jsmax,searchgrid
        USE GridData, ONLY: nx,ny,xll,yll,yllcorner,delxy
        
        IMPLICIT NONE

        INTEGER, INTENT(out) :: nxsrch,nysrch,nsrchpt
        INTEGER :: ios,i,j,readline,error,numsrch
        
        REAL(pr) :: xll2,yll2,delxy2,nodata2,xlls,ylls

        CHARACTER*220, INTENT(in) :: sfile
        CHARACTER*60 :: heading
        CHARACTER*12 :: stuff 
        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'reading search file'                  
        ios = 1
        nsrchpt = 0
        readline = 0
        numsrch = 0

        OPEN (10,STATUS = 'old',FILE = sfile,IOSTAT=ios)
        IF (ios.ne.0) THEN
           errmessage = 'error opening search file:'
           solution = 'check for existence of file with specified name and location'
           Call WriteError(1,errmessage,problemtype,solution,'ch ',0,sfile)         
        END IF         

        DO
          READ (10,1000,IOSTAT=ios) heading
          IF (ios.ne.0.or.readline.eq.6) EXIT

          solution = 'search file must be in Esri ASCII raster format, as described in Scoops3D manual'
          SELECT CASE (heading(1:5))
          CASE ('ncols','NCOLS')   
            BACKSPACE 10     
            READ (10,*,IOSTAT=ios) stuff,nxsrch
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "ncols" header line of search file '              
              CLOSE(10)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')               
            END IF

          CASE ('nrows','NROWS')   
            BACKSPACE 10     
            READ (10,*,IOSTAT=ios) stuff,nysrch
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "nrows" header line of search file '              
              CLOSE(10)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
            END IF

          CASE ('xllco','XLLCO')   
            BACKSPACE 10     
            READ (10,*,IOSTAT=ios) stuff,xll2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "xllcorner" header line of search file '              
              CLOSE(10)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
                           
          CASE ('yllco','YLLCO')   
            BACKSPACE 10     
            READ (10,*,IOSTAT=ios) stuff,yll2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "yllcorner" header line of search file '              
              CLOSE(10)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF 

          CASE ('cells','CELLS')   
            BACKSPACE 10     
            READ (10,*,IOSTAT=ios) stuff,delxy2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "cellsize" header line of search file '              
              CLOSE(10)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
            
          CASE ('NODAT','nodat')   
            BACKSPACE 10     
            READ (10,*,IOSTAT=ios) stuff,nodata2
            readline = readline + 1
            IF (ios.ne.0) THEN
              errmessage = 'error reading "nodata" header line of search file '              
              CLOSE(10)
              Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')              
            END IF
            EXIT        
            
          END SELECT
        END DO   
        
        IF (readline.ne.6) THEN
          errmessage = 'error reading search file -- missing or incorrect header data'              
          CLOSE(10)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')          
        END IF
          
        solution = 'check cellsize of input grid files'        
!     Compare DEM and search grid cell size and other parameters as given in header lines.
        IF (delxy.ne.delxy2) THEN                   
          errmessage = 'DEM file and search file must have equivalent cellsize'              
          CLOSE(10)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF     
        
        xlls = xll+REAL(ismin-1,pr)*delxy
        ylls = yll+REAL(jsmin-1,pr)*delxy
        
        solution = 'check horizontal alignment of input grid files'
        IF (xll2.ne.xlls) THEN
          errmessage = 'xllcorner of search grid does not match expected location calculated using xllcorner of DEM and ismin:'              
          Call WriteError(0,errmessage,problemtype,solution,'int',xlls,' ')        
        END IF
        
        IF (yll2.ne.ylls) THEN
          errmessage = 'yllcorner of search grid does not match expected location calculated using yllcorner of DEM and jsmin:'              
          Call WriteError(0,errmessage,problemtype,solution,'int',ylls,' ')             
        END IF

!     Calculate maximum i and j search node based on dimensions of search file.        
        ismax = ismin + (nxsrch-1)
        jsmax = jsmin + (nysrch-1)  
                   
        solution = 'reduce memory requirements; see "Practical Considerations" chapter of Scoops3D manual' 
        ALLOCATE (searchgrid(ismin:ismax,jsmin:jsmax), STAT=error)
        IF (error.ne.0) THEN
          errmessage = 'searchgrid data array not allocated successfully'                   
          CLOSE(10)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
        END IF  

        searchgrid = inull        

        solution = 'check format of input grid data, see "Program Input" chapter of Scoops3D manual'        
        DO j = jsmax,jsmin,-1
          READ (10,*,IOSTAT=ios) (searchgrid(i,j),i=ismin,ismax)
          IF (ios.ne.0.and.i.le.ismax) THEN
            errmessage = 'error reading search file on line number '              
            CLOSE(10)
            Call WriteError(1,errmessage,problemtype,solution,'int',ny-j+1,' ')          
          END IF
          numsrch = numsrch + (ismax-ismin)+1
        END DO
        
        nsrchpt = COUNT(searchgrid.gt.0) 

        solution = 'check dimensions of input data, number of data records should equal ncols*nrows'               
        IF (numsrch.ne.(nxsrch*nysrch)) THEN
          errmessage = 'search file dimensions are incorrect'              
          CLOSE(10)
          Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')         
        END IF

        CLOSE(10)
        RETURN
     
1000    FORMAT (A)        

        END
