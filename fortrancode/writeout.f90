       SUBROUTINE Writeout(iretro)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!      This subroutine creates output files for Program Scoops
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!    Called by Scoops3D
!
!     VARIABLES
!
!      colfile -- flag for opening error file reporting cols < limcol
!      delxy -- DEM grid resolution (delta x, delta y)
!      delz -- z resolution of search grid
!      failsurfdepth(nx,ny) -- array of failure surface depths below
!              the DEM surface
!      filtcount(nx,ny) -- counter of number of times each column was filtered by
!         minma.
!      filter -- indicates whether solution filters are applied (malpha and fos)
!      foscut -- factor of safety cutoff for removing slip surface from DEM
!      foslocal -- flag for writing local factors of safety to file, only
!         available for single failure option
!      fsangle(nx,ny) -- array of fall angles of minimum scoop at each DEM column
!      fsangle2d(nx,ny) -- array of fall angles of minimum 2d slice at each DEM column
!      fsarclength(nx,ny) -- array of arclengths of minimum 2d slice at each DEM column
!      fsarea2d(nx,ny) -- array of areas of minimum 2d slice at each DEM column
!      fsmin(nx,ny) -- array of minimum FOS at each DEM column
!      fsmin2d(nx,ny) -- 2-d minimum factor of safety at each DEM column
!      fsmin2d3d(nx,ny) -- 3-d minimum factor of safety associated with min 2d slice
!         at each DEM column
!      fsmin3d2d(nx,ny) -- 2-d minimum factor of safety associated with min 3d scoop
!         at each DEM column
!      fsrad(nx,ny) -- array of radius of minimum scoop at each DEM column
!      fsrad2d(nx,ny) -- array of radius of minimum 2d slice at each DEM column
!      fsvol(nx,ny) -- volume of minimum FOS scoop at each DEM column
!      fswidth2d(nx,ny) -- array of widths of minimum 2d slice at each DEM column
!      fsx(nx,ny),fsy(nx,ny),fsz(nx,ny) -- search sphere center of minimum FOS 
!         scoop at each DEM column
!      fsx2d(nx,ny),fsy2d(nx,ny),fsz2d(nx,ny) -- search sphere center of minimum FOS 
!         2d slice at each DEM column
!      icrit(nx,ny) -- array of volume or area flags indicating proximity to limits
!          of criteria range for each DEM column
!      icritlattice -- flag for creating file critfoslattice_out.3D for minimum FOS of 
!          critical surfaces at each search  grid node
!      ifailsurf -- flag for whether failure surface file is used
!      ifos2d -- flag for calculating 2-d factors of safety
!      ilattice -- flag for creating file foslattice_out.3D for minimum FOS at each search
!         grid node
!	     iretro -- counter for number of retrogressions
!      isubsurf -- flag for creating 3-d file of FOS below DEM surface&
!          1=create file in ijk format, 2=create file in xyz format,
!          3=create in .3D format
!      isrchmin,isrchmax -- array bounds of searchout output file
!      jsrchmin,jsrchmax -- array bounds of searchout output file
!      mcol(nx,ny) -- number of columns in minimum FOS scoop at each DEM column
!      method -- B for Bishops Simplified, O for Ordinary (Fellenius) method for
!         calculation of factor of safety.
!      mfile -- flag indicating whether m-alpha file was opened
!      nnset -- maximum number of subsets
!      nsetmax -- flag for exceeding number of subsets allowed nnset
!      ntry -- total number of slip surfaces for which FOS is calculated
!      nxsrchout,nysrchout -- number of search nodes in x and y directions in output
!         search grid file, includes intersection of search range with DEM range.
!      nx -- number of DEM cells in x direction
!      ny -- number of DEM cells in y direction
!      nz -- number of nodes used for 3-D FOS values when isubsurf=1 or 2 or 3
!      oangle -- slip angle associated with overall minimum factor of safety
!      oangle2d -- slip angle associated with overall minimum 2d factor of safety
!      oarclength -- arc length of overall minimum 2d factor of safety
!      oarclength3d2d -- arc length of 2d surface associated with overall min 3d scoop
!      oarea -- area of overall minimum scoop
!      oarea2d -- area of overall minimum 2d slice
!      oarea2d3d -- area of 3d surface associated with overall min 2d slice
!      oarea3d2d -- area of 2d surface associated with overall min 3d scoop
!      ocol -- number of columns in overall minimum scoop
!      ocol2d -- number of columns in overall minimum 2-D scoop
!      ocol2d3d -- number of columns in 3d surface associated with min 2d slice
!      ocol3d2d -- number of columns in 2d surface associated with min 3d scoop
!      ofos -- overall minimum factor of safety
!      ofos2d -- 2-d FOS associated with overall minimum 3-d factor of safety
!      ofos2d3d -- FOS of 3d slice associated with overall minimum 2d FOS
!      ofos3d2d -- FOS of 2d slice associated with overall minimum 3d FOS
!      ofsx,ofsy,ofsz -- search sphere center of minimum FOS scoop
!      ofsx2d,ofsy2d,ofsz2d -- search sphere center of minimum 2d FOS scoop
!      ominma -- overall min absolute value of m-alpha of all iterations for a surface
!      omset -- set number of overall minimum scoop
!      orad -- radius of sphere associated with overall minimum factor of safety
!      orad2d -- radius of arc associated with overall minimum 2d factor of safety
!      osliparea -- slip surface area of overall min 3d scoop
!      osliparea2d3d -- slip surface area of 3d slice associated with overall minimum 2d FOS
!      outputdir -- optional directory path to place output files into
!      ovangle -- slip angle of largest volume scoop with FOS < cutoff
!      ovarea -- area of largest volume scoop with FOS < cutoff
!      ovmset -- set number of largest volume scoop with FOS < cutoff
!      ovfos -- factor of safety of largest volume scoop with FOS < cutoff
!      ovrad -- radius of largest volume scoop with FOS < cutoff
!      ovfsx,ovfsy,ovfsz -- search sphere center of largest volume scoop with 
!         FOS < cutoff
!      ovol, oarea - volume and area of minimum  FOS scoop
!      ovol2d3d -- volume of 3d surface associated with overall min 2d slice
!      ovsliparea -- slip surface area of largest volume scoop with 
!         FOS < cutoff
!      ovvol,ovarea -- volume and area of largest volume scoop with FOS < cutoff
!      ovwt -- weight of largest volume scoop with FOS < cutoff
!      owidth2d -- wisth of overall min 2d FOS
!      owidth3d2d -- width of 2d surface associated with overall min 3d scoop
!      owt -- weight of overall min 3d scoop
!      owt2d3d -- weight of 3d slice associated with overall minimum 2d FOS
!      rad -- radius of search sphere
!      remove -- A= all scoops < FOS cutoff will be removed from new
!         DEM.  L= only largest volume scoop < FOS cutoff will be 
!         removed. Other= the scoop with the minimum FOS will be removed.
!      rnull -- real null value used throughout program 
!      searchout -- output array of area covered by search grid and indicating
!         which grid points were searched and whether they fell inside or outside
!         the DEM.
!      single -- flag for calculating single slip surface
!      vacriterion -- indicates primary control for intial radius (v or a
!         for volume or area)
!      xcolcount -- number of critical surfaces with total columns < limcol
!      xcount -- number of surfaces eliminated due to m-alpha limits or nonconvergence
!      xll,yll -- x and y origin of DEM grid, read in header lines
!      xllcorner,yllcorner -- character version of these header lines
!      xmax -- maximum x in zfos array
!      ymax -- maximum y in zfos array
!      zbot -- elevation of bottom of 3-d fos grid array
!      zdem(nx,ny) -- DEM elevations
!      zfos(nx,ny,nz) -- 3-d array of minimum factor of safety below DEM surface
!      zmax -- maximum DEM elevation
!	     zzmax -- maximum z in zfos array
!
!
!     OUTPUT FILES
!      unit filename
!        17 'filtergrid_out.asc' -- file of number of times each DEM column was
!             filtered due to absmina.  Opened and written in
!             Writeout if filter=1 and mflag=1.
!        19 'fos2d_out.asc' -- file of minimum 2-D factors of safety at each
!             DEM column. Opened and written in Writeout if ifos2d=1.
!        20 'inputfilename_out.txt' -- echo of input and results containing
!             overall minimum slip surface data, written in subroutines
!             Readin, Readpiezo, Readpressh, Readsearch, Readstrat,
!             Readstrength, and Writeout.
!        21 'spheres_out.okc' -- file of sphere center coordinates, radius, 
!             fall angle, # of columns in scoop, vol, and factor of safety 
!             of the least stable scoop at each DEM column. Opened and written
!             in Writeout.
!        22 'critcheck_out.asc' -- file of volume checks.  If volume of minimum
!             FOS scoop is < vmin+tol shows 1, if > vmax-tol shows 2,
!             if not close to volume boundary shows 0. Opened and written
!             in Writeout if isqout = 1.
!        23 'subsurffos_out.txt' or .3D-- file of minimum FOS in 3-D array, written
!             if isubsurf = 1,2 or 3. Opened and written in subroutine writeout.
!        24 'newDEM_out.asc' -- file of new DEM with certain scoops removed 
!             depending on value of flag remove. Opened and written in Writeout.
!        26 'fos3d_out.asc' -- file of minimum factors of safety at each
!             DEM column. Opened and written in Writeout.
!        27 'fosvol_out.asc' or 'fosarea_out.asc' -- file of volumes or areas (depending  
!             on primary criteria) for min FOS at each DEM column. Opened and written 
!             in Writeout.
!        28 'searchgrid_out.asc' -- search grid file.  Record of which nodes were searched
!             and whether they fell within the DEM bounds or not. Opened
!             and written in Scoops if isqout = 1.
!        29 'numcols_out.asc' -- file of # of columns for min FOS at each 
!             DEM column. Opened and written in Writeout if isqout = 1.
!        35 'arcs2d_out.okc' -- file of sphere center coordinates, radius, fall angle,
!             and 2-D factor of safety of the least stable 2-D solution at each DEM column. 
!             Opened and written in Writeout if ifos2d=1.
!        36 'slope_out.asc' -- slope of DEM at each cell
!        37 'fos3drel_out.asc' -- file of relative minimum factors of safety at each
!             DEM column. Relative factor of safety is defined and factor of safety divided
!             by the global minimum factor of safety (ofos).  Opened and written in Writeout 
!             if irelfos=1.
!        38 'fos2drel_out.asc' -- file of relative minimum 2-D factors of safety at each
!             DEM column. 2-D relative factor of safety is defined and factor of safety divided
!             by the global minimum 2-D factor of safety (ofos2d).Opened and written in Writeout 
!             if ifos2d=1 and irelfos=1.
!        41 'failsurfslope_out.asc' -- file of defined failure surface slopes
!        42 'failsurfdepth_out.asc' -- file of defined failure surface depths
!        43 'ordfos3d_out.asc' -- file of minimum factors of safety computed by the
!              Ordinary method at each DEM column. Opened and written in Writeout.
!              if method = 'b'
!        44 'boundcheck_out.asc' -- file of check for boundary limitations in
!              the search lattice. Opened and written in Writeout if isqout = 1
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        USE CommonData
        USE GridData, ONLY: nx,ny,nz,delxy,delz,xllcorner,yllcorner,rad,zmax,&
                            xll,yll,zdem,demflag,lengthunits,halfdelxy,xcenrot,ycenrot,zcenrot
        USE Strengthdata, ONLY: ceeunits,uwtunits                    
        USE FOSData, ONLY: remove,irelfos,isubsurf,ofos,orad,oangle,ofsx,ofsy,ofsz,ovol,oarea,&
                     ocol,ovfos,ovrad,ovangle,ovfsx,ovfsy,ovfsz,ovvol,ovarea,foscut,foslocal,&
                     fsmin,felminfos,ntry,ifos2d,ofos2d,filtcount,icritlattice,ilattice,ominma,single,method,mfile,&
                     ocol2d,oarea2d,oarclength,ofsx2d,ofsy2d,ofsz2d,owidth2d,orad2d,oangle2d,ocol2d,&
                     oarea3d2d,oarea2d3d,ocol3d2d,ocol2d3d,owidth3d2d,oarclength3d2d,ovol2d3d,icrit,&
                     osliparea,owt,owt2d3d,osliparea2d3d,ovwt,ovsliparea,ofos3d2d,ofos2d3d,xcount,zbot,&
                     fsx,fsy,fsz
        USE SearchData, ONLY: isqout,vacriterion,zsmin,zsmax,zsrchres,searchout,srchfile,&
                     isrchmin,jsrchmin,isrchmax,jsrchmax,nxsrchout,nysrchout,iincflag,irotcen,&
                     nsrchres,ismin,jsmin,ismax,jsmax
        USE SetData, ONLY: nsetmax,colfile,xcolcount,nnset,surfslope
        USE FailSurfData, ONLY: ifailsurf,failsurfdepth
        USE ResultsStats

        IMPLICIT NONE

        INTEGER, INTENT(in) :: iretro
        INTEGER :: i,j,jj,cardinality,numvars,ios

        REAL(pr) :: xmax,ymax,zzmax,x,y,z,f,ilatticemin
        REAL(pr) ::  xsmin,xsmax,ysmin,ysmax,zsmaxtrue       

        CHARACTER*2 :: wtunits
        CHARACTER*8 :: date
        CHARACTER*10 :: time
        CHARACTER*10 :: cnum,header
        CHARACTER*1 :: charetro
        CHARACTER*4 :: fosmeth
        
        cnum = '0123456789'
        charetro = cnum(iretro+1:iretro+1) 
        IF (method.eq.'b') then 
             fosmeth = 'Bish'
        ELSE
             fosmeth = 'Ord '
        END IF         

        CALL ArrayStats
        WRITE (20,1000)
        IF (iretro.gt.0) WRITE (20,4300) iretro
        WRITE (20,1100)
        IF (outputdir.ne.'') WRITE (20,1150) outputdir
        WRITE (20,1300) filin(bfile:nfile)//'_out.txt'
        WRITE (20,1300) filin(bfile:nfile)//'_errors_out.txt'
                                       
        CALL OpenFileWrite(36,20,'slope','asc',charetro)              
        WRITE (20,4600) minslope* 57.2957795_pr,maxslope* 57.2957795_pr 
        IF (ifailsurf.eq.1) THEN
          CALL OpenFileWrite(41,20,'failsurfslope','asc',charetro)              
          WRITE (20,4600) minfailslope,maxfailslope
          CALL OpenFileWrite(42,20,'failsurfdepth','asc',charetro)              
          WRITE (20,4600) minfaildepth,maxfaildepth
        END IF
                                 
        CALL OpenFileWrite(26,20,'fos3d','asc',charetro)
        WRITE (20,4600) minfsmin,maxfsmin
!!!   Fellenius comparison change        
        IF (method.eq.'b') CALL OpenFileWrite(43,20,'ordfos3d','asc',charetro)
         WRITE (20,4600) minfelfsmin,maxfelfsmin
!!!

        IF (single.ne.1) THEN                         
          IF (vacriterion.eq.'v') THEN
            CALL OpenFileWrite(27,20,'fosvol','asc',charetro)
            WRITE (20,4700) minfsvol,maxfsvol
          ELSE
            CALL OpenFileWrite(27,20,'fosarea','asc',charetro)
            WRITE (20,4700) minfsarea,maxfsarea
          END IF
        END IF
        CALL OpenFileWrite(21,20,'spheres','okc',charetro)
        numvars = 11
        if (ifos2d.eq.1) numvars = numvars + 1
        if (method.eq.'b') numvars = numvars +1
        WRITE(21,*) numvars,numdatapts,numvars*numdatapts  

        WRITE (20,1200)
        IF (irelfos.eq.1) THEN
          CALL OpenFileWrite(37,20,'fos3drel','asc',charetro)
          WRITE (20,4600) minfsrelmin,maxfsrelmin
        END IF          
        IF (demflag.eq.1.and.ofos.lt.foscut) THEN
          CALL OpenFileWrite(24,20,'newDEM','asc',charetro)
          WRITE (20,4800) minzdem,maxzdem
        END IF
        IF (single.ne.1) THEN
          IF (isqout.eq.1) THEN  
            CALL OpenFileWrite(29,20,'numcols','asc',charetro)            
            WRITE (20,4900) minmcol,maxmcol
            CALL OpenFileWrite(22,20,'critcheck','asc',charetro)
             IF (single.ne.1.and.srchfile.ne.1) CALL OpenFileWrite(44,20,'boundcheck','asc',charetro)
          END IF
          IF (isubsurf.eq.1.or.isubsurf.eq.2) CALL OpenFileWrite(23,20,'subsurffos','txt',charetro) 
          IF (isubsurf.eq.3) CALL OpenFileWrite(23,20,'subsurffos','3D',charetro)                     
        END IF  

        IF (icritlattice.eq.1) REWIND(15,iostat=ios)          
        IF (icritlattice.eq.1) THEN
          CALL OpenFileWrite(45,20,'critfoslattice','3D',charetro) 
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (45,8510) lengthunits,lengthunits,lengthunits,fosmeth
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (45,8520) lengthunits,lengthunits,lengthunits,fosmeth
            ELSE
              WRITE (45,8500) fosmeth
            END IF 
          END IF   
          READ (15,*,IOSTAT=ios) header              
          DO
            READ (15,6320,IOSTAT=ios) x,y,z,f
            IF (ios.ne.0) EXIT 
            IF (f.lt.100) THEN                 
              ilatticemin = 10000
            ELSE
              ilatticemin = f
            END IF    
            IF (f.lt.100) THEN
              DO j = 1,ny
                DO i = 1,nx      
                  IF (abs(fsx(i,j)-(x-xll)).lt.halfdelxy.and.abs(fsy(i,j)-(y-yll)).lt.halfdelxy&
                      .and.abs(fsz(i,j)-z).lt.(zsrchres/2)) THEN
                    IF (fsmin(i,j).lt.ilatticemin)  ilatticemin = fsmin(i,j)
                  END IF
                END DO
              END DO 
            END IF  
            write (45,6320) x,y,z,ilatticemin         
          END DO     
        END IF                
      
        IF (ilattice.eq.1) THEN
                 WRITE (20,1300) filin(bfile:nfile)//'_foslattice_out.3D'                 
                 PRINT *,outputdir(1:LEN_TRIM(outputdir))//&
                    filin(bfile:nfile)//'_foslattice_out.3D'                            
        END IF              
        IF (foslocal.eq.1.and.ofos.ne.nullhi)&
          WRITE (20,1300) filin(bfile:nfile)//'_foslocal_out.txt'
        IF (mfile.eq.1) THEN 
          CALL OpenFileWrite(17,20,'filtergrid','asc',charetro)
          WRITE (20,4900) minfiltcount,maxfiltcount
          WRITE (20,1300) filin(bfile:nfile)//'_filter_out.txt'                        
        END IF
        IF (colfile.eq.1.and.single.ne.1) WRITE (20,1300) filin(bfile:nfile)//'_ncolerr_out.txt'
        IF (ifos2d.eq.1.and.single.ne.1.and.ofos.ne.nullhi) THEN
          CALL OpenFileWrite(19,20,'fos2d','asc',charetro)
          WRITE (20,4600) minfsmin2d,maxfsmin2d
          IF (irelfos.eq.1) THEN 
            CALL OpenFileWrite(38,20,'fos2drel','asc',charetro)               
            WRITE (20,4600) minfsrelmin2d,maxfsrelmin2d                                                      
          END IF  
          CALL OpenFileWrite(35,20,'arcs2d','okc',charetro)              
          WRITE(35,*) 13,numdatapts,13*numdatapts
        END IF   
        IF (isqout.eq.1.and.iretro.eq.0) THEN
          CALL OpenFileWrite(28,20,'searchgrid','asc',charetro)
          WRITE (28,5500) nxsrchout
          WRITE (28,5600) nysrchout
          WRITE (28,5750) xll+REAL(isrchmin-1,pr)*delxy
          WRITE (28,5760) yll+REAL(jsrchmin-1,pr)*delxy
          WRITE (28,5800) delxy
          WRITE (28,5950) inull
          DO j = jsrchmax,jsrchmin,-1
            WRITE (28,*) (searchout(i,j),i=isrchmin,isrchmax)
          END DO 
          CLOSE (28)  !  Only write searchgrid file once for all retrogressions.
        END IF          
        
        WRITE (20,1000)  
        WRITE (20,1400)             
        IF (single.ne.1) WRITE (20,1500) ntry
        IF (remove.ne.'N') THEN
          IF (ofos.lt.foscut) THEN
            WRITE (20,1600)
          ELSE
            WRITE (20,1700)
          END IF
        END IF        
        IF (mfile.eq.1) THEN
          WRITE (20,1900) xcount 
          WRITE (20,2000)
        END IF        
        IF (colfile.eq.1) THEN
          IF (single.ne.1) THEN
            WRITE (20,2100) xcolcount
            WRITE (20,2200)
          ELSE
            WRITE (20,2300)
          END IF
        END IF
        IF (nsetmax.gt.0) THEN
          WRITE (20,2400)
          WRITE (20,2450)
        END IF
        IF (ifos2d.eq.1.and.single.ne.1) WRITE (20,2500)
! assign weight units for m/kN or ft/lb unit systems
        IF (lengthunits.eq.'m'.and.uwtunits.eq.'kN/m^3') THEN
          wtunits = 'kg'
        ELSE 
          IF (lengthunits.eq.'ft'.and.uwtunits.eq.'lb/ft^3') THEN
            wtunits = 'lb'
          ELSE 
            wtunits = '  '
          END IF
        END IF
                         
        WRITE (20,2600)
        IF (single.ne.1) THEN
          WRITE (20,2700)
        ELSE
          WRITE(20,2750)
        END IF   

        IF (ofos.eq.nullhi) THEN
          WRITE (20,2810)        
        ELSE
          IF (method.eq.'b') THEN
            WRITE (20,2800)"Bishop's ", ofos  
            WRITE (20,2800) 'Ordinary ',felminfos
          ELSE
            WRITE (20,2800) 'Ordinary ',ofos
          END IF    
          IF (lengthunits.eq.'  ') THEN
            WRITE (20,2905) ovol
            WRITE (20,3005) oarea
            WRITE (20,3015) osliparea
          ELSE
            WRITE (20,2900) lengthunits,ovol
            WRITE (20,3000) lengthunits,oarea
            WRITE (20,3010) lengthunits,osliparea
          END IF 
          IF (wtunits.eq.'  ') THEN
            WRITE (20,3025) owt
          ELSE        
            WRITE (20,3020) wtunits,owt
          END IF
          WRITE (20,3100) ocol
          WRITE (20,3200) 
          WRITE (20,3300) ofsx+xll,ofsy+yll,ofsz,orad
          IF (irotcen.eq.1) THEN
            WRITE (20,3201)
            WRITE (20,3301) xcenrot+xll,ycenrot+yll,zcenrot
          END IF
!          WRITE (20,3350)
          WRITE (20,3400) oangle
          IF (ifos2d.eq.1) THEN
            WRITE (20,3500)
            IF (method.eq.'b') THEN
              WRITE (20,3600) "Bishop's ",ofos3d2d
            ELSE
              WRITE (20,3600)  'Ordinary ',ofos3d2d
            END IF    
            IF (lengthunits.eq.'  ') THEN 
              WRITE (20,3705) oarea3d2d
              WRITE (20,3805) owidth3d2d
              WRITE (20,3905) oarclength3d2d
            ELSE                               
              WRITE (20,3700) lengthunits,oarea3d2d
              WRITE (20,3800) lengthunits,owidth3d2d
              WRITE (20,3900) lengthunits,oarclength3d2d
            END IF  
            WRITE (20,3100) ocol3d2d
            WRITE (20,2600)
            IF (single.ne.1) THEN
              WRITE (20,4000)              
              IF (method.eq.'b') THEN
                WRITE (20,3600) "Bishop's ",ofos2d
              ELSE
                WRITE (20,3600)  'Ordinary ',ofos2d
             END IF                                          
              IF (lengthunits.eq.'  ') THEN
                WRITE (20,3705) oarea2d
                WRITE (20,3805) owidth2d
                WRITE (20,3905) oarclength            
              ELSE
                WRITE (20,3700) lengthunits,oarea2d
                WRITE (20,3800) lengthunits,owidth2d
                WRITE (20,3900) lengthunits,oarclength
              END IF
              WRITE (20,3100) ocol2d
              WRITE (20,3200)
              WRITE (20,3300) ofsx2d+xll,ofsy2d+yll,ofsz2d,orad2d
!              WRITE (20,3350)
              WRITE (20,3400) oangle2d
              WRITE (20,4100)
               IF (method.eq.'b') THEN
                  WRITE (20,2800) 'Bishops  ', ofos2d3d              
              ELSE
                  WRITE (20,2800) 'Ordinary ', ofos2d3d
              END IF    
              IF (lengthunits.eq.'  ') THEN
                WRITE (20,2905) ovol2d3d
                WRITE (20,3005) oarea2d3d
                WRITE (20,3015) osliparea2d3d            
              ELSE
                WRITE (20,2900) lengthunits,ovol2d3d
                WRITE (20,3000) lengthunits,oarea2d3d
                WRITE (20,3010) lengthunits,osliparea2d3d
              END IF
              IF (wtunits.eq.'  ') THEN
                WRITE (20,3025) owt2d3d
              ELSE
                WRITE (20,3020) wtunits,owt2d3d
              END IF
              WRITE (20,3100) ocol2d3d
              WRITE (20,2600)
            END IF
          END IF
          IF (remove.eq.'L'.and.ovfos.lt.111.0_pr) THEN 
            WRITE (20,4200)
            IF (method.eq.'b') THEN
              WRITE (20,2800) 'Bishops  ', ovfos
            ELSE                     
              WRITE (20,2800) 'Ordinary ',ovfos
            END IF
            IF (lengthunits.eq.'  ') THEN
              WRITE (20,2905) ovvol
              WRITE (20,3005) ovarea
              WRITE (20,3015) ovsliparea          
            ELSE
              WRITE (20,2900) lengthunits,ovvol
              WRITE (20,3000) lengthunits,ovarea
              WRITE (20,3010) lengthunits,ovsliparea
            END IF
            IF (wtunits.eq.'  ') THEN
              WRITE (20,3025) ovwt        
            ELSE
              WRITE (20,3020) wtunits,ovwt
            END IF
            WRITE (20,3200) 
            WRITE (20,3300) ovfsx+xll,ovfsy+yll,ovfsz,ovrad
!            WRITE (20,3350)
            WRITE (20,3400) ovangle
            WRITE (20,2600)
          END IF
        END IF  

        CALL DATE_AND_TIME(date,time)
        WRITE (20,4400) date(5:6),date(7:8),date(1:4),time(1:2),time(3:4),time(5:6)
        
!     Write output to other output files        

!     Write out header for array of slope 
        CALL WriteHeader(36,5900,rnull)  
        IF (ifailsurf.eq.1) THEN
          CALL WriteHeader(41,5900,rnull)
          CALL WriteHeader(42,5900,rnull)
        END IF  
        
        IF (ifos2d.ne.1) THEN
          IF (method.eq.'o') THEN
             CALL writelengthunits(21,7020)
           ELSE
            CALL writelengthunits(21,7030)  
           END IF 
        ELSE
          IF (method.eq.'o') THEN
            CALL writelengthunits(21,7120)
          ELSE
            CALL writelengthunits(21,7130)
          END IF  
          IF (single.ne.1) CALL writelengthunits(35,7220)
        END IF
        
        cardinality=4        

        WRITE (21,*) minfsx+xll,maxfsx+xll,cardinality
        WRITE (21,*) minfsy+yll,maxfsy+yll,cardinality        
        WRITE (21,*) minfsz,maxfsz,cardinality
        WRITE (21,*) 1,nx,cardinality
        WRITE (21,*) 1,ny,cardinality        
        WRITE (21,*) minfsrad,maxfsrad,cardinality
        WRITE (21,*) minfsangle,maxfsangle,cardinality
        WRITE (21,*) minmcol,maxmcol,cardinality
        WRITE (21,*) minfsvol,maxfsvol,cardinality
        WRITE (21,*) minfsarea,maxfsarea,cardinality
        WRITE (21,*) minfsmin,maxfsmin,cardinality
        IF (method.eq.'b')  WRITE (21,*) minfelfsmin,maxfelfsmin,cardinality        
        IF (single.eq.1.or.ifos2d.eq.1)  WRITE (21,*) minfsmin3d2d,maxfsmin3d2d,cardinality           
 
        IF (ifos2d.eq.1.and.single.ne.1) THEN
          WRITE (35,*) minfsx2d+xll,maxfsx2d+xll,cardinality
          WRITE (35,*) minfsy2d+yll,maxfsy2d+yll,cardinality     
          WRITE (35,*) minfsz2d,maxfsz2d,cardinality            
          WRITE (35,*) 1,nx,cardinality
          WRITE (35,*) 1,ny,cardinality        
          WRITE (35,*) minfsrad2d,maxfsrad2d,cardinality
          WRITE (35,*) minfsangle2d,maxfsangle2d,cardinality
          WRITE (35,*) minmcol2d,maxmcol2d,cardinality
          WRITE (35,*) minfsarea2d,maxfsarea2d,cardinality
          WRITE (35,*) minfsarclength,maxfsarclength,cardinality
          WRITE (35,*) minfswidth2d,maxfswidth2d,cardinality
          WRITE (35,*) minfsmin2d,maxfsmin2d,cardinality
          WRITE (35,*) minfsmin2d3d,maxfsmin2d3d,cardinality 
        END IF        
                        

        IF (single.ne.1) THEN 
!     Write out i,j,k,fos or x,y,z,fos for 3-D array.        
          IF (isubsurf.eq.1.or.isubsurf.eq.2) THEN
            xmax = nx*delxy + xll
            ymax = ny*delxy + yll
            zzmax = zbot + (nz-1)*delz
            SELECT CASE (isubsurf)
! write out header for ijk format            
             CASE (1)
              WRITE (23,8000) fosmeth
! write out header for xyz format              
             CASE (2)
              WRITE (23,8010) fosmeth                                                
            END SELECT 
            IF (lengthunits.eq.'  ') THEN
              WRITE (23,8020) 
            ELSE  
              WRITE (23,8025) lengthunits  
            END IF                       
            WRITE (23,8100) nx,ny,nz
            WRITE (23,8200) xll,xmax
            WRITE (23,8300) yll,ymax
            WRITE (23,8400) zbot-0.5*delz,zzmax+0.5*delz                  
          END IF
          IF (isubsurf.eq.3) THEN
            IF (LEN_TRIM(lengthunits).eq.1) THEN
              WRITE (23,8510) lengthunits,lengthunits,lengthunits,fosmeth
            ELSE
              IF (LEN_TRIM(lengthunits).eq.2) THEN
                WRITE (23,8520) lengthunits,lengthunits,lengthunits,fosmeth
              ELSE 
                WRITE (23,8500) fosmeth
              END IF  
            END IF                             
          END IF          

!     Write out header for arrays of fosvol or fosarea 
          CALL WriteHeader(27,5900,rnull)                  

!     Write out headers for arrays of critcheck, numcols, and boundcheck
          IF (isqout.eq.1) THEN
            CALL WriteHeader(22,5950,REAL(inull,pr))
            CALL WriteHeader(29,5950,REAL(inull,pr))
            IF (single.ne.1.and.srchfile.ne.1) CALL WriteHeader(44,5950,REAL(inull,pr))
          END IF
        END IF
!     If creating file newDEM.out
        IF (demflag.eq.1.and.ofos.lt.foscut)&    
          CALL WriteHeader(24,5900,rnull)        
!     Write out header for array of minimum 3-D FOS        
        CALL WriteHeader(26,5900,nullhi)
        
!!!   Fellenius comparison change        
        IF (method.eq.'b') CALL WriteHeader(43,5900,nullhi)
!!!
        IF (irelfos.eq.1) CALL WriteHeader(37,5900,nullhi) 
!     Write out header for array of minimum 2-D FOS        
        IF (ifos2d.eq.1.and.single.ne.1) THEN
          CALL WriteHeader(19,5900,nullhi)
          IF (irelfos.eq.1) CALL WriteHeader(38,5900,nullhi)
        END IF         
!     If surfaces were filtered create file of count of filtered nodes
        IF (mfile.eq.1) CALL WriteHeader(17,5950,REAL(inull,pr)) 
   
        IF (single.ne.1.and.srchfile.ne.1) THEN    
          xsmin = REAL(ismin-1,pr)*delxy + halfdelxy
          xsmax = REAL( ismin+ (int((ismax-ismin)/nsrchres)*nsrchres ) - 1,pr ) * delxy+ halfdelxy
          ysmin = REAL(jsmin-1,pr)*delxy+halfdelxy
          ysmax = REAL( jsmin+ (int((jsmax-jsmin)/nsrchres)*nsrchres ) - 1,pr ) * delxy+ halfdelxy
          zsmaxtrue =  REAL( zsmin+ (int((zsmax-zsmin)/zsrchres)*zsrchres ) ,pr )
        END IF   

!     Convert slope from radians to degrees before writing values to file 
!      WHERE (surfslope.ne.rnull) surfslope= surfslope* 57.2957795_pr  ! this crashes with very large arrays
        DO j = ny,1,-1
          DO i = nx,1,-1
             if (surfslope(i,j).ne.rnull) surfslope(i,j) = surfslope(i,j) * 57.2957795_pr
         END DO  
        END DO      

!     Write array data                        
        DO j = ny,1,-1
          CALL writeslope(36,j)
          IF (ifailsurf.eq.1) THEN
            CALL writefailslope(41,j)
            CALL writefaildepth(42,j)
          END IF
          IF (ofos.ne.nullhi) THEN
            DO i = 1,nx
              jj = (ny-j)+1
              IF (fsmin(i,jj).ne.nullhi) THEN
                CALL writefs(i,jj)
                IF ((isubsurf.eq.1.or.isubsurf.eq.2.or.isubsurf.eq.3).and.single.ne.1)&
                 CALL writesubsurf(i,j)
              END IF
            END DO          
            Call writefsmin(26,j) 
!!!   Fellenius comparison change
            IF (method.eq.'b') Call writefelfsmin(43,j)
!!!            
            IF (irelfos.eq.1) CALL writefsrelmin(37,j) 
            IF (single.ne.1) CALL writefsvolarea(27,j)
            IF (ifos2d.eq.1.and.single.ne.1) THEN
               CALL writefsmin2d(19,j)
               IF (irelfos.eq.1) CALL writefsrelmin2d(38,j)
            END IF   
            IF (demflag.eq.1.and.ofos.lt.foscut) CALL writezdem(24,j)
            IF (isqout.eq.1.and.single.ne.1) THEN
              CALL writeicrit(22,j)
              CALL writemcol(29,j)
               IF (single.ne.1.and.srchfile.ne.1) CALL writeboundcheck(44,j,xsmin,xsmax,ysmin,ysmax,zsmaxtrue)
            END IF
            IF (mfile.eq.1) CALL writefiltcount(17,j)
          END IF
        END DO

        CLOSE (15)
        CLOSE (17)
        CLOSE (19)
        CLOSE (21)
        CLOSE (22)
        CLOSE (23)
        CLOSE (24)
        CLOSE (26)
        CLOSE (27)
        CLOSE (29)
        CLOSE (33)
        CLOSE (34)
        CLOSE (35)
        CLOSE (36)
        CLOSE (37)
        CLOSE (38)
        
 !!!   Fellenius comparison change
        IF (method.eq.'b') CLOSE (43)
!!! 

1000    FORMAT (/,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')        
1100    FORMAT ('III. OUTPUT FILES GENERATED: ')
1150    FORMAT (/,'     LOCATION FOR OUTPUT FILES:',/,'          ',A,/)
1200    FORMAT (/,'     Optional files generated:')
1300    FORMAT ('      ',A)
1400    FORMAT ('IV. RESULTS:')
1500    FORMAT ('Number of trial surfaces tried:                                     ',i10)
1600    FORMAT ('F < foscut found and newDEM_out file created?                              yes')
1700    FORMAT ('F < foscut found and newDEM_out file created?                               no')
1900    FORMAT ('Number of surfaces eliminated due to filters and/or nonconvergence: ',i10)
2000    FORMAT ('     Check file filter_out for detailed information.')
2100    FORMAT ('Number of surfaces with active column totals less than limcol:        ',i8)
2200    FORMAT ('     Check file ncolerr_out for detailed information.')
2300    FORMAT ('NOTE: Single surface has fewer columns than column limit.')
2400    FORMAT ('More subsets than nnset in some locations. Consider increasing nnset')
2450    FORMAT (' in Module SetData.')
2500    FORMAT ('2D factors of safety also computed')
2600    FORMAT ('------------')
2700    FORMAT ('3D POTENTIAL FAILURE - GLOBAL MINIMUM')
2750    FORMAT ('3D POTENTIAL FAILURE')
2800    FORMAT (A9,'3D factor of safety:                                       ',f10.4)
2810    FORMAT ('3D factor of safety:                                        NO VALID SURFACES')
2900    FORMAT ('Volume (',A2,'^3):                                                    ',es12.5)
2905    FORMAT ('Volume:                                                           ',es12.5)
3000    FORMAT ('Horizontal surface area (',A2,'^2):                                   ',es12.5)
3005    FORMAT ('Horizontal surface area:                                          ',es12.5)
3010    FORMAT ('Slip surface area (',A2,'^2):                                         ',es12.5)
3015    FORMAT ('Slip surface area:                                                ',es12.5)
3020    FORMAT ('Weight (',A2,'):                                                      ',es12.5)
3025    FORMAT ('Weight:                                                           ',es12.5)
3100    FORMAT ('Number of active columns:                                               ',i6)
3200    FORMAT ('        x-center         y-center          z-center     radius')
3201    FORMAT ('    x-rot. center    y-rot. center       z-rot. center')
3300    FORMAT (3(f16.4,'  '),es12.5)
3301    FORMAT (3(f16.4,'  '))
3350    FORMAT ('local coordinates are given as distance from lower left corner')
3400    FORMAT ('Slip direction, relative to search lattice:                         ',f10.4)
3500    FORMAT (/,'Associated 2D potential failure:')
3600    FORMAT (A9,'2D factor of safety:                                       ',f10.4)
3700    FORMAT ('Cross-sectional area ( ',A2,'^2):                                     ',es12.5)
3705    FORMAT ('Cross-sectional area:                                             ',es12.5)
3800    FORMAT ('Horizontal length ( ',A2,'):                                          ',es12.5)
3805    FORMAT ('Horizontal length:                                                ',es12.5)
3900    FORMAT ('Arc length ( ',A2,'):                                                 ',es12.5)
3905    FORMAT ('Arc length:                                                       ',es12.5)
4000    FORMAT ('2D POTENTIAL FAILURE - GLOBAL MINIMUM')
4100    FORMAT (/,'Associated 3D potential failure:')
4200    FORMAT ('LARGEST VOLUME FAILURE WITH F < FOSCUT')
4300    FORMAT ('Retrogression # ',i5,/)      
4400    FORMAT ('End date and time: ',A2,'/',A2,'/',A4,'  ',A2,':',A2,':',A2)
4600    FORMAT ('        Range: [',f10.4,',',f10.4,']')
4700    FORMAT ('        Range: [',es12.4,',',es12.4,']')
4800    FORMAT ('        Range: [',f11.4,',',f11.4,']')
4900    FORMAT ('        Range: [',i5,',',i5,']')   
5000    FORMAT (10000f11.4) 
5100    FORMAT (10000f13.4)
5200    FORMAT (10000i6)
5250    FORMAT (10000(1x,i11))
5300    FORMAT (10000es12.4) 
5500    FORMAT ('ncols         ',i5)
5600    FORMAT ('nrows         ',i5)
5700    FORMAT (A)
5750    FORMAT ('xllcorner     ',f15.6)
5760    FORMAT ('yllcorner     ',f15.6)
5800    FORMAT ('cellsize      ',f9.4)
5900    FORMAT ('NODATA_value  ',f10.4)
5950    FORMAT ('NODATA_value  ',i6)  
6000    FORMAT (5(f15.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',2(es10.3,' '),f10.4)
6100    FORMAT (5(f15.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',2(es10.3,' '),2f10.4)
6200    FORMAT (5(f15.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',3(es10.3,' '),2f10.4)
6300    FORMAT (3i5,f10.4)
6310    FORMAT (3f20.3,f10.4)
6320    FORMAT (3f30.4,f11.4)
8000    FORMAT ('# Comments: Scoops3D Factors of safety in the subsurface',/,&
        &'# Coordinate_system: ijk',/,'# Field: 1 i',/,'# Field: 2 j',/,&
        &'# Field: 3 k',/,'# Field: 4 F_',A4)
8010    FORMAT ('# Comments: Scoops3D Factors of safety in the subsurface',/,&
        &'# Coordinate_system: xyz',/,'# Field: 1 x',/,'# Field: 2 y',/,&
        &'# Field: 3 z',/,'# Field: 4 F_',A4)
8020    FORMAT ('# Length_units: UNDEFINED') 
8025    FORMAT ('# Length_units: ',A2)          
8100    FORMAT ('# Grid_size: ',i6,' x ',i6,' x ',i6)
8200    FORMAT ('# Grid_x_range: ',f13.3,' to ',f12.3)
8300    FORMAT ('# Grid_y_range: ',f12.3,' to ',f12.3)
8400    FORMAT ('# Grid_z_range: ',f12.3,' to ',f12.3)
8500    FORMAT ('x y z F_',A4)       
8510    FORMAT ('x_',A1,' y_',A1,' z_',A1,' F_',A4)               
8520    FORMAT ('x_',A2,' y_',A2,' z_',A2,' F_',A4)        

        END SUBROUTINE writeout
        

        
