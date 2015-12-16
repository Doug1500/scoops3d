        SUBROUTINE Arraystats
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine determines minimum and maximum output array values
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       Called by Writeout
!
!       
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
        USE CommonData
        USE ResultsStats
        USE GridData, ONLY: nx,ny,zdem
        USE FOSData, ONLY: irelfos,ifos2d,fsmin,fsvol,fsarea,mcol,&
                filtcount,fsmin2d,fsmin3d2d,fsmin2d3d,fsx2d,fsy2d,fsz2d,&
                fsrad2d,fsangle2d,fsarclength,fswidth2d,fsx,fsy,fsz,fsrad,&
                fsangle,ofos,ofos2d,mcol2d,method,fsarea2d,felfsmin
        USE SetData, ONLY: surfslope
        USE FailSurfData, ONLY: ifailsurf,failsurfslope,failsurfdepth
        
        IMPLICIT NONE
   	
        numdatapts = COUNT(fsmin.ne.nullhi)
        minslope = MINVAL(surfslope,MASK=surfslope.ne.rnull)
        maxslope = MAXVAL(surfslope,MASK=surfslope.ne.rnull)
        IF (ifailsurf.eq.1) THEN
          minfailslope = MINVAL(failsurfslope,MASK=failsurfslope.ne.rnull)
          maxfailslope = MAXVAL(failsurfslope,MASK=failsurfslope.ne.rnull)
          minfaildepth = MINVAL(failsurfdepth,MASK=failsurfdepth.ne.rnull)
          maxfaildepth = MAXVAL(failsurfdepth,MASK=failsurfdepth.ne.rnull)
        END IF
        minfsmin = MINVAL(fsmin,MASK=fsmin.ne.nullhi)
        maxfsmin = MAXVAL(fsmin,MASK=fsmin.ne.nullhi)
        IF (method.eq.'b') THEN
          minfelfsmin = MINVAL(felfsmin,MASK=felfsmin.ne.nullhi)
          maxfelfsmin = MAXVAL(felfsmin,MASK=felfsmin.ne.nullhi)
        END IF 
        IF (irelfos.eq.1) THEN
          minfsrelmin = MINVAL(fsmin/ofos,MASK=fsmin.ne.nullhi)
          maxfsrelmin = MAXVAL(fsmin/ofos,MASK=fsmin.ne.nullhi)
        END IF
        minfsvol = MINVAL(fsvol,MASK=fsvol.ne.rnull)
        maxfsvol = MAXVAL(fsvol,MASK=fsvol.ne.rnull)
        minfsarea = MINVAL(fsarea,MASK=fsarea.ne.rnull)
        maxfsarea = MAXVAL(fsarea,MASK=fsarea.ne.rnull)
        minzdem = MINVAL(zdem,MASK=zdem.ne.rnull)
        maxzdem = MAXVAL(zdem,MASK=zdem.ne.rnull)    
        minmcol = MINVAL(mcol,MASK=mcol.ne.inull)
        maxmcol = MAXVAL(mcol,MASK=mcol.ne.inull)
        minfiltcount = MINVAL(filtcount,MASK=fsmin.ne.inull)
        maxfiltcount = MAXVAL(filtcount,MASK=fsmin.ne.inull)
        IF (ifos2d.eq.1) THEN
          minfsmin2d = MINVAL(fsmin2d,MASK=fsmin2d.ne.nullhi)
          maxfsmin2d = MAXVAL(fsmin2d,MASK=fsmin2d.ne.nullhi)
          IF (irelfos.eq.1) THEN
            minfsrelmin2d = MINVAL(fsmin2d/ofos2d,MASK=fsmin2d.ne.nullhi)
            maxfsrelmin2d = MAXVAL(fsmin2d/ofos2d,MASK=fsmin2d.ne.nullhi)
          END IF
          minfsmin3d2d = MINVAL(fsmin3d2d,MASK=fsmin3d2d.ne.nullhi)
          maxfsmin3d2d = MAXVAL(fsmin3d2d,MASK=fsmin3d2d.ne.nullhi)
          minfsmin2d3d = MINVAL(fsmin2d3d,MASK=fsmin2d3d.ne.nullhi)
          maxfsmin2d3d = MAXVAL(fsmin2d3d,MASK=fsmin2d3d.ne.nullhi)
          minfsz2d = MINVAL(fsz2d,MASK=fsz2d.ne.rnull)
          maxfsz2d = MAXVAL(fsz2d,MASK=fsz2d.ne.rnull)
          minfsx2d = MINVAL(fsx2d,MASK=fsx2d.ne.rnull)
          maxfsx2d = MAXVAL(fsx2d,MASK=fsx2d.ne.rnull)
          minfsy2d = MINVAL(fsy2d,MASK=fsy2d.ne.rnull)
          maxfsy2d = MAXVAL(fsy2d,MASK=fsy2d.ne.rnull)
          minfsrad2d = MINVAL(fsrad2d,MASK=fsrad2d.ne.rnull)
          maxfsrad2d = MAXVAL(fsrad2d,MASK=fsrad2d.ne.rnull)
          minfsangle2d = MINVAL(fsangle2d,MASK=fsangle2d.ne.rnull)
          maxfsangle2d = MAXVAL(fsangle2d,MASK=fsangle2d.ne.rnull)
          minfsarclength = MINVAL(fsarclength,MASK=fsarclength.ne.rnull)
          maxfsarclength = MAXVAL(fsarclength,MASK=fsarclength.ne.rnull)
          minfswidth2d = MINVAL(fswidth2d,MASK=fswidth2d.ne.rnull)
          maxfswidth2d = MAXVAL(fswidth2d,MASK=fswidth2d.ne.rnull)
          minfsarea2d = MINVAL(fsarea2d,MASK=fsarea2d.ne.rnull)
          maxfsarea2d = MAXVAL(fsarea2d,MASK=fsarea2d.ne.rnull)
          minmcol2d = MINVAL(mcol2d,MASK=mcol2d.ne.rnull)
          maxmcol2d = MAXVAL(mcol2d,MASK=mcol2d.ne.rnull)
        END IF
        minfsz = MINVAL(fsz,MASK=fsz.ne.rnull)
        maxfsz = MAXVAL(fsz,MASK=fsz.ne.rnull)
        minfsx = MINVAL(fsx,MASK=fsx.ne.rnull)
        maxfsx = MAXVAL(fsx,MASK=fsx.ne.rnull)
        minfsy = MINVAL(fsy,MASK=fsy.ne.rnull)
        maxfsy = MAXVAL(fsy,MASK=fsy.ne.rnull)
        minfsrad = MINVAL(fsrad,MASK=fsrad.ne.rnull)
        maxfsrad = MAXVAL(fsrad,MASK=fsrad.ne.rnull)
        minfsangle = MINVAL(fsangle,MASK=fsangle.ne.rnull)
        maxfsangle = MAXVAL(fsangle,MASK=fsangle.ne.rnull)  
        
        END SUBROUTINE Arraystats
