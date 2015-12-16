
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx        
        SUBROUTINE OpenFileWrite(opennum,writenum,namefile,fileext,charetro)
 
        USE CommonData
               
        IMPLICIT NONE

        INTEGER, INTENT(in) :: opennum,writenum
        CHARACTER*1, INTENT(in) :: charetro 
        CHARACTER*3, INTENT(in) :: fileext
        CHARACTER*(*), INTENT(in) :: namefile  

        PRINT *,outputdir(1:LEN_TRIM(outputdir))//&
                    filin(bfile:nfile)//'_'//TRIM(namefile)//'_out.'//TRIM(fileext)
        IF (charetro.eq."0") THEN
          OPEN (opennum,STATUS = 'replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                    filin(bfile:nfile)//'_'//TRIM(namefile)//'_out.'//TRIM(fileext))
          WRITE (writenum,1300) filin(bfile:nfile)//'_'//TRIM(namefile)//'_out.'//TRIM(fileext)
        ELSE
          OPEN (opennum,STATUS = 'replace',file = outputdir(1:LEN_TRIM(outputdir))//&
                    filin(bfile:nfile)//'_'//TRIM(namefile)//charetro//'_out.'//TRIM(fileext))
          WRITE (writenum,1300) filin(bfile:nfile)//'_'//TRIM(namefile)//charetro//'_out.'//TRIM(fileext)
        END IF
          
1300    FORMAT ('      ',A)

        END SUBROUTINE OpenFileWrite

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx        
        SUBROUTINE writeHeader(filenum,nullfmt,nullnum)
 
        USE CommonData
        USE GridData, ONLY: nx,ny,delxy,xllcorner,yllcorner
        
        IMPLICIT NONE

        INTEGER, INTENT(in) :: filenum,nullfmt
        REAL(pr), INTENT (in) :: nullnum
        
        WRITE (filenum,5500) nx
        WRITE (filenum,5600) ny
        WRITE (filenum,5700) xllcorner
        WRITE (filenum,5700) yllcorner
        WRITE (filenum,5800) delxy
        IF (nullfmt.eq.5900) WRITE (filenum,5900) nullnum
        IF (nullfmt.eq.5950) WRITE (filenum,5950) NINT(nullnum)
        
5500    FORMAT ('ncols         ',i5)
5600    FORMAT ('nrows         ',i5)
5700    FORMAT (A)
5750    FORMAT ('xllcorner     ',f15.6)
5760    FORMAT ('yllcorner     ',f15.6)
5800    FORMAT ('cellsize      ',f9.4)
5900    FORMAT ('NODATA_value  ',f10.4)
5950    FORMAT ('NODATA_value  ',i6)       

        END SUBROUTINE writeHeader

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writeslope(filenum,j)
        
        USE CommonData
        USE GridData, ONLY: nx
        USE SetData, ONLY: surfslope        
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
       
        WRITE (filenum,5000) (surfslope(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writeslope
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writefailslope(filenum,j)
        
        USE CommonData
        USE GridData, ONLY: nx
        USE FailSurfData, ONLY: failsurfslope  
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
 
        WRITE (filenum,5000) (failsurfslope(i,j)*57.2957795_pr,i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefailslope
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writefaildepth(filenum,j)
        
        USE CommonData
        USE GridData, ONLY: nx
        USE FailSurfData, ONLY: failsurfdepth  
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
 
        WRITE (filenum,5000) (failsurfdepth(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefaildepth
         

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writefsmin(filenum,j)
        
        USE GridData, ONLY: nx
        USE FOSData, ONLY: fsmin
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
        WRITE (filenum,5000) (fsmin(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefsmin
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writefsrelmin(filenum,j)

        USE CommonData        
        USE GridData, ONLY: nx
        USE FOSData, ONLY: ofos,fsmin
        USE SetData, ONLY: zb  ! used as dummy array to save memory
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j

          DO i = 1,nx
            IF (fsmin(i,j).ne.nullhi) THEN
              zb(i,j) = fsmin(i,j)/ofos
            ELSE
              zb(i,j) = nullhi  
            END IF
          END DO
        
          WRITE (filenum,5000) (zb(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefsrelmin         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writefsvolarea(filenum,j)
        
        USE GridData, ONLY: nx
        USE FOSData, ONLY: fsvol,fsarea
        USE SearchData, ONLY: vacriterion
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
        IF (vacriterion.eq.'v') THEN  
          WRITE (filenum,5300) (fsvol(i,j),i=1,nx)
        ELSE
          WRITE (filenum,5300) (fsarea(i,j),i=1,nx)
        END IF
         
5300    FORMAT (10000es12.4) 

        END SUBROUTINE writefsvolarea  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE writefsmin2d(filenum,j)
        
        USE GridData, ONLY: nx
        USE FOSData, ONLY: fsmin2d
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
        WRITE (filenum,5000) (fsmin2d(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefsmin2d        
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  

        SUBROUTINE writefsrelmin2d(filenum,j)
        
        USE CommonData
        USE GridData, ONLY: nx
        USE FOSData, ONLY: ofos2d,fsmin2d
        USE SetData, ONLY: zmid  ! used as dummy array to save memory
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
          DO i = 1,nx
            IF (fsmin2d(i,j).ne.nullhi) THEN
              zmid(i,j) = fsmin2d(i,j)/ofos2d
            ELSE
              zmid(i,j) = nullhi  
            END IF
          END DO
          WRITE (filenum,5000) (zmid(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefsrelmin2d         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    

       SUBROUTINE writezdem(filenum,j)
        
        USE GridData, ONLY: nx,zdem
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
        WRITE (filenum,5100) (zdem(i,j),i=1,nx)
         
5100    FORMAT (10000f13.4)

        END SUBROUTINE writezdem         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                

       SUBROUTINE writeicrit(filenum,j)
       
        USE CommonData 
        USE GridData, ONLY: nx
        USE FOSData, ONLY: icrit,fsmin
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
        DO i = 1,nx
          IF (fsmin(i,j).eq.nullhi) icrit(i,j) = inull
        END DO
        
        WRITE (filenum,5200) (icrit(i,j),i=1,nx)
         
5200    FORMAT (10000i6)

        END SUBROUTINE writeicrit         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                

       SUBROUTINE writemcol(filenum,j)
        
        USE GridData, ONLY: nx
        USE FOSData, ONLY: mcol
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j

        WRITE (filenum,5250) (mcol(i,j),i=1,nx)
         
5250    FORMAT (10000(1x,i11))

        END SUBROUTINE writemcol         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  
     SUBROUTINE writeboundcheck(filenum,j,xsmin,xsmax,ysmin,ysmax,zsmaxtrue)
 
        USE CommonData       
        USE GridData, ONLY: nx,delxy,halfdelxy,xll,yll
        USE FOSData, ONLY: fsx,fsy,fsz,fsmin
        USE SearchData, ONLY: zsmin        
        IMPLICIT NONE
        
        INTEGER :: i,bound(nx)
        INTEGER, INTENT(in) :: filenum,j
        REAL(pr), INTENT(in) :: xsmin,xsmax,ysmin,ysmax,zsmaxtrue

        bound = 0
        DO i = 1,nx 
          IF (fsmin(i,j).ne.nullhi) THEN
            IF (fsx(i,j).eq.xsmin) bound(i) = 100
            IF (fsx(i,j).eq.xsmax) bound(i) = 900
            IF (fsy(i,j).eq.ysmin) bound(i) = bound(i)+10
            IF (fsy(i,j).eq.ysmax) bound(i) = bound(i)+90
            IF (fsz(i,j).eq.zsmin) bound(i) = bound(i)+1
            IF (fsz(i,j).eq.zsmaxtrue) bound(i) = bound(i)+9
          ELSE  
            bound(i) = inull 
          END IF  
        END DO
 
        WRITE (filenum,5200) (bound(i),i=1,nx)
         
5200    FORMAT (10000i6)

        END SUBROUTINE writeboundcheck         

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                 

       SUBROUTINE writefiltcount(filenum,j)
 
        USE CommonData       
        USE GridData, ONLY: nx
        USE FOSData, ONLY: filtcount,fsmin
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
       
        DO i = 1,nx
          IF (fsmin(i,j).eq.nullhi) filtcount(i,j) = inull
        END DO 
        
        WRITE (filenum,5200) (filtcount(i,j),i=1,nx)
         
5200    FORMAT (10000i6)

        END SUBROUTINE writefiltcount         

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     



       SUBROUTINE writelengthunits(filenum,formatnum)
        
        USE GridData, ONLY: lengthunits
        
        IMPLICIT NONE
        
        INTEGER, INTENT(in) :: filenum,formatnum
        
        IF (formatnum.eq.7020) THEN
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (filenum,7021) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                  &lengthunits
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (filenum,7022) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                    &lengthunits
            ELSE
              WRITE (filenum,7020) 
            END IF  
          END IF
        END IF
 
        IF (formatnum.eq.7030) THEN
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (filenum,7031) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                  &lengthunits
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (filenum,7032) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                    &lengthunits
            ELSE
              WRITE (filenum,7030) 
            END IF  
          END IF
        END IF
        
        IF (formatnum.eq.7220) THEN
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (filenum,7221) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                  &lengthunits,lengthunits
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (filenum,7222) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                    &lengthunits,lengthunits
            ELSE
              WRITE (filenum,7220) 
            END IF  
          END IF
        END IF
        
        IF (formatnum.eq.7120) THEN
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (filenum,7121) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                  &lengthunits
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (filenum,7122) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                    &lengthunits
            ELSE
              WRITE (filenum,7120) 
            END IF  
          END IF
        END IF
        
        IF (formatnum.eq.7130) THEN
          IF (LEN_TRIM(lengthunits).eq.1) THEN
            WRITE (filenum,7131) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                  &lengthunits
          ELSE
            IF (LEN_TRIM(lengthunits).eq.2) THEN
              WRITE (filenum,7132) lengthunits,lengthunits,lengthunits,lengthunits,lengthunits,&
                    &lengthunits
            ELSE
              WRITE (filenum,7130) 
            END IF  
          END IF
        END IF          
         
7020    FORMAT ('xcen',/,'ycen',/,'zcen',/,&
                &'i',/,'j',/,'radius',/,'angle',/,'cols',/,&
                & 'vol',/,'area',/,'F_Ord')    
7021    FORMAT ('xcen_',A1,/,'ycen_',A1,/,'zcen_',A1,/,&
                &'i',/,'j',/,'radius_',A1,/,'angle',/,'cols',/,&
                & 'vol_',A1,'^3',/,'area_',A1,'^2',/,'F_Ord')                               
7022    FORMAT ('xcen_',A2,/,'ycen_',A2,/,'zcen_',A2,/,&
                &'i',/,'j',/,'radius_',A2,/,'angle',/,'cols',/,&
                & 'vol_',A2,'^3',/,'area_',A2,'^2',/,'F_Ord')  
                
7030    FORMAT ('xcen',/,'ycen',/,'zcen',/,&
                &'i',/,'j',/,'radius',/,'angle',/,'cols',/,&
                & 'vol',/,'area',/,'F_Bish',/,'F_Ord')    
7031    FORMAT ('xcen_',A1,/,'ycen_',A1,/,'zcen_',A1,/,&
                &'i',/,'j',/,'radius_',A1,/,'angle',/,'cols',/,&
                & 'vol_',A1,'^3',/,'area_',A1,'^2',/,'F_Bish',/,'F_Ord')                               
7032    FORMAT ('xcen_',A2,/,'ycen_',A2,/,'zcen_',A2,/,&
                &'i',/,'j',/,'radius_',A2,/,'angle',/,'cols',/,&
                & 'vol_',A2,'^3',/,'area_',A2,'^2',/,'F_Bish',/,'F_Ord')  
                
7120    FORMAT ('xcen',/,'ycen',/,'zcen',/,'i',/,&
                &'j',/,'radius',/,'angle',/,'cols',/,&
                &'vol',/,'area',/,'F_Ord_3D',/,'F_Ord_2D')                                   
7121    FORMAT ('xcen_',A1,/,'ycen_',A1,/,'zcen_',A1,/,&
                &'i',/,'j',/,'radius_',A1,/,'angle',/,'cols',/,&
                & 'vol_',A1,'^3',/,'area_',A1,'^2',/,'F_Ord_3D',/,'F_Ord_2D')     
7122    FORMAT ('xcen_',A2,/,'ycen_',A2,/,'zcen_',A2,/,&
                &'i',/,'j',/,'radius_',A2,/,'angle',/,'cols',/,&
                & 'vol_',A2,'^3',/,'area_',A2,'^2',/,'F_Ord_3D',/,'F_Ord_2D')    

7130    FORMAT ('xcen',/,'ycen',/,'zcen',/,'i',/,&
                &'j',/,'radius',/,'angle',/,'cols',&
                & /,'vol',/,'area',/,'F_Bish_3D',/,'F_Ord_3D',/,'F_Bish_2D')                                   
7131    FORMAT ('xcen_',A1,/,'ycen_',A1,/,'zcen_',A1,/,&
                &'i',/,'j',/,'radius_',A1,/,'angle',/,'cols',&
                & /,'vol_',A1,'^3',/,'area_',A1,'^2',/,'F_Bish_3D',/,'F_Ord_3D',/,'F_Bish_2D')     
7132    FORMAT ('xcen_',A2,/,'ycen_',A2,/,'zcen_',A2,/,&
                &'i',/,'j',/,'radius_',A2,/,'angle',/,'cols',&
                & /,'vol_',A2,'^3',/,'area_',A2,'^2',/,'F_Bish_3D',/,'F_Ord_3D',/,'F_Bish_2D')        

7220    FORMAT ('xcen',/,'ycen',/,'zcen',/,'i',/,&
                &'j',/,'radius',/,'angle',/,'cols',/,'area',&
                &/,'arclength',/,'width',/,'F_2D',/,'F_3D')                             
7221    FORMAT ('xcen_',A1,/,'ycen_',A1,/,'zcen_',A1,/,&
                &'i',/,'j',/,'radius_',A1,/,'angle',/,'cols',/,'area_',&
                &A1,'^2',/,'arclength_',A1,/,'width_',&
                &A1,/,'F_2D',/,'F_3D')   
7222    FORMAT ('xcen_',A2,/,'ycen_',A2,/,'zcen_',A2,/,&
                &'i',/,'j',/,'radius_',A2,/,'angle',/,'cols',/,'area_',&
                &A2,'^2',/,'arclength_',A2,/,'width_',&
                &A2,/,'F_2D',/,'F_3D')                   

        END SUBROUTINE writelengthunits         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   

        SUBROUTINE writesubsurf(i,j)
       
        USE CommonData
        USE GridData, ONLY: xll,yll,delxy,delz,nz
        USE FOSData, ONLY: isubsurf,zfos,zbot
        
        IMPLICIT NONE
        
        INTEGER :: k
        INTEGER, INTENT(in) :: i,j
        REAL :: x,y,z
        
        DO k = 1,nz 
          IF (zfos(i,j,k).ne.rnull) THEN
            SELECT CASE (isubsurf)
              CASE (1)
                WRITE (23,6300) i,j,k,zfos(i,j,k)                   
              CASE (2,3)
                x = REAL(i-0.5,pr)*delxy + xll
                y = REAL(j-0.5,pr)*delxy + yll
                z = zbot + (REAL(k-1,pr)*delz)
                WRITE (23,6310) x,y,z,zfos(i,j,k)
             END SELECT    
          END IF       
        END DO
         
6300    FORMAT (3i5,f10.4)
6310    FORMAT (3f20.3,f10.4)

        END SUBROUTINE writesubsurf         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     
        SUBROUTINE writefs(i,jj)
        
        USE CommonData
        USE GridData, ONLY: xll,yll
        USE FOSData, ONLY: method,fsx,fsy,fsz,fsx,fsy,fsrad,fsangle,mcol,&
                        fsvol,fsarea,fsmin,fsmin3d2d,fsx2d,fsy2d,& 
                        fsz2d,fsx2d,fsy2d,fsrad2d,fsangle2d,mcol2d,&
                        fsarea2d,fsarclength,fswidth2d,&
                        fsmin2d,felfsmin,fsmin2d3d,ifos2d,single
                        
        IMPLICIT NONE
        
        INTEGER, INTENT(in) :: i,jj
        
        IF (ifos2d.ne.1) THEN  
          IF (method.eq.'o') THEN
              WRITE (21,6000) fsx(i,jj)+xll,fsy(i,jj)+yll,fsz(i,jj),i,jj,&                 
                  fsrad(i,jj),fsangle(i,jj),mcol(i,jj),&
                  fsvol(i,jj),fsarea(i,jj),fsmin(i,jj)     
            ELSE
              WRITE (21,6100) fsx(i,jj)+xll,fsy(i,jj)+yll,fsz(i,jj),i,jj,&                 
                  fsrad(i,jj),fsangle(i,jj),mcol(i,jj),&
                  fsvol(i,jj),fsarea(i,jj),fsmin(i,jj),felfsmin(i,jj)     
            END IF 
        ELSE
           IF (method.eq.'o') THEN         
              WRITE (21,6100) fsx(i,jj)+xll,fsy(i,jj)+yll,fsz(i,jj),i,jj,&                  
                  fsrad(i,jj),fsangle(i,jj),mcol(i,jj),&
                  fsvol(i,jj),fsarea(i,jj),fsmin(i,jj),fsmin3d2d(i,jj)      
            ELSE  
              WRITE (21,6150) fsx(i,jj)+xll,fsy(i,jj)+yll,fsz(i,jj),i,jj,&                  
                  fsrad(i,jj),fsangle(i,jj),mcol(i,jj),&
                  fsvol(i,jj),fsarea(i,jj),fsmin(i,jj),felfsmin(i,jj),&
                  fsmin3d2d(i,jj)      
            END IF
          IF (fsmin2d(i,jj).ne.nullhi.and.single.ne.1) THEN
              WRITE (35,6200) fsx2d(i,jj)+xll,fsy2d(i,jj)+yll,fsz2d(i,jj),i,jj,&                  
                  fsrad2d(i,jj),fsangle2d(i,jj),mcol2d(i,jj),&
                  fsarea2d(i,jj),fsarclength(i,jj),fswidth2d(i,jj),&
                  fsmin2d(i,jj),fsmin2d3d(i,jj)
          END IF
        END IF
         
6000    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',2(es10.3,' '),f10.4)
6100    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',2(es10.3,' '),2f10.4)
6150    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',2(es10.3,' '),3f10.4)
6200    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',f8.2,' ',i11,' ',3(es10.3,' '),2f10.4)
6005    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',3(f16.4,' '),' ',f8.2,' ',i11,' ',2(es10.3,' '),f10.4)
6105    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',3(f16.4,' '),' ',f8.2,' ',i11,' ',2(es10.3,' '),2f10.4)
6155    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',3(f16.4,' '),' ',f8.2,' ',i11,' ',2(es10.3,' '),3f10.4)
6205    FORMAT (3(f16.4,' '),i4,' ',i4,' ',es10.3,' ',3(f16.4,' '),' ',f8.2,' ',i11,' ',3(es10.3,' '),2f10.4)

        END SUBROUTINE writefs         
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  
!!!   Fellenius comparison addition

        SUBROUTINE writefelfsmin(filenum,j)
        
        USE GridData, ONLY: nx
        USE FOSData, ONLY: felfsmin
        
        IMPLICIT NONE
        
        INTEGER :: i
        INTEGER, INTENT(in) :: filenum,j
        
        WRITE (filenum,5000) (felfsmin(i,j),i=1,nx)
         
5000    FORMAT (10000f11.4)

        END SUBROUTINE writefelfsmin
         
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                       
