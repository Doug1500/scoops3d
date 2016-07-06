        SUBROUTINE Subsets (newrad,nset,kset,zavg,ierr)
        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       This subroutine determines all contiguous sets of columns in the 
!       DEM that are intersected by the trial slip surface.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!      Called by Scoops3D
!
!      VARIABLES
!
!      allin -- flag to indicate if include columns are contained by slip surface
!      colcount -- number of columns intersected by sphere
!      delxy -- DEM grid resolution (delta x, delta y)
!      goodset -- flag for finding valid set for single surface
!      ierr -- error flag indicating initial radius problem (2) 
!              or too many subsets (1)
!      ifailsurf -- flag for whether failure surface file is used
!      in(nx,ny) -- number of contiguous nodes of each column bounded by sphere,
!         equals -1 if truncated surface, -2 if column contains boundary node.
!      includegrid -- Defines columns of the DEM that must be included in all
!           potential failure surfaces (1 if included, 0 if not)
!      insphere(nx,ny) -- indicates whether DEM column node is bounded by sphere
!      kset -- number of valid subsets
!      mini(nnset),maxi(nnset) -- minimum and maximum i bounds of each subset
!      minj(nnset),maxj(nnset) -- minimum and maximum j bounds of each subset
!      newrad -- indicates whether valid initial range subset has been found.
!         (1=no, 0=yes, 2=no valid sets after 10 radius adjustments.)
!      nnset -- maximum number of subsets
!      nset -- number of subsets found
!      nx -- number of DEM cells in x direction
!      ny -- number of DEM cells in y direction
!      rnull -- value for zdem when not in problem domain
!      omini,omaxi,ominj,omaxj -- minimum and maximum bounds of all sets included
!         in search sphere
!      rad -- radius of trial sphere
!      rnull -- real null value used throughout program 
!      set(nnset) -- number of columns in a contiguous set
!      setflag(nnset) -- indicates valid subset (0), or invalid by containing
!         truncated column (1), or adjacent to DEM boundary (2), or empty (-1)
!      subset(nx,ny) -- indicates set membership of DEM columns
!      xdem,ydem -- x and y locations of column lower left node
!      xrad,yrad,zrad -- x,y, and z coordinates of equation for search sphere radius
!      xcen,ycen,zcen -- center of search sphere relative to DEM origin
!      zavg(nnset) -- average distance of surface of columns in subset to sphere center
!      zmid(nx,ny) -- base of slip surface at midpoint of intersected column
!      zb(nx,ny) -- base of slip surface at each column node
!      zdem(nx,ny) -- DEM elevations
!      zz -- square of zrad
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
        USE CommonData
        USE GridData, ONLY: nx,ny,xcen,ycen,zcen,rad,zdem,delxy,halfdelxy,radsq
        USE SetData, ONLY: omini,omaxi,ominj,omaxj,mini,minj,maxi,maxj,set,&
                           setflag,subset,insphere,in,nnset,zb,xdem,ydem,zmid
        USE SearchData, ONLY:  iincflag,allin,includegrid,minclude,nincpt
        USE FosData, ONLY: single
        USE FailSurfData, ONLY: ifailsurf        
        
        IMPLICIT NONE

        INTEGER, INTENT(out) :: nset,kset,ierr
        INTEGER, INTENT(in) :: newrad
        INTEGER :: nn,i,j,m,nodecount
        INTEGER :: j1,j2,i1,i2
        INTEGER :: foundset,zcount(nnset),goodset,incltot
        REAL(pr), INTENT(out) :: zavg(nnset)
        REAL(pr) :: xrad,yrad,zz,zrad

        CHARACTER*70 :: problemtype
        CHARACTER*120 :: errmessage,solution
      
        errmessage = ' '
        solution = ' '
        problemtype = 'determining subsets'   
                
        allin = 0 
        minclude = 0
        incltot = 0
        nset = 0
        kset = 0
        omini = 999999
        omaxi = 0
        ominj = 999999
        omaxj = 0   
        mini = 0
        minj = 0
        maxi = 0
        maxj = 0
        set = 0
        ierr = 0
        nodecount = 0
        zavg = rnull
        goodset = 0
        
!    Calculate bounds so that the starting and ending columns always
!    fall at least two columns outside of the search sphere "shadow"
!    to initialize arrays adequately without initializing full array 
!    (which has dimensions of DEM) each time to save on run time.        
        j1 = INT((ycen-rad)/delxy)-3
        j2 = INT((ycen+rad)/delxy)+3
        i1 = INT((xcen-rad)/delxy)-3
        i2 = INT((xcen+rad)/delxy)+3
        
        j1 = MAX(j1,1)
        j2 = MIN(j2,ny+1)
        i1 = MAX(i1,1)
        i2 = MIN(i2,nx+1)
        
        DO j = j1,j2
          DO i = i1,i2
            IF (i.lt.nx+1.and.j.lt.ny+1) THEN
              subset(i,j) = 0
              in(i,j) = 0
            END IF
!    Determine if DEM node is intersected by search sphere and track if
!    intersection results in truncated surface. Calculate slip base elevation and
!    save min and max bounds of all sets and number of intersected columns.
            CALL checkinsphere(i,j,nodecount)
          END DO
        END DO

!    If no intersected columns and not doing single run, return to Main to
!    recalculate radius.
        IF (nodecount.eq.0) THEN
          IF (single.eq.1) THEN
            errmessage = 'radius for single surface does not intersect grid'
            solution = 'check coordinates for center of sphere'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
          ELSE
            RETURN
          END IF
        END IF
        
!    Expand bounds of all intersected nodes to take column calculations
!    into account. 
        j1 = MAX(ominj-1,1)
        j2 = MIN(omaxj,ny)
        i1 = MAX(omini-1,1)
        i2 = MIN(omaxi,nx)
        
        DO j = j1,j2
          DO i = i1,i2
!    Determine set membership of each column with two or more adjacent intersected nodes 
!    and number of subsets. Track min and max array bounds of columns in each set.
            CALL checksets(i,j,i1,j1,nset,ierr)
!    If too many subsets, return to Main and go to next search node.
            IF (ierr.eq.1) THEN
              IF (single.eq.1) THEN
                errmessage = 'exceeded maximum number of subsets'
                solution = 'decrease DEM resolution or increase minimum size criteria'
                Call WriteError(2,errmessage,problemtype,solution,'no ',0,' ')               
                PRINT *,'         maximum number of subsets ("nnset") = ',nnset
                WRITE (20,*) '        maximum number of subsets ("nnset") = ',nnset
                CLOSE (20)
                STOP
              ELSE
                RETURN
              END IF
            END IF
          END DO
        END DO
        
!    If some nodes were intersected but no columns were intersected return to Main 
!    to recalculate radius.        
        IF (nset.eq.0) RETURN
        
!     Determine min and max bounds of columns in union of all sets.       
        omini = MINVAL(mini,MASK=(mini.gt.0))
        omaxi = MAXVAL(maxi)
        ominj = MINVAL(minj,MASK=(minj.gt.0))
        omaxj = MAXVAL(maxj)
 
!     Check if more sets counted than valid sets from readjusting set
!     membership during search.  Remove empty sets if at end of
!     set line (higher numbered sets with no good sets after them)
!     and count how many good sets.
        nn = nset
        foundset = 0

        DO m = nn,1,-1
          IF (set(m).eq.0) THEN
            IF (m.eq.nn) THEN
              nset = nset - 1
            ELSE
              IF ((set(m+1).eq.0).and.(foundset.eq.0))&
                  nset = nset - 1
            END IF
              setflag(m) = -1  ! empty set
          ELSE
            foundset = 1
            setflag(m) = 0
            kset = kset + 1
          END IF
        END DO
          
!    If no valid sets return to Main to recalculate radius.
        IF (kset.eq.0) THEN
          nset = 0
          IF (single.eq.1) THEN
            errmessage = 'no valid sets for single surface'
            solution = 'check coordinates for center of sphere'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
          ELSE
            RETURN
          END IF
        END IF
          
        zavg = 0.0_pr
        zcount = 0         
!     Determine whether any set is adjacent to DEM boundary (in=-2) or 
!     contains truncated column (in=-1) and flag it.
        DO j = ominj,omaxj
          DO i = omini,omaxi
            m = subset(i,j)
            IF (m.gt.0)  THEN
              IF (in(i,j).eq.-1.and.setflag(m).ne.2) setflag(m) = 1  ! Set contains truncated column.
              IF (in(i,j).eq.-2) setflag(m) = 2  ! Set contains boundary.

!    If inclusion file, check if include columns are contained in set.            
              IF (iincflag.eq.1) THEN
                IF (includegrid(i,j).eq.1) THEN
                  IF (setflag(m).eq.1.and.allin.ne.-2) allin = -1
                  IF (setflag(m).eq.2) allin = -2
                  IF (minclude.eq.0.or.minclude.eq.m) THEN
                    minclude = m
                    IF (in(i,j).eq.4) incltot = incltot + 1
                  ELSE
                    minclude = -1  ! include area is in more than one set
                  END IF
                END IF
              END IF
                
!     Find circular slip base elevation at column center to compare to defined failure surface elevation
!     in volume and FOS calculations.  The higher base elevation is used in calculations.                       
              IF (ifailsurf.eq.1.and.in(i,j).eq.4) THEN        
                xrad = xdem(i,j) + halfdelxy - xcen
                yrad = ydem(i,j) + halfdelxy - ycen
                zz = radsq-xrad*xrad-yrad*yrad
                zrad = SQRT(zz)
                zmid(i,j) = zcen - zrad
              END IF
              IF (newrad.eq.1) THEN 
!     Find average elevation of set columns for radius adjustment. 
                zavg(m) = zdem(i,j)+zavg(m)
                zcount(m) = zcount(m) + 1
              END IF  
            END IF
          END DO
        END DO

!    Determine if any good sets and complete average DEM distance calculation.          
        DO m = 1,nset
          IF (setflag(m).eq.0) THEN 
            goodset = 1
          END IF
          IF (zcount(m).ne.0.and.newrad.eq.1) THEN
            zavg(m) = zavg(m)/REAL(zcount(m),pr)
          END IF
        END DO
        PRINT *, goodset
!    If include area is contained in one set and set does not contain boundary or truncated surface, then good set.         
        IF (iincflag.eq.1.and.incltot.eq.nincpt.and.allin.eq.0) allin = 1
        IF (single.eq.1) THEN
          IF (goodset.eq.0) THEN
            errmessage = 'no valid sets for single surface (contains boundary or truncated surface)'
            solution = 'check coordinates for center of sphere'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')             
          END IF
          IF (minclude.eq.-1) THEN
            errmessage = 'include area is contained by more than one set'
            solution = 'check coordinates for center of sphere and region defined by inclusion area'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
          IF (allin.eq.-1) THEN
            errmessage = 'include area is contained by invalid set (contains truncated surface)'
            solution = 'check coordinates for center of sphere; inclusion area should not extend to DEM boundary'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
          IF (allin.eq.-2) THEN
            errmessage = 'include area is contained by invalid set (contains boundary)'
            solution = 'check coordinates for center of sphere; inclusion area should not extend to DEM boundary'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
          IF (incltot.lt.nincpt) THEN
            errmessage = 'include area is not fully contained'
            solution = 'check coordinates for center of sphere and extent of inclusion area'
            Call WriteError(1,errmessage,problemtype,solution,'no ',0,' ')           
          END IF
        END IF

        END SUBROUTINE Subsets
           
          
