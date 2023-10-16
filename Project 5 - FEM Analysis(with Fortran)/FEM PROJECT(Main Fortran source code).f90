
    
      MODULE FEM1
          INTEGER(4),SAVE :: NPOIN,NELEM,NMATS,NBOUN,NPROP,NNODE,&
                             NEVAB,NSVAB,NDOFN,NDIME,NSTRE,NCASE,NTYPE
      END MODULE FEM1
      
      MODULE FEM2
           Integer(4), Allocatable :: LNODS(:,:),MATNO(:),IFPRE(:),TRACODE(:),POINTCODE(:),&
                                      DISTCODE(:),CONCODE(:)
           Real(8), Allocatable :: PROPS(:,:),COORD(:,:),FIXED(:),& 
                               RLOAD(:,:),ELOAD(:,:),ASTIF(:,:),ASLOAD(:),&
                               XDISP(:),TDISP(:,:),REACT(:),TREAC(:,:),&
                               STRESS(:,:),FACTOR(:,:),SIDENODE(:,:),TRACLOAD(:,:),&
                               POINTCOORD(:,:),POINTFORCE(:,:),&
                               DISTLOAD(:,:),PROXDIST(:),CONCLOAD(:,:),M(:,:) 
     END MODULE FEM2
    
     MODULE FEM3
     INTEGER(4),SAVE ::LOADNO
     REAL(8),ALLOCATABLE::U(:),UD(:),UDD(:),U2(:),U2D(:),U2DD(:),&
         &      KEFF(:,:),LEFF(:),DYNOLOAD(:,:), &  
         &      FRE(:,:),AMP(:,:)           
     REAL(8),SAVE::T,DT,TTOT,TLOAD,DL,ALPHA  
     REAL(8),SAVE::A0,A1,A2,A3,A4,A5,A6,A7
    END MODULE FEM3
    
          
 
    MODULE FEM4
          INTEGER(4),SAVE :: NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,&
                             NEVAB,NSVAB,NDOFN,NDIME,NSTRE,NCAS,Option,AType,thick
      END MODULE FEM4
      
      MODULE FEM5

          Integer,Allocatable :: LNODS(:,:),MATNO(:),IFPRE(:)
          Real,Allocatable :: PROPS(:,:),COORD(:,:),FIXED(:),RLOAD(:,:),ELOAD(:,:),ASLOD(:),XDISP(:),REACT(:),ASTIF(:,:),TDISP(:,:),TREAC(:,:),STRES(:,:),ESTIF(:,:)
          
      END MODULE FEM5
    
    
    
    
       
!     ****************************************************************    
!     PROGRAM FOR TRUSS STRUCTURES ANALISYS WITH DYNAMIC LOAD
!     USING NEWMARK METHOD
!     ****************************************************************
!	  THIS SUBROUTINE CONTROLS THE CALLING,IN ORDER,OF ALL SUBROUTINES
!     ****************************************************************
      PROGRAM ANALYSIS
      
      USE FEM1
      USE FEM2
      USE FEM3
      USE FEM4
      USE FEM5

      
      Print*,"***************** Welcome *********************"
      Print*,"This Program is Developed by SEYED ALI KHODAEI"
      Print*,"                              MEHRAB ZAMANIAN"
      Print*,"                               ARMIN MAHDAVI"
      Print*,"************* Winter 2020 @ UT ****************"
      Print*," "
      Print*," *Please Choose One of Options Below "
      Print*,"     1. Defining My Structure in BlackBox(only available for CST)"
      Print*,"     2. Defining My Structure in an Input File"
      Print*,">> PRESS '0' TO EXIT"
      Read*,Option
     
      IF(OPTION==0)THEN 
          PRINT*," BA-BYE :D "
          PAUSE
      STOP
      END IF
      Print*," *Please Choose Type of Structure"
      Print*," "
      Print*," 1) 1D Truss     2) 2D Truss      3) Frame"
      Print*," 4) 2D Elasticity"
      Print*," 5)Dynamic Analysis"
      READ*,ATYPE
      
      
  IF(ATYPE==5) THEN
              PRINT*,"1) Truss   2) Frame "
              READ*,S
              
              IF(S==1) THEN
 
      OPEN(7,FILE="C:\Users\asus\Desktop\DYN_TRUSS\DYN_TRUSS\IMPORT_DYN_TRUSS.txt",STATUS="OLD",&
     	 FORM="FORMATTED",ACTION="READWRITE")            
      OPEN(8,FILE="C:\Users\asus\Desktop\DYN_TRUSS\DYN_TRUSS\EXPORT_DYN_TRUSS.txt",STATUS="REPLACE",&
     	 FORM="FORMATTED",ACTION="READWRITE")   
      OPEN(1,FILE="TMP.tmp",STATUS="SCRATCH",RECL=256,FORM="UNFORMATTED")
      OPEN(2,FILE="TMP2.tmp",STATUS="SCRATCH",RECL=256,FORM="UNFORMATTED")  

      CALL DATA_DYM_TRUSS
      CALL LOADTRUSS
      CALL STIFFB
      CALL ASSEMBATRUSS
      CALL MASSTRUSS
      CALL MASSASTRUSS
      CALL PARAMETERSFRAME
      DO T=DT,TTOT,DT
        CALL STEP
        IF(T==DT) THEN
             CALL GREDUC
        ELSE
            CALL RESOLVE
        ENDIF
        CALL BAKSUB
        CALL RESULTS
      END DO                                                
      CLOSE (1)
      CLOSE (2)
      CLOSE(8,STATUS="KEEP")                                                          
      CLOSE(7,STATUS="KEEP")  
      PAUSE                                                  
      STOP     
              ELSE IF (S==2) THEN
                       OPEN(7,FILE="C:\Users\asus\Desktop\DYN_FRAME\DYN_FRAME\IMPORT_DYN_FRAME.txt",STATUS="OLD",&
     	 FORM="FORMATTED",ACTION="READWRITE")      
      OPEN(8,FILE="C:\Users\asus\Desktop\DYN_FRAME\DYN_FRAME\EXPORT_DYN_FRAME.txt",STATUS="REPLACE",&
     	 FORM="FORMATTED",ACTION="READWRITE")   
      OPEN(1,FILE="TMP.tmp",STATUS="SCRATCH",RECL=256,FORM="UNFORMATTED")
      OPEN(2,FILE="TMP2.tmp",STATUS="SCRATCH",RECL=256,FORM="UNFORMATTED")  

      CALL DATA_DYM_FRAME
      CALL LOADFRAME
      CALL STIFFC 
      CALL ASSEMBAFRAME
      CALL MASSFRAME
      CALL MASSASFRAME
      CALL PARAMETERSFRAME
      DO T=DT,TTOT,DT
        CALL STEP
        IF(T==DT) THEN
             CALL GREDUC
        ELSE
            CALL RESOLVE
        ENDIF
        CALL BAKSUB
        CALL RESULTS
      END DO                                                
      CLOSE (1)
      CLOSE (2)
      CLOSE(8,STATUS="KEEP")                                                          
      CLOSE(7,STATUS="KEEP")  
      PAUSE                                                  
      STOP    
            End IF           
          END IF     
      
! ***********************************************************************       
! ********************* CST ANALYSIS ALGORITHM   ***********************    
          IF (ATYPE==4) then
          PRINT*,"1)PLANE STRAIN  2)PLANE STRESS"
          READ*,S
     
    IF(S==1) THEN  
    thick=1
    ELSE 
    Print*,"Thickness of Elements"
    Read*,thick
    END IF
     END IF         
      OPEN(7,FILE="E:\First Semester\Finite Element\Task1\IMPORT_CST.txt",STATUS="OLD",&
     	 FORM="FORMATTED",ACTION="READWRITE")    
      OPEN(8,FILE="E:\First Semester\Finite Element\Task1\EXPORT_CST.txt",STATUS="REPLACE",&
     	 FORM="FORMATTED",ACTION="READWRITE")  
      
    
      OPEN(1,FILE="TMP.TMP",STATUS="SCRATCH",RECL=450,FORM="UNFORMATTED")
      Open(2,FILE="TMP2.TMP",STATUS="SCRATCH",RECL=450,FORM="UNFORMATTED")    
CALL DATA_CST      
CALL STIFF_CST
 CALL ASSEMBA_CST
 DO ICASE=1,NCAS
CALL LOAD_ELEM_CST 
CALL ASSEMBB_CST
IF (ICASE==1) THEN
 CALL GREDUC_CST
ELSE
 CALL RESOLVE_CST
END IF
CALL BAKSUB_CST
CALL FORCE_CST                                               
CALL RESULT_CST
END DO
CLOSE (1)                                                         
 CLOSE(8,STATUS="KEEP")                                                          
CLOSE(7,STATUS="KEEP") 
 STOP  
! ***********************************************************************       
! ******************* DYNAMIC ANALYSIS ALGORITHM ***********************       
         
          
          
          
119 END

!     ********************************************************************
!     SUBROUTINE DATA_FRAME:THIS SUBROUTINE IS USED TO GET THE BASIC INFORMATION
!                     OF THE FRAME STRUCTURE FROM FILE NO.7
!     ********************************************************************  
          SUBROUTINE DATA_DYM_FRAME
     
      USE FEM1
      USE FEM2
      USE FEM3
      
      Integer(2), Allocatable :: ICODE(:)
      Real(4), Allocatable :: PRESC(:)
      CHARACTER*20  ::  TITLE   
      READ(7,*) TITLE                                                 
      WRITE(8,*) TITLE                                                         
!
!     READ AND WRITE THE CONTROL DATA                                   
!                                    
      READ(7,*)  NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,NDOFN,NDIME,& 
                 NSTRE,TTOT,DT,TLOAD,DL,ALPHA                                                                                               
      WRITE(8,905) NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,NDOFN,NDIME,&
                   NSTRE,TTOT,DT                                               
905                FORMAT(//1X,'NPOIN =',I5,3X,'NELEM =',I5,3X,'NBOUN =',I5,3X,'NMATS =',I5,3X,&
                               'NPROP =',I5,3X,'NNODE =',I5,3X,'NDOFN =',I5,3X,'NDIME =',I5,3X,&
                               'NSTRE =',I5,3X,'TOTAL TIME=',F6.2,3X,'DELTA T=',F6.2)                   
      
      NSVAB=NPOIN*NDOFN                                                 
      NEVAB=NNODE*NDOFN      
      
      Allocate(LNODS(NELEM,NNODE),MATNO(NELEM),IFPRE(NSVAB),&
        PROPS(NMATS,NPROP),COORD(NPOIN,NDIME),FIXED(NSVAB),&
        RLOAD(NPOIN,NDOFN),ELOAD(NELEM,NEVAB),ASTIF(NSVAB,NSVAB),&
        ASLOAD(NSVAB),XDISP(NSVAB),TDISP(NPOIN,NDOFN),REACT(NSVAB),&
        TREAC(NPOIN,NDOFN),STRESS(NELEM,NSTRE),ICODE(NDOFN),PRESC(NDOFN),&
        FACTOR(NSVAB,NSVAB),U(NSVAB),UD(NSVAB),UDD(NSVAB),U2(NSVAB),&
        U2D(NSVAB),U2DD(NSVAB),KEFF(NSVAB,NSVAB),LEFF(NSVAB),DYNOLOAD(NSVAB,NSVAB))
                                               
!     READ AND WRITE THE MATERIAL PROPERTIES                            
!                                   
      ELOAD(:,:)=0.0
      WRITE(8,950)                                                      
950   FORMAT('MATERIAL PROPERTIES') 
      DO  IMATS=1,NMATS                                               
          READ(7,*)  JMATS,PROPS(JMATS,:)  
          WRITE(8,910) PROPS(JMATS,:)  
910       FORMAT(4F20.5) 
      END DO     
! 
!     READ AND WRITE THE ELEMENT(MEMBER) NODAL                           
!     (JOINT) CONNECTIONS                                          
!                           
      WRITE(8,33)                                                      
33   FORMAT('ELEMENT',5X,'NODES',18X,'MAT') 
      DO  IELEM=1,NELEM                                               
          READ(7,*)  JELEM,LNODS(JELEM,:),MATNO(JELEM)
          WRITE(8,44) JELEM,LNODS(JELEM,:),MATNO(JELEM)        
44       FORMAT(15I10)
      END DO     
!                                                                       
!     READ AND WRITE THE NODAL(JOINT) COORDINATES                       
!                               
      WRITE(8,111)                                                      
111   FORMAT(5X,'NODE',5X,'COORD')
      DO  IPOIN=1,NPOIN                                               
          READ(7,*) JPOIN,COORD(JPOIN,:)                
          WRITE(8,55) JPOIN,COORD(JPOIN,:)             
55       FORMAT(I10,2F15.5)
      END DO    
!                                                                       
!     READ AND WRITE THE BOUNDARY CONDITIONS                            
!     AND STORE IN GLOBAL VECTORS                                       
!                               
      IFPRE(:)=0                                                    
      FIXED(:)=0.0  
      WRITE(8,980)                                                      
980   FORMAT(1X,'RESTRAINED NODES,FIXITY CODE AND PRESCRIBED VALUES')        
      IF (NBOUN/=0) THEN
              DO IBOUN=1,NBOUN                                               
                  READ(7,*) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)        
                  WRITE(8,940) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)      
940               FORMAT(I10,3(I5,1X,F10.5))                                           
                  DO  IDOFN=1,NDOFN                                               
                      INDEX=(NODFX-1)*NDOFN+IDOFN                                       
                      IFPRE(INDEX)=ICODE(IDOFN)                                         
                      FIXED(INDEX)=PRESC(IDOFN) 
                  END DO
              END DO 
      END IF
      RETURN                                                              	 
    END      

    
    
      SUBROUTINE DATA_DYM_TRUSS
     
      USE FEM1
      USE FEM2
      USE FEM3
      
      Integer(2), Allocatable :: ICODE(:)
      Real(4), Allocatable :: PRESC(:)
      CHARACTER*20  ::  TITLE   
      READ(7,*) TITLE                                                 
      WRITE(8,*) TITLE                                                         
!
!     READ AND WRITE THE CONTROL DATA                                   
!                                    
      READ(7,*)  NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,NDOFN,NDIME,& 
                 NSTRE,TTOT,DT,TLOAD,DL,ALPHA                                                                                               
      WRITE(8,905) NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,NDOFN,NDIME,&
                   NSTRE,TTOT,DT                                               
905                FORMAT(//1X,'NPOIN =',I5,3X,'NELEM =',I5,3X,'NBOUN =',I5,3X,'NMATS =',I5,3X,&
                               'NPROP =',I5,3X,'NNODE =',I5,3X,'NDOFN =',I5,3X,'NDIME =',I5,3X,&
                               'NSTRE =',I5,3X,'TOTAL TIME=',F6.2,3X,'DELTA T=',F6.2)                   
      
      NSVAB=NPOIN*NDOFN                                                 
      NEVAB=NNODE*NDOFN      
      
      Allocate(LNODS(NELEM,NNODE),MATNO(NELEM),IFPRE(NSVAB),&
        PROPS(NMATS,NPROP),COORD(NPOIN,NDIME),FIXED(NSVAB),&
        RLOAD(NPOIN,NDOFN),ELOAD(NELEM,NEVAB),ASTIF(NSVAB,NSVAB),&
        ASLOAD(NSVAB),XDISP(NSVAB),TDISP(NPOIN,NDOFN),REACT(NSVAB),&
        TREAC(NPOIN,NDOFN),STRESS(NELEM,NSTRE),ICODE(NDOFN),PRESC(NDOFN),&
        FACTOR(NSVAB,NSVAB),U(NSVAB),UD(NSVAB),UDD(NSVAB),U2(NSVAB),&
        U2D(NSVAB),U2DD(NSVAB),KEFF(NSVAB,NSVAB),LEFF(NSVAB),DYNOLOAD(NSVAB,NSVAB))
                                               
!     READ AND WRITE THE MATERIAL PROPERTIES                            
!                                   
      ELOAD(:,:)=0.0
      WRITE(8,950)                                                      
950   FORMAT('MATERIAL PROPERTIES') 
      DO  IMATS=1,NMATS                                               
          READ(7,*)  JMATS,PROPS(JMATS,:)  
          WRITE(8,910) PROPS(JMATS,:)  
910       FORMAT(4F20.5) 
      END DO     
! 
!     READ AND WRITE THE ELEMENT(MEMBER) NODAL                           
!     (JOINT) CONNECTIONS                                          
!                           
      WRITE(8,33)                                                      
33   FORMAT('ELEMENT',5X,'NODES',18X,'MAT') 
      DO  IELEM=1,NELEM                                               
          READ(7,*)  JELEM,LNODS(JELEM,:),MATNO(JELEM)
          WRITE(8,44) JELEM,LNODS(JELEM,:),MATNO(JELEM)        
44       FORMAT(15I10)
      END DO     
!                                                                       
!     READ AND WRITE THE NODAL(JOINT) COORDINATES                       
!                               
      WRITE(8,111)                                                      
111   FORMAT(5X,'NODE',5X,'COORD')
      DO  IPOIN=1,NPOIN                                               
          READ(7,*) JPOIN,COORD(JPOIN,:)                
          WRITE(8,55) JPOIN,COORD(JPOIN,:)             
55       FORMAT(I10,2F15.5)
      END DO    
!                                                                       
!     READ AND WRITE THE BOUNDARY CONDITIONS                            
!     AND STORE IN GLOBAL VECTORS                                       
!                               
      IFPRE(:)=0                                                    
      FIXED(:)=0.0  
      WRITE(8,980)                                                      
980   FORMAT(1X,'RESTRAINED NODES,FIXITY CODE AND PRESCRIBED VALUES')        
      IF (NBOUN/=0) THEN
              DO IBOUN=1,NBOUN                                               
                  READ(7,*) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)        
                  WRITE(8,940) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)      
940               FORMAT(I10,2(I5,1X,F10.5))                                           
                  DO  IDOFN=1,NDOFN                                               
                      INDEX=(NODFX-1)*NDOFN+IDOFN                                       
                      IFPRE(INDEX)=ICODE(IDOFN)                                         
                      FIXED(INDEX)=PRESC(IDOFN) 
                  END DO
              END DO 
      END IF
      RETURN                                                              	 
    END      
    
        SUBROUTINE LOADFRAME
         USE FEM1
         USE FEM2
         USE FEM3
         ALLOCATE(FRE(NPOIN,NDOFN),AMP(NPOIN,NDOFN))
               WRITE(8,77)                                                      
77   FORMAT(5X,'NODE',5X,'X FREQUENCY',5X,'Y FREQUENCY',3X,'TETHA FREQUENCY',3X,'X AMPLITUDE',3X,'Y AMPLITUDE','TETHA AMPLITUDE')                           
      FRE(:,:)=0.0
      AMP(:,:)=0.0  
      IPOIN=0    
      DO WHILE (IPOIN < NPOIN)
            READ(7,*) IPOIN,FRE(IPOIN,:),AMP(IPOIN,:)               
            WRITE(8,88) IPOIN,FRE(IPOIN,:),AMP(IPOIN,:)   
88         FORMAT(1X,I5,7X,6F13.5)                
      END DO                                  
    RETURN
    END      

!	  *****************************************************************
!     SUBROUTINE LOAD:THIS SUBROUTINE IS USED TO GET DYNAMIC LOAD
!                         OF THE STRUCTURE FROM FILE NO.7
!     *****************************************************************
    
    
    SUBROUTINE LOADTRUSS
         USE FEM1
         USE FEM2
         USE FEM3
         ALLOCATE(FRE(NPOIN,NDOFN),AMP(NPOIN,NDOFN))
               WRITE(8,77)                                                      
77   FORMAT(5X,'NODE',5X,'X FREQUENCY',5X,'Y FREQUENCY',3X,'X AMPLITUDE',3X,'Y AMPLITUDE')                           
      FRE(:,:)=0.0
      AMP(:,:)=0.0  
      IPOIN=0    
      DO WHILE (IPOIN < NPOIN)
            READ(7,*) IPOIN,FRE(IPOIN,:),AMP(IPOIN,:)               
            WRITE(8,88) IPOIN,FRE(IPOIN,:),AMP(IPOIN,:)   ! ATTENTION:THE LAST POINT'S LOAD MUST BE WRITTEN
88         FORMAT(1X,I5,7X,4F13.5)                !RESPECTIVELY:FREQUENCY FOR EACH DEGREE OF FREEDOM,
      END DO                                  !    AMPLITUDE FOR EACH DEGREE OF FREEDOM
    RETURN
    END      
!     ***************************************************************
!     SUOUTINE STIFFB:IT CALCULATES THE LOCAL MATRIX FOR EACH ELEMENT
!                     WHEN NDOFN=2	
!     ***************************************************************	
	  SUBROUTINE STIFFB
     
	  USE FEM1
      USE FEM2
      
      REAL(8),ALLOCATABLE :: ESTIF(:,:)
              ALLOCATE(ESTIF(NEVAB,NEVAB))
!
!     EVALUATION OF MEMBER STIFFNESS MATRICES                           
!     FOR PIN JOINTED PLANE FRAMEWORKS

      REWIND 1
      
      DO  IELEM=1,NELEM   !LOOP OVER EACH ELEMENT                                            
          LPROP=MATNO(IELEM)                                                
          YOUNG=PROPS(LPROP,1)                                              
          XAREA=PROPS(LPROP,2)                                              
          NODE1=LNODS(IELEM,1)                                              
          NODE2=LNODS(IELEM,2)
	      D1=COORD(NODE2,1)-COORD(NODE1,1)                                  
          D2=COORD(NODE2,2)-COORD(NODE1,2)                                  
          ELENG=SQRT(D1*D1+D2*D2)                                     
          SINTH=D2/ELENG                                                    
          COSTH=D1/ELENG                                                    
          FMULT=YOUNG*XAREA/ELENG                                           
          ESTIF(1,1)=FMULT*COSTH*COSTH                                      
          ESTIF(1,2)=FMULT*SINTH*COSTH                                      
          ESTIF(2,1)=FMULT*SINTH*COSTH                                      
          ESTIF(2,2)=FMULT*SINTH*SINTH									  
	      DO  INODE=1,NNODE                                               
                  DO  JNODE=1,NNODE                                               
                      KOUNT=(-1)**INODE*(-1)**JNODE                                     
                      DO  KNODE=1,NNODE                                               
                          DO  LNODE=1,NNODE                                               
                              INDEX=(INODE-1)*NNODE+KNODE                                       
                              JNDEX=(JNODE-1)*NNODE+LNODE                                       
                              ESTIF(INDEX,JNDEX)=KOUNT*ESTIF(KNODE,LNODE) 
                          END DO
                      END DO
                  END DO 
          END DO
          WRITE(1) ESTIF          
      END DO !END LOOP OVER EACH ELEMENT
      RETURN                                                            
    END
     
       
!     *************************************************************
!     SUBROUTINE MASS:IT ASSEMBLES THE GLOBAL MASS MATRIXS
!     ************************************************************* 
    SUBROUTINE STIFFC
      USE FEM1
      USE FEM2
      
      REAL(8)::T(6,6)
      REAL(8),Allocatable :: ESTIF(:,:),K(:,:),KK(:,:),T1(:,:)
      REAL(8):: INERTIA
      NSVAB=NPOIN*NDOFN                                                 
      NEVAB=NNODE*NDOFN      
      Allocate(ESTIF(NEVAB,NEVAB),K(NEVAB,NEVAB),KK(NEVAB,NEVAB),T1(NEVAB,NEVAB))   
      REWIND 1
      DO  IELEM=1,NELEM   !LOOP OVER EACH ELEMENT                                            
          LPROP=MATNO(IELEM)                                                
          YOUNG=PROPS(LPROP,1)                                              
          XAREA=PROPS(LPROP,2)
          INERTIA=PROPS(LPROP,3)
          NODE1=LNODS(IELEM,1)                                              
          NODE2=LNODS(IELEM,2)
	      D1=COORD(NODE2,1)-COORD(NODE1,1)                                  
          D2=COORD(NODE2,2)-COORD(NODE1,2)                                  
          ELENG=SQRT(D1*D1+D2*D2)                                     
          SINTH=D2/ELENG                                                    
          COSTH=D1/ELENG                                                    
          FMULT=YOUNG*XAREA/ELENG 
          MMULT=YOUNG*INERTIA/ELENG
          K=0.0
          T=0.0
          K(1,1)=(FMULT*COSTH**2)+(12*SINTH**2*MMULT/ELENG**2)
          K(1,2)=(FMULT*COSTH*SINTH)+(12*SINTH*COSTH*MMULT/ELENG**2)
          K(1,3)=(-6*MMULT*SINTH/ELENG)
          K(1,4)=-(FMULT*COSTH**2)-(12*MMULT*SINTH/ELENG**2)
          K(1,5)=(12*MMULT*SINTH*COSTH/ELENG**2)-(FMULT*COSTH*SINTH)
          K(1,6)=-6*MMULT*SINTH/ELENG
          K(2,1)=K(1,2)
          K(2,2)=(12*MMULT*COSTH**2/ELENG**2)+(FMULT*SINTH**2)
          K(2,3)=6*MMULT*COSTH/ELENG
          K(2,4)=(12*MMULT*COSTH*SINTH/ELENG**2)-(FMULT*SINTH*COSTH)
          K(2,5)=(-12*MMULT*COSTH**2/ELENG**2)-(FMULT*SINTH**2)
          K(2,6)=6*MMULT*COSTH/ELENG
          K(3,1)=K(1,3)
          K(3,2)=K(2,3)
          K(3,3)=4*MMULT
          K(3,4)=6*MMULT*SINTH/ELENG
          K(3,5)=-6*MMULT*COSTH/ELENG
          K(3,6)=2*MMULT
          K(4,1)=K(1,4)
          K(4,2)=K(2,4)
          K(4,3)=K(3,4)
          K(4,4)=(FMULT*COSTH**2)+(12*MMULT*SINTH**2/ELENG**2)
          K(4,5)=(FMULT*COSTH*SINTH)-(12*MMULT*SINTH*COSTH/ELENG**2)
          K(4,6)=6*MMULT*SINTH/ELENG
          K(5,1)=K(1,5)
          K(5,2)=K(2,5)
          K(5,3)=K(3,5)
          K(5,4)=K(4,5)
          K(5,5)=(12*MMULT*COSTH**2/ELENG**2)+(FMULT*SINTH**2)
          K(5,6)=-6*MMULT*COSTH/ELENG
          K(6,1)=K(1,6)
          K(6,2)=K(2,6)
          K(6,3)=K(3,6)
          K(6,4)=K(4,6)
          K(6,5)=K(5,6)
          K(6,6)=4*MMULT
          ESTIF=K
          WRITE(1) ESTIF  
      END DO 
      RETURN
    END  
          SUBROUTINE MASSFRAME
      USE FEM1
      USE FEM2
      USE FEM3
      REAL(8) MPART1,MPART2,MPART3
      REAL(8),ALLOCATABLE :: EMASS(:,:),RMAT(:,:)
      ALLOCATE(EMASS(NEVAB,NEVAB),RMAT(NEVAB,NEVAB))
      DO IELEM=1,NELEM
          LPROP=MATNO(IELEM)
          XAREA=PROPS(LPROP,2)
          RO=PROPS(LPROP,4)
          NODE1=LNODS(IELEM,1)
          NODE2=LNODS(IELEM,2)
          D1=COORD(NODE2,1)-COORD(NODE1,1)
          D2=COORD(NODE2,2)-COORD(NODE1,2)
          ELENG=SQRT(D1**2+D2**2)
          SINTH=D2/ELENG
          COSTH=D1/ELENG
          EMASS(:,:)=0.0
          RMAT(:,:)=0.0
          
          MPART1=RO*XAREA*ELENG/420
          MPART2=RO*XAREA*ELENG**2/420
          MPART3=RO*XAREA*ELENG**3/420
          
          DO INODE=1,NNODE
              EMASS(3*(INODE-1)+1,3*(INODE-1)+1)=140*MPART1
              EMASS(3*(INODE-1)+2,3*(INODE-1)+2)=156*MPART1
              EMASS(3*(INODE-1)+2,3*(INODE-1)+3)=22*MPART2*(-1)**(IN0DE+1)
              EMASS(3*(INODE-1)+3,3*(INODE-1)+2)=22*MPART2*(-1)**(INODE+1)
              EMASS(3*(INODE-1)+3,3*(INODE-1)+3)=4*MPART3
          END DO
          
          T=1
          DO INODE=1,NNODE
              IND=3*T+1
              EMASS(IND,3*(INODE-1)+1)=70*MPART1
              EMASS(IND+1,3*(INODE-1)+2)=54*MPART1
              EMASS(IND+1,3*(INODE-1)+3)=13*MPART2*(-1)**(IND)
              EMASS(IND+2,3*(INODE-1)+2)=13*MPART2*(-1)**(IND+1)
              EMASS(IND+2,3*(INODE-1)+3)=70*MPART1
              T=0
          END DO
          RMAT(3,3)=1
          RMAT(6,6)=1
          RMAT(1,1)=COSTH
          RMAT(2,2)=COSTH
          RMAT(4,4)=COSTH
          RMAT(5,5)=COSTH
          RMAT(1,2)=SINTH
          RMAT(2,1)=-SINTH
          RMAT(4,5)=SINTH
          RMAT(5,4)=-SINTH
        EMASS=MATMUL(MATMUL(TRANSPOSE(RMAT),EMASS),RMAT)
        WRITE(2) EMASS
        END DO       
      RETURN
      END

      SUBROUTINE MASSTRUSS
      USE FEM1
      USE FEM2
      USE FEM3
      REAL(8) MPART1,MPART2,MPART3
      REAL(8),ALLOCATABLE :: EMASS(:,:),RMAT(:,:)
      ALLOCATE(EMASS(NEVAB,NEVAB),RMAT(NEVAB,NEVAB))
      DO IELEM=1,NELEM
          LPROP=MATNO(IELEM)
          XAREA=PROPS(LPROP,2)
          RO=PROPS(LPROP,4)
          NODE1=LNODS(IELEM,1)
          NODE2=LNODS(IELEM,2)
          D1=COORD(NODE2,1)-COORD(NODE1,1)  !bedst avardane x har eleman
          D2=COORD(NODE2,2)-COORD(NODE1,2)  !bedst avardane  Y har eleman
          ELENG=SQRT(D1**2+D2**2)
          SINTH=D2/ELENG
          COSTH=D1/ELENG
          EMASS(:,:)=0.0
          RMAT(:,:)=0.0
          
          MPART=RO*XAREA*ELENG/6
         ! MPART2=RO*XAREA*ELENG**2/420
         ! MPART3=RO*XAREA*ELENG**3/420
          
          DO INODE=1,NEVAB
              EMASS(INODE,INODE)=2
              RMAT(INODE,INODE)=COSTH
          END DO
          EMASS(1,3)=1
          EMASS(2,4)=1
          EMASS(3,1)=1
          EMASS(4,2)=1
          EMASS=EMASS*MPART
          RMAT(1,2)=SINTH
          RMAT(3,4)=SINTH
          RMAT(2,1)=-SINTH
          RMAT(4,3)=-SINTH
        EMASS=MATMUL(MATMUL(TRANSPOSE(RMAT),EMASS),RMAT)
        WRITE(2) EMASS
        END DO       
      RETURN
    END
          SUBROUTINE ASSEMBAFRAME
      
      USE FEM1
      USE FEM2           
      REAL(8),ALLOCATABLE :: ESTIF(:,:) 
              ALLOCATE(ESTIF(NEVAB,NEVAB))
!     THIS ROUTINE ASSEMBLES THE ELEMENT(MEMBER)                        
!     STIFFNESSES                           
!                                                                       
      REWIND 1                                                                                                       
      ASTIF(:,:)=0.0                                                     
!
!     ASSEMBLE THE STIFFNESS MATRIX
!     
      DO  IELEM=1,NELEM                                       
          
          READ(1) ESTIF
      
          DO  INODE=1,NNODE                                               
              NODEI=LNODS(IELEM,INODE)                                          
              DO  IDOFN=1,NDOFN                                               
                  NROWS=(NODEI-1)*NDOFN+IDOFN                                     
                  NROWE=(INODE-1)*NDOFN+IDOFN                                                         
!                                                                    
!     ASSEMBLE THE ELEMENT STIFFNESS MATRICES                           
!                                                                       
                  DO  JNODE =1,NNODE                                             
                      NODEJ=LNODS(IELEM,JNODE)                                          
                      DO  JDOFN =1,NDOFN                                              
                          NCOLS=(NODEJ-1)*NDOFN+JDOFN                                     
                          NCOLE=(JNODE-1)*NDOFN+JDOFN                                     
                          ASTIF(NROWS,NCOLS)=ASTIF(NROWS,NCOLS)+ESTIF(NROWE,NCOLE) 
                          
                      END DO
                  END DO
              END DO
          END DO
      END DO  
      RETURN          
    END                                                               

!     *************************************************************
!     SUBROUTINE ASSEMB:IT ASSEMBLES THE ELEMENT STIFFNESS MATRICES
!                       TO MAKE THE GLOBAL STIFFNESS MATRIX
!     *************************************************************     
      SUBROUTINE ASSEMBATRUSS
      
      USE FEM1
      USE FEM2           
      
      REAL(8),ALLOCATABLE :: ESTIF(:,:) 
              ALLOCATE(ESTIF(NEVAB,NEVAB))
      
!     THIS ROUTINE ASSEMBLES THE ELEMENT(MEMBER)                        
!     STIFFNESSES                           
!                                                                       
      REWIND 1                                                      
      
                                                        
      ASTIF(:,:)=0.0

                                                            
!
!     ASSEMBLE THE STIFFNESS MATRIX
!
              
      DO  IELEM=1,NELEM                                       
          
          READ(1) ESTIF
      
          DO  INODE=1,NNODE                                               
              NODEI=LNODS(IELEM,INODE)                                          
              DO  IDOFN=1,NDOFN                                               
                  NROWS=(NODEI-1)*NDOFN+IDOFN                                     
                  NROWE=(INODE-1)*NDOFN+IDOFN                                                         
!                                                                    
!     ASSEMBLE THE ELEMENT STIFFNESS MATRICES                           
!                                                                       
                  DO  JNODE =1,NNODE                                             
                      NODEJ=LNODS(IELEM,JNODE)                                          
                      DO  JDOFN =1,NDOFN                                              
                          NCOLS=(NODEJ-1)*NDOFN+JDOFN                                     
                          NCOLE=(JNODE-1)*NDOFN+JDOFN                                     
                          ASTIF(NROWS,NCOLS)=ASTIF(NROWS,NCOLS)+ESTIF(NROWE,NCOLE) 
                          
                      END DO
                  END DO
              END DO
          END DO
      END DO  
      RETURN          
    END                                                               
!*************************************************************
!     SUBROUTINE MASSEMB:IT ASSEMBLES THE ELEMENT MASS MATRICES
!                       TO MAKE THE GLOBAL MASS MATRIX
!     *************************************************************    
        SUBROUTINE MASSASFRAME
        USE FEM1
        USE FEM2
        USE FEM3
        REAL(8),ALLOCATABLE :: EMASS(:,:)
        ALLOCATE(EMASS(NEVAB,NEVAB),M(NSVAB,NSVAB))
        REWIND 2
        M(:,:)=0
        DO IELEM=1,NELEM
            READ(2) EMASS
            DO INODE=1,NNODE
                NODEI=LNODS(IELEM,INODE)
                DO IDOFN=1,NDOFN
                    NROWS=(NODEI-1)*NDOFN+IDOFN
                    NROWE=(INODE-1)*NDOFN+IDOFN
                    DO JNODE =1,NNODE
                        NODEJ=LNODS(IELEM,JNODE)
                        DO JDOFN =1,NDOFN
                            NCOLS=(NODEJ-1)*NDOFN+JDOFN
                            NCOLE=(JNODE-1)*NDOFN+JDOFN
                            M(NROWS,NCOLS)=M(NROWS,NCOLS)+EMASS(NROWE,NCOLE)
                        END DO       
                    END DO
              END DO
          END DO
        END DO
        RETURN
        END

 
!*************************************************************
!     SUBROUTINE MASSEMB:IT ASSEMBLES THE ELEMENT MASS MATRICES
!                       TO MAKE THE GLOBAL MASS MATRIX
!     *************************************************************    
        SUBROUTINE MASSASTRUSS
        USE FEM1
        USE FEM2
        USE FEM3
        REAL(8),ALLOCATABLE :: EMASS(:,:)
        ALLOCATE(EMASS(NEVAB,NEVAB),M(NSVAB,NSVAB))
        REWIND 2
        M(:,:)=0
        DO IELEM=1,NELEM
            READ(2) EMASS
            DO INODE=1,NNODE
                NODEI=LNODS(IELEM,INODE)
                DO IDOFN=1,NDOFN
                    NROWS=(NODEI-1)*NDOFN+IDOFN
                    NROWE=(INODE-1)*NDOFN+IDOFN
                    DO JNODE =1,NNODE
                        NODEJ=LNODS(IELEM,JNODE)
                        DO JDOFN =1,NDOFN
                            NCOLS=(NODEJ-1)*NDOFN+JDOFN
                            NCOLE=(JNODE-1)*NDOFN+JDOFN
                            M(NROWS,NCOLS)=M(NROWS,NCOLS)+EMASS(NROWE,NCOLE)
                        END DO       
                    END DO
              END DO
          END DO
        END DO
        RETURN
    END
    !     *************************************************************
!     SUBROUTINE INITIAL:IT CALCULATS THE INITIAL VALUE AND MAKES 
!             EFFECTIVE STIFFNESS MATRIX
!     ************************************************************* 
      SUBROUTINE PARAMETERSFRAME
        USE FEM1
        USE FEM2
        USE FEM3
        U(:)=0
        UD(:)=0
        UDD(:)=0         
        U2(:)=0
        U2D(:)=0
        U2DD(:)=0
        A0=1/(ALPHA*DT**2)
        A1=DL/(ALPHA*DT)
        A2=1/(ALPHA*DT)
        A3=1/(2*ALPHA)-1
        A4=DL/ALPHA-1
        A5=DT/2*(DL/ALPHA-2)
        A6=DT*(1-DL)
        A7=DL*DT
        DO I=1,NSVAB
            DO J=1,NSVAB
            KEFF(I,J)=ASTIF(I,J)+A0*M(I,J)
            END DO
        END DO
        RETURN
        END

!     *************************************************************
!     SUBROUTINE INITIAL:IT CALCULATS THE INITIAL VALUE AND MAKES 
!             EFFECTIVE STIFFNESS MATRIX
!     ************************************************************* 
      SUBROUTINE PARAMETERSTRUSS
        USE FEM1
        USE FEM2
        USE FEM3
        U(:)=0
        UD(:)=0
        UDD(:)=0         !INITIAL VALUE IN ZERO TIME STEP
        U2(:)=0.0
        U2D(:)=0.0
        U2DD(:)=0.0
        A0=1/(ALPHA*DT**2)
        A1=DL/(ALPHA*DT)
        A2=1/(ALPHA*DT)
        A3=1/(2*ALPHA)-1
        A4=DL/ALPHA-1
        A5=DT/2*(DL/ALPHA-2)
        A6=DT*(1-DL)
        A7=DL*DT
        DO I=1,NSVAB
            DO J=1,NSVAB
            KEFF(I,J)=ASTIF(I,J)+A0*M(I,J)
            END DO
        END DO
        RETURN
        END
!     *************************************************************
!     SUBROUTINE ITERATION:ADD DYNAMIC LOAD IN EACH POINT AND EACH
!      "DOF" TO GLOBALL LOAD MATRIX AND MAKES EFFECTIVE LOAD MATRIX
!     *************************************************************     
        SUBROUTINE STEP
        
        USE FEM1
        USE FEM2
        USE FEM3
        
        REAL(8),ALLOCATABLE :: UX(:)
        REAL(8) LP
        ALLOCATE(UX(NSVAB))
        NSVAB=NPOIN*NDOFN                                                 
        NEVAB=NNODE*NDOFN      
        ASLOAD (:) =0.0
        UX(:)=0.0
        LEFF(:)=0.0
        LP=0 
        IF(TLOAD==1) THEN
        DO  IPOIN=1,NPOIN      
          DO  IDOFN=1,NDOFN                                               
              NROWS=(IPOIN-1)*NDOFN+IDOFN   
              ASLOAD(NROWS)=ASLOAD(NROWS)+AMP(IPOIN,IDOFN)*COS(FRE(IPOIN,IDOFN)*T) !T DAR LOOPE PROGRAM TAARIF MISHAVAD
          END DO
        END DO
        ELSE IF (TLOAD==2) THEN
                DO  IPOIN=1,NPOIN      
          DO  IDOFN=1,NDOFN                                               
              NROWS=(IPOIN-1)*NDOFN+IDOFN   
              ASLOAD(NROWS)=ASLOAD(NROWS)+AMP(IPOIN,IDOFN)*SIN(FRE(IPOIN,IDOFN)*T) !T DAR LOOPE PROGRAM TAARIF MISHAVAD
          END DO
                END DO
                ENDIF
        UX=A0*U+A2*UD+A3*UDD
        LP=0.0
        DO I=1,NSVAB
            DO J=1,NSVAB
                LP=LP+M(I,J)*UX(J) 
            END DO
            LEFF(I)= ASLOAD(I) + LP
            LP=0
        END DO !bordare R eff mohasebe shod va dar LEFF qarar gereft
        RETURN
        END                                                                  
!     *************************************************************
!     SUBROUTINE GREDUC:THIS ROUTINE REDUCES THE GLOBAL STIFFNESS                         
!                       EQUATIONS BY DIRECT GAUSSIAN ELIMINATION  
!     *************************************************************
      SUBROUTINE GREDUC
      
      USE FEM1
      USE FEM2
      USE FEM3
      FACTOR=0.0
      NEQNS=NSVAB                                                       
      DO IEQNS=1,NEQNS                                             
          IF(IFPRE(IEQNS) /= 1) THEN    !CONTORLS THAT IF THE DEGREE OF FREEDOM HAS A PRESCRIBED VALUE
!                                                                       
!     REDUCE EQUATIONS                                                  
!                                                                       
          PIVOT=KEFF(IEQNS,IEQNS)     
              IF(ABS(PIVOT)<1.0E-10) THEN
                  WRITE(8,900) PIVOT,IEQNS                                          
900               FORMAT(5X,'INCORRECT PIVOT = ',E20.6,5X,'EQUATION NO. ',I5)     
                  STOP
              END IF  
              IF(IEQNS<NEQNS) THEN
                  IEQN1=IEQNS+1                                                     
                  DO IROWS=IEQN1,NEQNS                                           
                      FACTR=KEFF(IROWS,IEQNS)/PIVOT 
                      IF(FACTR/=0.0) THEN
                          DO ICOLS=IEQNS,NEQNS                                  
                              KEFF(IROWS,ICOLS)=KEFF(IROWS,ICOLS)-FACTR*KEFF(IEQNS,ICOLS)
                          END DO
                      END IF    
                      LEFF(IROWS)=LEFF(IROWS)-FACTR*LEFF(IEQNS)
                      FACTOR(IROWS,IEQNS)=FACTR
                  END DO
              END IF
!
!     ADJUST RHS(LOADS) FOR PRESCRIBED DISPLACEMENTS
!
          ELSE
              DO IROWS=IEQNS,NEQNS
                  LEFF(IROWS)=LEFF(IROWS)-KEFF(IROWS,IEQNS)*FIXED(IEQNS)
                  KEFF(IROWS,IEQNS)=0.0
              END DO
          END IF
      END DO          
      RETURN 
      END
!     *************************************************************
!     SUBROUTINE RESOLVE:THIS ROUTINE REDUCES THE LOAD MATRICES                        
!                        BY DIRECT GAUSSIAN ELIMINATION  
!     *************************************************************
      SUBROUTINE RESOLVE  
      
      USE FEM1
      USE FEM2
      USE FEM3
      NEQNS=NSVAB
      DO IEQNS=1,NEQNS                                             
          IF(IFPRE(IEQNS) /= 1) THEN                                                              
              IF(IEQNS<NEQNS) THEN
                  IEQN1=IEQNS+1                                                     
                  DO IROWS=IEQN1,NEQNS 
                      FACTR=FACTOR(IROWS,IEQNS)
                      LEFF(IROWS)=LEFF(IROWS)-FACTR*LEFF(IEQNS)
                  END DO
              END IF
!    
!      ADJUST RHS(LOADS) FOR PRESCRIBED DISPLACEMENTS
!     
          ELSE
              DO IROWS=IEQNS,NEQNS
                  LEFF(IROWS)=LEFF(IROWS)-KEFF(IROWS,IEQNS)*FIXED(IEQNS)
              END DO
          END IF
      END DO    
      RETURN 
      END
!     **********************************************************************
!     SUBROUTINE BAKSUB:THIS ROUTINE PERFORMS THE BACK SUBSTITUTION PHASE                                                                         
!     **********************************************************************      
	  SUBROUTINE BAKSUB
      
      USE FEM1
      USE FEM2
      USE FEM3
       
      NEQNS=NSVAB
      NEQN1=NEQNS+1 
      DO  IEQNS=1,NEQNS                                               
          NBACK=NEQN1-IEQNS                                                 
          PIVOT=KEFF(NBACK,NBACK)                                          
          RESID=LEFF(NBACK)                                             
          IF(NBACK/=NEQNS) THEN                                     
              NBAC1=NBACK+1                                                     
              DO  ICOLS=NBAC1,NEQNS                                           
                  RESID=RESID-KEFF(NBACK,ICOLS)*U(ICOLS)                       
              END DO
          END IF
          IF(IFPRE(NBACK)==0) U2(NBACK)=RESID/PIVOT         
          IF(IFPRE(NBACK)==1) U2(NBACK)=FIXED(NBACK)
      END DO
      DO I=1,NSVAB
        U2DD(I)=A0*(U2(I)-U(I))-A2*UD(I)-A3*UDD(I)
        U2D(I)=UD(I)+A6*UDD(I)+A7*U2DD(I)
      END DO
      RETURN                                                            
      END    
!     *******************************************************************
!     SUBROUTINE RESULT:WRITE THE NODAL(JOINT) DISPLACEMENTS AND REACTIONS
!                   IN EACH TIME STEP AND EACH NODE
!     *******************************************************************
      SUBROUTINE RESULTS
      
      USE FEM1
      USE FEM2
      USE FEM3
      IF (T==DT) THEN
        WRITE(8,10)
10      FORMAT("TIME",3X,"NUMBER OF DOF:",3X,"DISPLACEMENT(U(X))")
        ENDIF
        DO I=1,NSVAB
            IF (IFPRE(I).EQ.0) THEN
            WRITE(8,11) T,I,U2(I)
11          FORMAT(F6.2,'|',3X,I3,10X,'|',E15.4)           
            U(I)=U2(I)
            UD(I)=U2D(I)
            UDD(I)=U2DD(I)
            END IF
        END DO
        WRITE(8,13)
13      FORMAT(6X,'|',16X,'|')    
        RETURN
    END
    
    
    
    
    
    
    ! CST
    
      

      
!     ********************************************************************
!     SUBROUTINE DATA:THIS SUBROUTINE IS USED TO GET THE BASIC INFORMATION
!                     OF THE STRUCTURE FROM FILE NO.7
!     ********************************************************************
      SUBROUTINE DATA_CST
     
      USE FEM4
      USE FEM5
      
      Integer,Allocatable ::  ICODE(:)
      Real,Allocatable    ::  PRESC(:)
      CHARACTER(5)  ::  TITLE
!	                                                               
!     DATA INPUT SUBROUTINE
!     READ AND WRITE THE PROBLEM TITLE  

      READ(7,*) TITLE                                                 
      WRITE(8,*) TITLE                                                    
      
!
!     READ AND WRITE THE CONTROL DATA   
      
     Select Case (Option)
         
     Case(1)
     Print*,"======================================================================"
     Print*,"                        Introducing Structure"
     Print*,"======================================================================"

     
     
     Print*,"Number Of Global Nodes: "; Read*,NPOIN
     Print*,"Number Of Elements: "; Read*,NELEM
     Print*,"Number Of Points With Boundary Conditions: "; Read*,NBOUN
     Print*,"Number Of Elements with Exclusive Properties: "; Read*,NMATS
     Print*,"Number Of Properties (E/V/I) : "; Read*,NPROP
     Print*,"Number Of Local Nodes: "; Read*,NNODE
     Print*,"Degrees of Freedom in Each Point: "; Read*,NDOFN
     Print*,"Dimension of Structure: "; Read*,NDIME
     Print*,"Number Of Force Components in Each Point: "; Read*,NSTRE   
     Print*,"Number Of Load Cases: ";Read*,NCAS

     
   
     Case(2)
     READ(7,*)  NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,NDOFN,NDIME,& 
                 NSTRE,NCAS


     End Select
     
     
     
!  900 FORMAT(16I5)                                                  
      WRITE(8,905) NPOIN,NELEM,NBOUN,NMATS,NPROP,NNODE,NDOFN,NDIME,&
                   NSTRE,NCAS                                              
905   FORMAT(//1X,'NPOIN =',I5,3X,'NELEM =',I5,3X,'NBOUN =',I5,3X,&
              'NMATS =',I5,//1X,'NPROP =',I5,3X,'NNODE =',I5,3X,&
               'NDOFN =',I5,3X,'NDIME =',I5,//1X,'NSTRE =',I5,3X,'NCASE =',I5)                   
      
      NSVAB=NPOIN*NDOFN                                                 
      NEVAB=NNODE*NDOFN      
      Allocate(LNODS(NELEM,NNODE),MATNO(NELEM),IFPRE(NSVAB),ICODE(NSVAB),ASTIF(NSVAB,NSVAB),STRES(NELEM,NSTRE),TDISP(NPOIN,NDOFN),TREAC(NPOIN,NDOFN),ESTIF(NEVAB,NEVAB))
      Allocate(PROPS(NMATS,NPROP),COORD(NPOIN,NDIME),FIXED(NSVAB),RLOAD(NPOIN,NDOFN),ASLOD(NSVAB),XDISP(NSVAB),REACT(NSVAB),PRESC(NSVAB),ELOAD(NELEM,NEVAB))
  
      
      
!                                    -+-                                
!     READ AND WRITE THE MATERIAL PROPERTIES                            
!                                   

      WRITE(8,950)                                                      
950   FORMAT('MATERIAL PROPERTIES') 
      DO  IMATS=1,NMATS   
          
          IF(Option==1) Then
          Print*,"Enter 1.Material's Type Label  2.E (Young Model) / 3.A (Area) / 4.I (Moment of Inertia)"
          READ(*,*)  JMATS,PROPS(JMATS,:)             

          Else
          READ(7,*)  JMATS,PROPS(JMATS,:)  
          End IF
          
          WRITE(8,910) PROPS(JMATS,:)  
910       FORMAT(4F20.5) 
      END DO
      
! 
!     READ AND WRITE THE ELEMENT(MEMBER) NODAL                          
!     (JOINT) CONNECTIONS                                               
!                            

      WRITE(8,960)                                                      
960   FORMAT('0',2X,'ELEMENT',3X,'NODES',3X,'MAT.') 
      DO  IELEM=1,NELEM               
          
          IF(Option==1) Then
          Print*,"Enter 1.Element's Label  2.Introduce Nodes that This Element is Consists of 3.Material's Type Number"
          READ(*,*) JELEM,LNODS(JELEM,:),MATNO(JELEM)
          
          Else
          READ(7,*)  JELEM,LNODS(JELEM,:),MATNO(JELEM)  
          End IF
          
          WRITE(8,920) JELEM,LNODS(JELEM,:),MATNO(JELEM)
920       FORMAT(5I5)
      END DO
 
      
 
      
      
      
!                                                                       
!     READ AND WRITE THE NODAL(JOINT) COORDINATES                       
!                               

      WRITE(8,970)                                                      
970   FORMAT(5X,'NODE',5X,'COORD')
      
      DO  IPOIN=1,NPOIN          
          IF(Option==1) Then
                        
          Print*,"Enter 1.Node's Label  2.Coordinations 'X' and 'Y' "
          Read(*,*) JPOIN,COORD(JPOIN,:)   
          
          Else
          READ(7,*) JPOIN,COORD(JPOIN,:)     
          End IF
          
          WRITE(8,930) JPOIN,COORD(JPOIN,:)             
930       FORMAT(I10,2F15.5)
      END DO
      
!                                                                       
!     READ AND WRITE THE BOUNDARY CONDITIONS                            
!     AND STORE IN GLOBAL VECTORS                                       
!                               

      IFPRE(:)=0                                                    
      FIXED(:)=0.0
      
      WRITE(8,980)                                                      
980   FORMAT(1X,'RESTRAINED NODES,FIXITY CODE AND PRESCRIBED VALUES')        
      IF (NBOUN/=0) THEN
              DO IBOUN=1,NBOUN   
                  
                  
          IF(Option==1) Then              
          Print*,"Enter 1.Label for Intended fixed Node  2. Icode  3. Presc  4.Icode  5. Presc"
          Read(*,*) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)        
          Else
                   READ(7,*) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)        
          End IF
          
                  WRITE(8,940) NODFX,(ICODE(IDOFN),PRESC(IDOFN),IDOFN=1,NDOFN)      
940               FORMAT(I10,2(I5,X,F10.5))                                           
                  DO  IDOFN=1,NDOFN                                               
                      INDEX=(NODFX-1)*NDOFN+IDOFN                                       
                      IFPRE(INDEX)=ICODE(IDOFN)                                         
                      FIXED(INDEX)=PRESC(IDOFN) 
                  END DO
              END DO 
      END IF
      RETURN                                                            	 
    END

    
!   ****************************************************                                                               
!     READ AND WRITE THE NODAL(JOINT) APPLIED LOADS                     
!   ****************************************************   
    
    
    
!   ****************************************************                                                               
!     READ AND WRITE THE ELEMENTAL APPLIED LOADS                     
!   **************************************************** 
    
!	  *****************************************************************
!     SUBROUTINE STIFFA:IT CALCULATES THE LOCAL MATRIX FOR EACH ELEMENT
!                       WHEN NDOFN=1
!     *****************************************************************
      Subroutine LOAD_ELEM_CST

      USE FEM4
      USE FEM5
      
      REAL :: ELEMLOAD(NELEM,NEVAB),Body_Force(NELEM,NEVAB)
      
      ELOAD(:,:)=0.0                
      Body_Force(:,:)=0.0

      Do IELEM=1,NELEM
          
          ELEMLOAD(:,:)=0.0
          
     IF(Option==1) Then
              
     Print*,"======================================================================"
     Print*,"                        Introducing Loading"
     Print*,"======================================================================"

     
     
    ! ******************* Defining and Calculating Body Loads *******************
     
          PRINT*,"PLEASE Enter 1. Fx  2.Fy for Body Force"
          Read(*,*) fx,fy
          
     Else
                   Read(7,*) fx,fy
     End IF
     
            Body_Force(IELEM,1) = (thick*ELAREA/3)*fx
            Body_Force(IELEM,2) = (thick*ELAREA/3)*fy
            Body_Force(IELEM,3) = (thick*ELAREA/3)*fx
            Body_Force(IELEM,4) = (thick*ELAREA/3)*fy
            Body_Force(IELEM,5) = (thick*ELAREA/3)*fx
            Body_Force(IELEM,6) = (thick*ELAREA/3)*fy
     
     ELOAD(IELEM,:)=Body_Force(IELEM,:)
     
    
              
            ! ******************* Defining and Calculatin Traction Loads + ELOAD *******************
  
          IF(Option==1) Then
            
          Print*,"1.Node1  2.Node2  3.Node3"
          Print*,"*** Traction is on Node1-Node2 Side"
          READ(*,*) NODE1,NODE2,NODE3
          
          Print*,"1.Tx1 2.Tx2 3.Ty1 4.Ty2"
          
          READ(*,*) Tx1,Tx2,Ty1,Ty2
          
          ELSE
              
           READ(7,*) NODE1,NODE2,NODE3
           READ(7,*) Tx1,Tx2,Ty1,Ty2
           
          END IF          
          
    
       D1=COORD(NODE2,1)-COORD(NODE1,1)                                  
       D2=COORD(NODE2,2)-COORD(NODE1,2) 
     
          ELENG=SQRT(D1*D1+D2*D2)
          VALUE1=(thick*ELENG/6)
          

IF(NODE1.lt.NODE2 .and. NODE1.lt.NODE3) NODE1=1
IF(NODE2.lt.NODE1 .and. NODE2.lt.NODE3) NODE2=1
IF(NODE3.lt.NODE2 .and. NODE3.lt.NODE1) NODE3=1
IF(NODE1.gt.NODE2 .and. NODE1.lt.NODE3) NODE1=2
IF(NODE1.gt.NODE3 .and. NODE1.lt.NODE2) NODE1=2
IF(NODE2.gt.NODE1 .and. NODE2.lt.NODE3) NODE2=2
IF(NODE2.gt.NODE3 .and. NODE2.lt.NODE1) NODE2=2
IF(NODE3.gt.NODE2 .and. NODE1.lt.NODE1) NODE3=2
IF(NODE3.gt.NODE1 .and. NODE1.lt.NODE2) NODE3=2
IF(NODE1.gt.NODE2 .and. NODE1.gt.NODE3) NODE1=3
IF(NODE2.gt.NODE3 .and. NODE2.gt.NODE1) NODE2=3
IF(NODE3.gt.NODE2 .and. NODE3.gt.NODE1) NODE3=3

          
          
ELOAD(IELEM,2*(NODE1)-1)= VALUE1*(2*Tx1 + Tx2)
ELOAD(IELEM,2*(NODE1)) = VALUE1*(2*Ty1 + Ty2)
ELOAD(IELEM,2*(NODE2)-1) = VALUE1*(Tx1 + 2*Tx2)
ELOAD(IELEM,2*(NODE2)) = VALUE1*(Ty1 + 2*Ty2)     
ELOAD(IELEM,2*(NODE3)-1)=0
ELOAD(IELEM,2*(NODE3))=0

          ELOAD = ELEMLOAD + ELOAD
          
             ELEMLOAD = ELOAD

             
             
          
       
                        WRITE(8,*) IELEM,ELOAD(IELEM,:) 

      END DO
      
      
         WRITE(8,990)                                                      
990   FORMAT(5X,'NODE',7X,'LOADS') 
                                 
      RLOAD(:,:)=0.0
          
      IPOIN=0    
      DO WHILE (IPOIN < NPOIN)
          IF(Option==1) Then
     Print*,"======================================================================"
     Print*,"                        Introducing Loading"
     Print*,"======================================================================"
                 Print*,"Enter 1.Number of Node  2.Fx   3.Fy "

            Read(*,*) IPOIN,RLOAD(IPOIN,:) 
            Else
            READ(7,*) IPOIN,RLOAD(IPOIN,:)     
            End IF
            WRITE(8,930) IPOIN,RLOAD(IPOIN,:)     
            
930         FORMAT(I10,2F15.5)
            
      END DO     
     
      
Return
    End
    
      
      
!     ***************************************************************
!     SUOUTINE STIFFB:IT CALCULATES THE LOCAL MATRIX FOR EACH ELEMENT
!                     WHEN NDOFN=2	
!     ***************************************************************	

    
    
    	SUBROUTINE STIFF_CST                                                         
     
	 USE FEM4
      USE FEM5
                Real :: VALUE1,VALUE2,VALUE3,VALUE4

      Real :: D(3,3),B(6,3),MULTI(6,3)
      REWIND 1
      


      DO  IELEM=1,NELEM   !LOOP OVER EACH ELEMENT       
          
          LPROP=MATNO(IELEM)                                                
          YOUNG=PROPS(LPROP,1)     !ELASTICITY MODULE                                           
          V=PROPS(LPROP,2)         !POISSON RATIO
          
          NODE1=LNODS(IELEM,1)                                              
          NODE2=LNODS(IELEM,2)
          NODE3=LNODS(IELEM,3)

          BETTA1=COORD(NODE2,2)-COORD(NODE3,2)
          BETTA2=COORD(NODE3,2)-COORD(NODE1,2)
          BETTA3=COORD(NODE1,2)-COORD(NODE2,2)

          GAMMA1=COORD(NODE3,1)-COORD(NODE2,1)
          GAMMA2=COORD(NODE1,1)-COORD(NODE3,1)
          GAMMA3=COORD(NODE2,1)-COORD(NODE1,1)
        
          
          
          ! ELEMENT 1
	      D1=COORD(NODE2,1)-COORD(NODE1,1)                                  
          D2=COORD(NODE2,2)-COORD(NODE1,2) 
          ! ELEMENT 2          
          D3=COORD(NODE3,1)-COORD(NODE2,1)
          D4=COORD(NODE3,2)-COORD(NODE2,2)
          ! ELEMENT 3
          D5=COORD(NODE3,1)-COORD(NODE1,1)
          D6=COORD(NODE3,2)-COORD(NODE1,2)
          
          ELENG1=SQRT(D1*D1+D2*D2)
          ELENG2=SQRT(D3*D3+D4*D4)
          ELENG3=SQRT(D5*D5+D6*D6)
                    
          SEMIPERIMETER=(ELENG1+ELENG2+ELENG3)/2
          ELAREA=SQRT(SEMIPERIMETER*(SEMIPERIMETER-ELENG3)*(SEMIPERIMETER-ELENG2)*(SEMIPERIMETER-ELENG1))
          
          VALUE1=YOUNG/(1-V*V)
          VALUE2=YOUNG/((1+V)*(1-2*V))
          VALUE3=(2*ELAREA)
          VALUE4=1/VALUE3
          
          If (S==1) Then !(It's Plane Strain)
          ! MATRIX D PLANE STRAIN
          D(1,1)=VALUE2*(1-V)
          D(1,2)=VALUE2*V
          D(1,3)=0
          D(2,1)=VALUE2*V
          D(2,2)=VALUE2*(1-V)
          D(2,3)=0
          D(3,1)=0
          D(3,2)=0
          D(3,3)=VALUE2*(1-2*V)/2
          
          Else !(It's Plane Stress)
          ! MATRIX D PLANE STRESS
          D(1,1)=VALUE1
          D(1,2)=VALUE1*V
          D(1,3)=0
          D(2,1)=VALUE1*V
          D(2,2)=VALUE1
          D(2,3)=0
          D(3,1)=0
          D(3,2)=0
          D(3,3)=VALUE1*(1-V)/2
          
          End If
          
          
          !MATRIX B(TRANSPOSE)
          B(1,1)=VALUE4*BETTA1
          B(2,1)=0
          B(3,1)=VALUE4*BETTA2
          B(4,1)=0
          B(5,1)=VALUE4*BETTA3
          B(6,1)=0
        
          B(1,2)=0
          B(2,2)=VALUE4*GAMMA1
          B(3,2)=0
          B(4,2)=VALUE4*GAMMA2
          B(5,2)=0
          B(6,2)=VALUE4*GAMMA3
        
          B(1,3)=VALUE4*GAMMA1
          B(2,3)=VALUE4*BETTA1
          B(3,3)=VALUE4*GAMMA2
          B(4,3)=VALUE4*BETTA2
          B(5,3)=VALUE4*GAMMA3
          B(6,3)=VALUE4*BETTA3
          
          
          
          
          MULTI=MATMUL(B,D)
          ESTIF=thick*ELAREA*MATMUL(MULTI,TRANSPOSE(B))
      

          WRITE(1) ESTIF       
          
      END DO !END LOOP OVER EACH ELEMENT
      RETURN                                                            
    END
    
!     *************************************************************
!     SUBROUTINE ASSEMBA:IT ASSEMBLES THE ELEMENT STIFFNESS MATRICES
!                       TO MAKE THE GLOBAL STIFFNESS MATRIX
!     *************************************************************
       SUBROUTINE ASSEMBA_CST                                  
      
      USE FEM4
      USE FEM5           
      
      
!     THIS ROUTINE ASSEMBLES THE ELEMENT(MEMBER)                        
!     STIFFNESSES AND APPLIED LOADS TO FORM THE                         
!     GLOBAL STIFFNESS MATRIX AND FORCE VECTOR                          
!                                                                       
      REWIND 1                                                      
                                                     
      ASTIF(:,:)=0.0
                                               

!
!     ASSEMBLE THE STIFFNESS MATRIX AND IF THERE IS ANY ELEMENTAL LOAD
!
              
      DO  IELEM=1,NELEM      
          
          READ(1) ESTIF
      
          DO  INODE=1,NNODE                                               
              NODEI=LNODS(IELEM,INODE)                                          
              DO  IDOFN=1,NDOFN                                               
                  NROWS=(NODEI-1)*NDOFN+IDOFN                                     
                  NROWE=(INODE-1)*NDOFN+IDOFN                                               
!                                                                    
!     ASSEMBLE THE ELEMENT STIFFNESS MATRICES                           
!                                                                       
                  DO  JNODE =1,NNODE                                             
                      NODEJ=LNODS(IELEM,JNODE)                                          
                      DO  JDOFN =1,NDOFN                                              
                          NCOLS=(NODEJ-1)*NDOFN+JDOFN                                     
                          NCOLE=(JNODE-1)*NDOFN+JDOFN                                     
                          ASTIF(NROWS,NCOLS)=ASTIF(NROWS,NCOLS)+ESTIF(NROWE,NCOLE) 
                          
                      END DO
                  END DO
              END DO
          END DO
      END DO     
      
      RETURN          
    END                                                                    
!     ********************************************************
!     SUBROUTINE ASSEMB :THIS ROUTINE ASSEMBLES APPLIED LOADS                        
!                        TO FORM THE FORCE VECTOR                           
!     ********************************************************
          SUBROUTINE ASSEMBB_CST
      USE FEM4
      USE FEM5           
      
                                                                    
      REWIND 1                                                      
      
      ASLOD(:)=0.0                                                    
                                 
!                                                                       
!     ASSEMBLE THE ELEMENT LOADS                                        
!                                                                       
      DO  IPOIN=1,NPOIN                                               
          DO  IDOFN=1,NDOFN                                               
              NROWS=(IPOIN-1)*NDOFN+IDOFN                                       
              ASLOD(NROWS)=ASLOD(NROWS)+RLOAD(IPOIN,IDOFN)
                            
          END DO
      END DO
!
!      ELEMENTAL LOAD
!
              
      DO  IELEM=1,NELEM                                       
          
      
          DO  INODE=1,NNODE                                               
              NODEI=LNODS(IELEM,INODE)                                          
              DO  IDOFN=1,NDOFN                                               
                  NROWS=(NODEI-1)*NDOFN+IDOFN                                     
                  NROWE=(INODE-1)*NDOFN+IDOFN                                     
                  ASLOD(NROWS)=ASLOD(NROWS)+ELOAD(IELEM,NROWE)
                  
              END DO
          END DO
      END DO    
      RETURN          
    END                                
    
!     *************************************************************
!     SUBROUTINE GREDUC:THIS ROUTINE REDUCES THE GLOBAL STIFFNESS                         
!                EQUATIONS BY DIRECT GAUSSIAN ELIMINATION  
!     *************************************************************
          
      SUBROUTINE GREDUC_CST
      
      USE FEM4
      USE FEM5     
      
      NEQNS=NSVAB                                                       
      DO IEQNS=1,NEQNS                                             
          IF(IFPRE(IEQNS) /= 1) THEN    !CONTORLS THAT IF THE DEGREE OF FREEDOM HAS A PRESCRIBED VALUE
!                                                                       
!     REDUCE EQUATIONS                                                  
!                                                                       
          PIVOT=ASTIF(IEQNS,IEQNS)     
          
              IF(ABS(PIVOT)<1.0E-10) THEN
                  WRITE(8,900) PIVOT,IEQNS                                          
900               FORMAT(5X,'INCORRECT PIVOT = ',E20.6,5X,'EQUATION NO. ',I5)     
                  STOP
              END IF
              
              IF(IEQNS<NEQNS) THEN
                  IEQN1=IEQNS+1                                                     
                  DO IROWS=IEQN1,NEQNS                                           
                      FACTR=ASTIF(IROWS,IEQNS)/PIVOT
              WRITE(2) FACTR
                      IF(FACTR/=0.0) THEN
                          DO ICOLS=IEQNS,NEQNS                                  
                              ASTIF(IROWS,ICOLS)=ASTIF(IROWS,ICOLS)-FACTR*ASTIF(IEQNS,ICOLS)
                          END DO
                      END IF    
                      ASLOD(IROWS)=ASLOD(IROWS)-FACTR*ASLOD(IEQNS)
                  END DO
              END IF
!
!     ADJUST RHS(LOADS) FOR PRESCRIBED DISPLACEMENTS
!
          ELSE
              DO IROWS=IEQNS,NEQNS
                  ASLOD(IROWS)=ASLOD(IROWS)-ASTIF(IROWS,IEQNS)*FIXED(IEQNS)
                  ASTIF(IROWS,IEQNS)=0.0
              END DO
          END IF
      END DO     
            
      RETURN 
    END
    
     SUBROUTINE RESOLVE_CST 
      
      USE FEM4
      USE FEM5   
      
     REWIND 2
     
      NEQNS=NSVAB                                                       
      DO IEQNS=1,NEQNS                                             
          IF(IFPRE(IEQNS) /= 1) THEN    !CONTORLS THAT IF THE DEGREE OF FREEDOM HAS A PRESCRIBED VALUE
!                                                                       
!     REDUCE EQUATIONS                                                  
!                                                                       
          PIVOT=ASTIF(IEQNS,IEQNS)     
          
              IF(ABS(PIVOT)<1.0E-10) THEN
                  WRITE(8,900) PIVOT,IEQNS                                          
900               FORMAT(5X,'INCORRECT PIVOT = ',E20.6,5X,'EQUATION NO. ',I5)     
                  STOP
              END IF
              
              IF(IEQNS<NEQNS) THEN
                  IEQN1=IEQNS+1                                                     
                  DO IROWS=IEQN1,NEQNS          
                      
              READ(2) FACTR
                      ASLOD(IROWS)=ASLOD(IROWS)-FACTR*ASLOD(IEQNS)
                      
                  END DO
              END IF
!
!     ADJUST RHS(LOADS) FOR PRESCRIBED DISPLACEMENTS
!
          ELSE
              DO IROWS=IEQNS,NEQNS
                  ASLOD(IROWS)=ASLOD(IROWS)-ASTIF(IROWS,IEQNS)*FIXED(IEQNS)
              END DO
          END IF
      END DO     
            
      RETURN 
      END
!     **********************************************************************
!     SUBROUTINE BAKSUB:THIS ROUTINE PERFORMS THE BACK-                                   
!     SUBSTITUTION PHASE                                                                         
!     **********************************************************************      
	  SUBROUTINE BAKSUB_CST
      
      USE FEM4
      USE FEM5                                                                        
       
      NEQNS=NSVAB                                                       

      REACT(:)=0.0                                                  

      
      
      NEQN1=NEQNS+1 
      DO  IEQNS=1,NEQNS                                               
          NBACK=NEQN1-IEQNS                                                 
          PIVOT=ASTIF(NBACK,NBACK)                                          
          RESID=ASLOD(NBACK)                                             
          IF(NBACK/=NEQNS) THEN                                     
              NBAC1=NBACK+1                                                     
              DO  ICOLS=NBAC1,NEQNS                                           
                  RESID=RESID-ASTIF(NBACK,ICOLS)*XDISP(ICOLS)                       
              END DO
          END IF
          IF(IFPRE(NBACK)==0) XDISP(NBACK)=RESID/PIVOT         !!
          IF(IFPRE(NBACK)==1) XDISP(NBACK)=FIXED(NBACK)                   
          IF(IFPRE(NBACK)==1) REACT(NBACK)=-RESID
      END DO
!
!     PUT THE RESULT DATA IN TDISP AND TREAC IN ORDER TO SHOW TO USER
!      
      KOUNT=0
      DO  IPOIN=1,NPOIN
          DO  IDOFN=1,NDOFN                                               
              KOUNT=KOUNT+1                                                     
              TDISP(IPOIN,IDOFN)=XDISP(KOUNT)                                   
              TREAC(IPOIN,IDOFN)=REACT(KOUNT)
          END DO
      END DO
      RETURN                                                            
      END 

!    **************************************************************
!     SUBROUTINE FORCE:MEMBER FORCE CALCULATIONS FOR BOTH AXIAL                          
!                      BAR AND PIN-JOINTED PLANE FRAME PROBLEMS   
!     **************************************************************
      SUBROUTINE FORCE_CST                               
                                    
      USE FEM4
      USE FEM5 
      REAL(8) :: FOMEM(8)
      
      REWIND 1                                                          
      DO  IELEM=1,NELEM ! LOOP OVER EAC ELEMENT
          READ(1) ESTIF
!                                                                       
!     EVALUATE THE MEMBER END FORCES                                    
!                                                                       
	      DO  IEVAB=1,NEVAB                                               
              FOMEM(IEVAB)=0.0                                                  
              KOUNT=0                                                           
              DO  INODE=1,NNODE                                               
                  LOCAL=LNODS(IELEM,INODE)                                          
                  DO  IDOFN=1,NDOFN                                               
                      KOUNT=KOUNT+1                                                     
                      FOMEM(IEVAB)=FOMEM(IEVAB)+ESTIF(IEVAB,KOUNT)*TDISP(LOCAL,IDOFN)
                  END DO 
              END DO 
          END DO
!                                                                       
!     EVALUATE THE AXIAL FORCE                                          
!                                                                       
          IF(NDOFN==1) STRES(IELEM,1)= ABS(FOMEM(1))                      
	      IF(NDOFN==2) STRES(IELEM,1)=SQRT(FOMEM(1)*FOMEM(1)+FOMEM(2)*FOMEM(2))
	
!
!     DETERMINE THE CORRECT SIGN                                        
!                                                                       
          IF (FOMEM(1)<0.0) STRES(IELEM,1)=-STRES(IELEM,1)
          
      END DO  !END LOOP OVER EACH ELEMENT
      RETURN                                                            
      END
      
      
!     *******************************************************************
!     SUBROUTINE RESULT:WRITE THE NODAL(JOINT) DISPLACEMENTS AND REACTIONS  
!     *******************************************************************
      SUBROUTINE RESULT_CST
      
      USE FEM4
      USE FEM5
!
!     OUTPUT OF RESULTS 
!
      WRITE(8,900)
900   FORMAT(5X,'NODE',1X,'DISPLACEMENTS',6X,'REACTIONS')
      
      
     Print*,"======================================================================"
     Print*,"                             RESULTS"
     Print*,"======================================================================"
    
                Print*,"         Node           Displacements               Reactions"

      DO  IPOIN=1,NPOIN
          Print*, IPOIN,TDISP(IPOIN,:),TREAC(IPOIN,:)
          WRITE(8,*) IPOIN,TDISP(IPOIN,:),&
              TREAC(IPOIN,:)
      END DO
!                                                                       
!     WRITE THE ELEMENT(MEMBER) STRESSES(FORCES)                        
!             
                Print*,"    Element Internal Force"

      IF(NSTRE/=0) THEN                                           
          WRITE(8,920)                                                      
920       FORMAT(1X,'ELEMENT',6X,'STRESSES')                            
          DO  IELEM=1,NELEM 
              Print*, IELEM,STRES(IELEM,:)
              WRITE(8,*) IELEM,STRES(IELEM,:) 
          END DO
      END IF     
      Pause
      RETURN                                                            
      END	            	 	       	 	 
