PROGRAM SCFT

INTEGER :: i,j,k
REAL*8 :: x,y,z,dens(100,100,100),ascale

OPEN(UNIT=101,file='XE2.dat')
OPEN(UNIT=102,file='out1.xyz')

ascale = 0.80d0 ! This parameter scales the model (i.e. packs the points to different densities)

DO k=1,36
 DO j=1,64
  DO i=1,28

   READ (101,*) dens(i,j,k) 

  END DO
 END DO
END DO

DO k=1,36
 DO j=1,64
  DO i=1,28

    IF(dens(i,j,k).GE.0.75d0) WRITE (102,*) 'H',ascale*DBLE(i),ascale*DBLE(j),ascale*DBLE(k)

  END DO
 END DO
END DO

! Once XYZ file is written, must add two lines to the top of file.
! First line is integer number of points in the file.
! The second line is a comment line.

END PROGRAM SCFT
