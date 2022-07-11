program Lab2
      INTEGER i,j,z
      integer, parameter ::NDIM=5, N=5 ,STEPS=3, REALNESS=4
      REAL(REALNESS) :: R(N,N),Original(N,N), NOTA(N,N), A(N,N), U(N,N), L(N,N), X, COND, IPVT(N), WORK(N), DET, B(N), E(N,N)
      do z=1,STEPS
         print*, "STEP",Z
         if (z==1) x=1.01
         if (z==2) x=1.001
         if (z==3) x=1.0001
         print*,"Value of X", X
         forall(i=1:N,j=1:N) A(i,j)=1
         forall (i=1:N,j=1:N,i==j) A(i,j)=x
         Original=A
         print*, "Matrica A"
         print 101, ((A(i,j),j=1,N),i=1,N)
         E=0
         forall(i=1:N,j=1:N,i==j)E(i,j)=1
         call DECOMP(NDIM,N,A,COND,IPVT,WORK)
         print*, "Chislo obyslovlennosti matrici A", COND
         U=0
         forall(i=1:N,j=1:N,j>=i) U(i,j)=A(i,j)
         L=0
         forall(i=1:N,j=1:N,i==j)L(i,j)=1
         forall(i=1:N,j=1:N,J<i) L(i,j)=A(i,j)
         forall(i=1:N,j=1:N,j<i) L(i,J)=-L(i,j)
         DET=IPVT(N)
         do i=1,N
            DET=DET*A(i,i)
         end do
         do i=1,N
            B=0
            B(i)=1
            call SOLVE(NDIM,N,A,B,IPVT)
            NOTA(1:N,i)=B
         end do
         print*, "Obratnaya matrica A^-1"
         print 101,((NOTA(i,j),j=1,N),i=1,N)
         R=matmul(Original,NOTA)
         R=R-E
         print*, "Vichislennaya matrica R"
         print 101,((R(i,j),j=1,N),i=1,N)
         X=norm2(R)
         print*, "Norma matrici R"
         print 102,X
         print *, ''
         call DECOMP(NDIM,N,R,COND,IPVT,WORK)
      end do

 101  FORMAT((5ES14.7))
 102  FORMAT(1ES12.6)
end program Lab2
