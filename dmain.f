      program dmain
      
      integer nmax
      parameter(nmax=500)
      integer nread, nwrite
      parameter(nread=1,nwrite=2)
c     **********
c
c     Driver for unconstrained problems.
c
c     Subprograms called
c
c       USER ........... dminfg, dminhs
c
c       MINPACK-2 ...... dnmtr, wallclock
c
c       Level 1 BLAS ... dnrm2
c
c     MINPACK-2 Project. March 2000.
c     Argonne National Laboratory.
c     Jorge J. More'.
c
c     **********
      double precision zero, one
      parameter(zero=0.0d0,one=1.0d0)

      character*60 task
      logical scale
      integer n
      integer isave(6)
      double precision delta, f
      double precision x(nmax), g(nmax), fhes(nmax,nmax)
      double precision wa(7*nmax), dsave(7)

c     Tolerances.

      double precision frtol, fatol, fmin

c     Summary information.

      integer  nfev, ngev, nhev

c     Test problems.

      logical search
      character*6 prob
      integer maxfev, nx, ny
      double precision par
      double precision s(nmax), wa1(2*nmax)
      double precision gnorm
      double precision dnrm2

c     Timing.

      double precision wctime1, wctime2, time1, time2
      double precision fgtime, htime, wctime

      external dminfg, dminhs, dnmtr, wallclock, dnrm2

      open (nread,file='nmtr.dat',status='old')
      open (nwrite,file='nmtr.info')

c     Read tolerances.

      read (nread,*) frtol, fatol, fmin
      write (*,*) frtol, fatol, fmin

      do while (1. eq. 1)

c       Read problem data.
      
         read (nread,*) prob, n, nx, ny, par
         write (*,*)    prob, n, nx, ny, par
         
         if (prob(1:4) .eq. 'STOP') then
            close(nread)
            stop
         end if

c        Generate the starting point.

         call dminfg(n,nx,ny,x,f,g,'XS',prob,par,wa1)

c        Initialize variables.

         nfev = 0 
         ngev = 0
         nhev = 0
         fgtime = zero
         htime = zero 

c        Set parameters.

         scale = .false.         
         maxfev = 200
         
c        Start of search.

         write (nwrite,1000) prob, n
         call wallclock(wctime1)
         task = 'START'
         search = .true.
         
         do while (search)

c           Function evaluation.

            if (task .eq. 'F' .or. task .eq. 'START') then
               call wallclock(time1)
               call dminfg(n,nx,ny,x,f,g,'F',prob,par,wa1)
               nfev = nfev + 1
               call wallclock(time2)
               fgtime = fgtime + (time2 - time1)
            end if

c           Evaluate the gradient and the Hessian matrix.

            if (task .eq. 'GH' .or. task .eq. 'START') then
               call wallclock(time1)
               call dminfg(n,nx,ny,x,f,g,'G',prob,par,wa1)
               ngev = ngev + 1
               call wallclock(time2)
               fgtime = fgtime + (time2 - time1)

c              Evaluate the Hessian matrix.

               call wallclock(time1)
               do i = 1, n
                  s(i) = zero
               end do
               do j = 1, n
                  s(j) = one 
                  call dminhs(n,nx,ny,x,s,fhes(1,j),prob,par,wa1)
                  s(j) = zero
               end do
               nhev = nhev + 1
               call wallclock(time2)
               htime = htime + (time2 - time1)

            end if

c           Initialize the trust region bound.

            if (task .eq. 'START') then
               gnorm = dnrm2(n,g,1) 
               delta = zero
               write (nwrite,2000) f, gnorm
            end if

c           Call the optimizer.

            if (search) then
               call dnmtr(n,x,f,g,fhes,nmax,frtol,fatol,fmin,task,delta,
     +              wa(5*n+1),scale,isave,dsave,wa(1),wa(n+1),wa(2*n+1),
     +              wa(3*n+1),wa(4*n+1))
            end if

c           Test for convergence or termination.

            if (nfev .gt. maxfev) then
               search = .false.
               task = 'ERROR: NFEV > MAXFEV'
            end if
            if (task(1:4) .eq. 'CONV') search = .false. 
            if (task(1:4) .eq. 'WARN' .or. task(1:5) .eq. 'ERROR') 
     +      search = .false.

         end do

c        End of search.

         call wallclock(wctime2)
         write (*,*) task

c        Summary information.
         
         call dminfg(n,nx,ny,x,f,g,'G',prob,par,wa1)
         gnorm = dnrm2(n,g,1)
         write (nwrite,3000) nfev, ngev, nhev, f, gnorm
         
c        Timing information.
         
         wctime = wctime2 - wctime1
         fgtime = 100*(fgtime/wctime)
         htime  = 100*(htime/wctime)
         write (nwrite,4000) wctime, fgtime, htime, task

      end do

 1000 format (' Problem ', a6,                                   //,
     +        ' Number of variables                         ',i12)

 2000 format (
     +        ' Function value at initial iterate        '   ,d15.8,/,
     +        ' Gradient norm at initial iterate         '   ,d15.3)

 3000 format (
     +        ' Number of function evaluations              ',i12,/,
     +        ' Number of gradient evaluations              ',i12,/,
     +        ' Number of Hessian evaluations               ',i12,/,
     +        ' Function value at final iterate          '   ,d15.8,/,
     +        ' Gradient norm at final iterate           '   ,d15.3,/)

 4000 format (' Total execution time                        ',f12.2,/,
     +        ' Percentage in function evaluations          ',f12.2,/,
     +        ' Percentage in Hessian evaluations           ',f12.2,//,
     +        ' Exit message     '                           ,a60,/)

      end

