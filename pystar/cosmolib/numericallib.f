C file:lcdmlib.f
C**********************************************************************************
C   LCDM Cosmological Library:
C
C**********************************************************************************


        subroutine about()

        write(*,*) "A FORTRAN wrapper library for numerical analisys in"
        write(*,*) "Python."
        write(*,*) "DEPENDENCE:"
         write(*,*) "python-numpy"

        write(*,*) "INSTALL:"
        write(*,*) "$ python setup.py install"

        write(*,*) "DISCLAIMER"

        write(*,*) "cosmolib is free software: you can redistribute"
         write(*,*) "it and/or modify"
        write(*,*) "it under the terms of the GNU General Public"
        write(*,*) " License as published by"
        write(*,*) "the Free Software Foundation, either "
        write(*,*) "version 3 of the License."
        write(*,*) "cosmolib is distributed in the hope that"
        write(*,*) " it will be useful,"
        write(*,*) "but WITHOUT ANY WARRANTY; without even"
        write(*,*) " the implied warranty of"
        write(*,*) "MERCHANTABILITY or FITNESS FOR A PARTICULAR "
        write(*,*) "PURPOSE.  See the"
        write(*,*) "GNU General Public License for more details."
        write(*,*) " "
        write(*,*) "You should have received a copy of the GNU "
        write(*,*) "General Public License"
        write(*,*) "   along with Foobar.  If not,"
        write(*,*) "   see <http://www.gnu.org/licenses/>"


        write(*,*) "AUTHOR"
        write(*,*) "Eduardo S. Pereira"
        write(*,*) "email: pereira.somoza@gmail.com"

        write(*,*) ""
        write(*,*) " NUMERICAL METHOD FUNCTION:"
        write(*,*) "rk4_in: "
        write(*,*) "4th-order Runge-Kutta method for solving the"
        write(*,*) " initial value problem { y}\' = { F(x,{ y} )}"
        write(*,*) " where { y} = { y[0],y[1],...y[n-1]} ."
        write(*,*) " ARGUMENTS:"
        write(*,*) "y: nitial conditions."
        write(*,*) "Xarray : Array with the x values"
        write(*,*) "Yarray: Output Array with the solutions"
        write(*,*) "fun: user-supplied function that returns"
        write(*,*) "the array F(x,y)"
        write(*,*) " RETURN:"
        write(*,*) "  The integrated numerical function"
        write(*,*) " romberg:"
        write(*,*) " Romberg Integration"
        write(*,*) " ARGUMENTS:"
        write(*,*) "  func : Function to be integrated"
        write(*,*) "  a : start point"
        write(*,*) "  b : end point"
        write(*,*) "  tol : tolerance"
        write(*,*) ""
        write(*,*) " RETURN:"
        write(*,*) "  The integrated value"
        write(*,*) "  "
        write(*,*) "locate:"
        write(*,*) " Localiza a posicao de dado ponto"
        write(*,*) " a partir de dois adjacentes"
        write(*,*) " ARGUMENTS:  "
        write(*,*) " func --- function or entry table"
        write(*,*) " xx   --- entry table"
        write(*,*) " n--- n point in the table"
        write(*,*) " x---  value in x that is related to y"
        write(*,*) " RETURN"
        write(*,*) " j---  x,y position"
        write(*,*) " "
        write(*,*) "dfridr:"
        write(*,*) " Gives the derivate of func with respect to x"
        write(*,*) " Arguments:"
        write(*,*) "func - Function to be derived"
        write(*,*) "x - point in x where the derivative "
        write(*,*) "h - step for diferential"
        write(*,*) "error - internal parameter for function error"
        write(*,*) "RETURN:"
        write(*,*) "The derived value of the funtion in the x point"
        write(*,*) " "

        end



C**********************************************************************************
C  4th-order Runge-Kutta method for solving the
C   initial value problem { y} ’ = { F(x,{ y} )} , where
C   { y} = { y[0],y[1],...y[n-1]} .
C    y    = initial conditions.
C   Xarray => Array with the x value for all range
C   Yarray => Output Array with the solutions of Y'
C    fun      => user-supplied function that returns the
C             array F(x,y) = { y’[0],y’[1],...,y’[n-1]} .
C**********************************************************************************

        subroutine rk4_int(fun,y,Harray,Xarray,Yarray,n)

          external fun
          real*8 fun
          real*8 y(1)
          real*8 x,h
          real*8 k0,k1,k2,k3,r
          integer i,n
          real*8 Xarray(n),Yarray(n),Harray(n)

cf2py real(DP) a
cf2py a = fun(a)
Cf2py real intent(in) :: y
Cf2py real intent(in) :: Xarray
Cf2py real intent(out) :: Yarray
Cf2py iteger intent(in) :: n
Cf2py real inten(in) :: Harray
Cf2py integer intent(hide), depend(Xarray) :: n


          do i=1,n
              Yarray(i)=y(1)
              x=Xarray(i)
              h = Harray(i)


              k0=h*fun(x,y(1))
              k1=h*fun(x+h/2.0d+0,y(1)+k0/2.0d+0)
              k2=h*fun(x+h/2.0d+0,y(1)+k1/2.0d+0)
              k3=h*fun(x+h,y(1)+k2)
              r= (k0+2.0d+0*k1+2.0d+0*k2+k3)/6.0d+0

              y(1)=y(1)+r



          enddo

        end

C**************************************************************************************
C Romberg Integration
C**************************************************************************************
        subroutine romberg(func,a,b,tol,res)
           external func
           real*8 func
           real*8 a,b,tol,res,lnew
           real*8 r(21)
           real*8 r_old
           real*8 tmp, tmp2
           integer k
cf2py real(DP) a
cf2py a = func(a)
cf2py real intent(in) :: func
cf2py real, intent(in) :: a,b,tol
cf2py real, intent(out) :: res
           call trapezoid(func,a,b,0.0d+0,1,lnew)
           r(1) = lnew
           r_old = r(1)
           do k =2,20
               call trapezoid(func,a,b,r(k-1),k,lnew)

               r(k) = lnew

               call richardson(r,k)
               tmp = dabs(r(1)-r_old)
               tmp2 = tol*max(dabs(r(1)),1.0d+0)
               if (tmp .lt. tmp2 )then
                   res = r(1) !,2**(k-1)
                  return
               endif
               r_old = r(1)
           enddo

           print*, "Romberg quadrature did not converge"
       end

       subroutine  richardson(r,k)
cf2py real, intent(in,out) :: r
cf2py integer, intent(in) :: k
           real*8 r(21)
           real*8 const
           integer j,k
           do j=k-1,0,-1
               const = 4.0**(k-j)
               r(j) = (const*r(j+1) - r(j))/(const - 1.0)
           enddo
       end

       subroutine trapezoid(func,a,b,Iold,k,lnew)
          external func
          real*8 func
          real*8 a,b,Iold
          real*8 lnew,h,x,sum
          integer k,n
cf2py real, intent(in) :: a,b,lold
cf2py integer, intent(in) :: k
cf2py real,intent(out) :: lnew
cf2py real(DP) a
cf2py a = func(a)
cf2py real intent(in) :: func
          if( k == 1) then
              lnew = (func(a) + func(b))*(b - a)/2.0d+0
          else
             n = 2**(k -2 )
            ! Number of new points
            h = (b - a)/float(n)
            !Spacing of new points
            x = a + h/2.0d+0
            !Coord. of 1st new point
            sum = 0.0d+0
            do i =0,n-1
                sum = sum + func(x)
                x = x + h
            enddo
            lnew = (Iold + h*sum)/2.0d+0
          endif
       end

C**************************************************************************************
C       Localiza a posicao de dado ponto a partir de dois adjacentes.
C      argumentos:  func --- funcao ou tabela de entrada
C       xx   --- tabela de entrada
C      n    --- numero de pontos da tabela
C       x    --- valor de x que se deseja determinar y
C      j    --- posicao de saida
C**************************************************************************************

         subroutine locate(xx,x,jl,n)
           integer jl,ju,jm
           integer n
           real*8 xx(n)
           real*8 x
Cf2py real intent(in) :: xx,x
Cf2py integer intent(hide), depend(xx) :: n
Cf2py integer intent(out) :: jl
           jl=0
           ju=n+1
             do while(ju-jl .lt. 1)
                 jm=int((ju+jl)/2)
                 if(xx(n).lt. xx(1) .and. x .lt. xx(jm)) then
                        jl=jm
                else
                       ju=jm
                endif
            end do
            jl = jl+1
            return
         end


C**************************************************************************************
C  Gives the derivate of func with respect to x
C   Arguments: func - Function to be derived
C                          x - point in x where the derivative is analysed
C                          h - step for diferential
C                          error - internal parameter for function error
C**************************************************************************************

        REAL*8 FUNCTION dfridr(func,x,h,err)
             EXTERNAL func
              real*8 func
              INTEGER NTAB
              REAL*8 err,h,x,CON,CON2,BIG,SAFE
              PARAMETER (CON=1.4d+0,CON2=CON*CON)
              PARAMETER (BIG=1.d+30,NTAB=10,SAFE=2.0d+0)
              INTEGER i,j
              REAL*8 errt,fac,hh,a(NTAB,NTAB)

cf2py real(DP) a
cf2py a = func(a)
cf2py real intent(in) :: func
cf2py real intent(in) :: x,h,err
cf2py real intent(out) :: dfridr

              if(h.eq.0) then
                  print*,'h must be nonzero in dfridr'
                  dfridr = 0.0d+0
                  return
              endif

              hh=h
              a(1,1)=(func(x+hh)-func(x-hh))/(2.0d+0*hh)
              err=BIG
              do i=2,NTAB
                hh=hh/CON
                a(1,i)=(func(x+hh)-func(x-hh))/(2.0d+0*hh)
                fac=CON2
                do  j=2,i
                  a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.0d+0)
                  fac=CON2*fac
                  errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                  if (errt.le.err) then
                    err=errt
                    dfridr=a(j,i)
                  endif
                enddo
                if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err) then
                    return
                endif
              enddo
              return
        end
