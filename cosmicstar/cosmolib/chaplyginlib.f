C file:lcdmlib.f
C**********************************************************************************
C   LCDM Cosmological Library:
C
C**********************************************************************************


        subroutine about()

        write(*,*) "A FORTRAN wrapper library for Chaplygin Gas"
        write(*,*) "cosmology analisys in Python."
        write(*,*) ""
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


        write(*,*) "Initialization: "

        write(*,*) "  import cosmolib"
        write(*,*) "  omegab = 0.04"
        write(*,*) "  omegam = 0.24"
        write(*,*) "  omegal = 0.7"
        write(*,*) "  h = 0.7"
        write(*,*) "myUniverse = cosmolib.init(omegab,omegam,omegal,h)"

        write(*,*) "FUNCTIONS:"
        write(*,*) "dtdzCG: Time and z Relation for Chaplygin Gas"
        write(*,*) "rzGC(z) : comove distancy  for  Chaplygin Gas"
        write(*,*) "drGC_dz(z) : variation of comove distance with z"
        write(*,*) " for  Chaplygin Gas"
        write(*,*) " ageCG(z) : Age of the Universe for  Chaplygin Gas"
        write(*,*) ""
        write(*,*) " "

        end



C**********************************************************************************
C      Inicializa esse modulo com os valores de Omega (b,m e Lambda) e h
C**********************************************************************************

        subroutine init(omegab1,omegam1,omegal1,h1)
C      Inicializa esse modulo com os valores de Omega (b,m e Lambda) e h

            real*8 omegab1,omegam1,omegal1,h1
            real*8 omegab,omegam,omegal,h
            common/cparam/ omegab,omegam,omegal,h
cf2py real intent(in) :: omegab1,omegam1,omegal1,h1
            omegab = omegab1
            omegam = omegam1
            omegal = omegal1
            h = h1
        end


C**********************************************************************************
C Retorna a variacao do tempo com o redshift para o Gas de Chaplygin
C**********************************************************************************

        real*8 function  dtdzCG(z)
          external E
          real*8 E
          real*8 ct3,z
          real*8 omegab,omegam,omegal,h
          common/cparam/ omegab,omegam,omegal,h

cf2py real(DP) a
cf2py a = E(a)
Cf2py intent(in) :: z
Cf2py intent(out) :: dtdz

          ct3=9.78d+09/h

          dtdzCG=ct3/((1.0d+0+z)*E(z))
          return
        end

C**************************************************************************************
C  Comove Distancy Chaplygin Gas
C**************************************************************************************
        real*8 function rzGC(z)
           external E
           real*8 E
           real*8 z,g2,vl
           real*8 omegab,omegam,omegal,h
           common/cparam/ omegab,omegam,omegal,h
cf2py real(DP) a
cf2py a = E(a)
Cf2py intent(in) :: z
Cf2py intent(out) :: rzGC
           vl  = 3.0d+5 ! speed of light in km/s
           call romberg(E,0.0d+0,z,1.0d-4,g2)
           rzGC = (vl/(100.0d+0*h))*g2
        end

C**************************************************************************************
C Expansion Function
C**************************************************************************************
        real*8 function E(z)
             external DEstate
             real*8 DEstate
             real*8 iDES,z,a
             real*8 omegab,omegam,omegal,h
             common/cparam/ omegab,omegam,omegal,h

cf2py real(DP) a
cf2py a = DEstate(a)
Cf2py real intent(in) :: z

             a = 1.0d+0/(1.0d+0+z)
             call romberg(DEstate,1.0d+0,a,1.0d-4,IDES)
             IDES = dexp(-3.0d+0*IDES)

             E = dsqrt(omegam*(1.0d+0)**3.0d+0+omegal*IDES)
             return

        end
C**************************************************************************************
C  derivative Comove Distancy Chaplygin Gas
C**************************************************************************************
        real*8 function drGC_dz(z)
            external E
            real*8 E
            real*8 vl,z
            real*8 omegab,omegam,omegal,h
            common/cparam/ omegab,omegam,omegal,h
            vl  = 3.0d+5 ! speed of light in km/s
            drGC_dz = (vl/(100.0d+0*h))/E(z)

        end

C**************************************************************************************
C (Dark-Energy Equation-of-State +1) /(scala factor) Generalized Chaplygin Gas (alpha = 0.2)
C**************************************************************************************
        real*8 function DEstate(a)
            real*8 a
            real*8 alpha,w0,w,e1

Cf2py real intent(in) ::a
Cf2py real intent(out) :: DEstate

            alpha = 0.2d+0
            w0 = -0.8d+0
            e1 = ((1.0d+0/w0)+1.0d+0)*a**(-3.0d+0*(alpha+1))
            w = -1.0d+0/(1.0d+0-e1)
            DEstate = (1.0d+0+w)/a
            return
        end
C**************************************************************************************
C Funcao que fornece a idade do universo para um dado redshift
C Para Universo com Gas de Chaplygin
C**************************************************************************************
        real*8 function ageCG(z)
          external dtdzCG
          real*8 dtdzCG
          real*8 z,tu
          real*8 omegab,omegam,omegal,h
          common/cparam/ omegab,omegam,omegal,h
cf2py real(DP) a
cf2py a = dtdzCG(a)
Cf2py real intent(in) :: z
Cf2py real intent(out) :: ageCG

          call romberg(dtdzCG,1000.0d+0,z,1.0d-4,tu)
          ageCG = tu
        end




C******************************************************************************
C******************************************************************************
C******************************************************************************
C                  Numerical Method Section
C******************************************************************************
C******************************************************************************
C******************************************************************************


C******************************************************************************
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
