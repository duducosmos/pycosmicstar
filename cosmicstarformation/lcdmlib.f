C file:lcdmlib.f
C**********************************************************************************
C   LCDM Cosmological Library:
C
C**********************************************************************************


        subroutine about()

        write(*,*) "A FORTRAN wrapper library for LCDM "
        write(*,*) "cosmology analisys in Python."
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
        write(*,*) "dtdz(z) : Time and Redshift Relation for "
        write(*,*) "  \LambdaCDM Universe"
        write(*,*) "rz(z): comove distancy  for \LambdaCDM Universe "
        write(*,*) "dr_dz(z) : variation of comove distance with z for"
        write(*,*) "   \LambdaCDM Universe"
        write(*,*) "dV_dz(z) : variation of comove volume with redshift"
        write(*,*) "   \LambdaCDM Universe"
        write(*,*) "age(z):  Age of the Universe for LCDM Universe"
        write(*,*) ""

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
C Retorna a variacao do tempo com o redshift
C**********************************************************************************

        real*8 function  dtdz(z)
          real*8 ct3,z
          real*8 omegab,omegam,omegal,h
          common/cparam/ omegab,omegam,omegal,h
Cf2py intent(in) :: z
Cf2py intent(out) :: dtdz

          ct3=9.78d+09/h

          dtdz=ct3/((1.0d+0+z)*dsqrt(omegal+omegam*(1.0d+0+z)**3.0d+0))
          return
        end

C**************************************************************************************
C  Comove Distancy LCDM
C**************************************************************************************
        real*8 function rz(z)
           external G
           real*8 G
           real*8 z,g2,vl
           real*8 omegab,omegam,omegal,h
           common/cparam/ omegab,omegam,omegal,h
cf2py real(DP) a
cf2py a = G(a)
Cf2py intent(in) :: z
Cf2py intent(out) :: rz
           vl  = 3.0d+5 ! speed of light in km/s
           call romberg(G,0.0d+0,z,1.0d-4,g2)
           rz = (vl/(100.0d+0*h))*g2
        end

       real*8 function G(z)
           real*8 z
           real*8 omegab,omegam,omegal,h
           common/cparam/ omegab,omegam,omegal,h
Cf2py intent(in) :: z
Cf2py intent(out) :: G
           G =1.0d+0/dsqrt(omegal+omegam*(1.0d+0+z)**3.0d+0)
        end

C**************************************************************************************
C       Return the comove distance derivative with respect to z
C**************************************************************************************

        real*8 function dr_dz(z)
           real*8 z,dr_dz1,vl
           real*8 omegab,omegam,omegal,h
           common/cparam/ omegab,omegam,omegal,h

cf2py real intent(in) ::   z
cf2py real intent(out) :: dr_dz

           vl  = 3.0d+5 ! speed of light in km/s

           dr_dz1 = (1.0d+0/dsqrt(omegal+omegam*(1.0d+0+z)**3.0d+0))
           dr_dz = (vl/(100.0d+0*h))*dr_dz1
        end

C**************************************************************************************
C   The comove volume variation
C**************************************************************************************

        real*8 function dV_dz(z)
            external rz,dr_dz
            real*8 rz,dr_dz
            real*8 z
            real*8 pi
cf2py real(DP) a
cf2py a =  dr_dz(a)
cf2py real(DP) b
cf2py b  = rz(b)
cf2py real intent(in) :: z
cf2py real intent(out) :: dV_dz

            pi = 3.1415d+0
            dV_dz = 4.0d+0*pi*dr_dz(z)*(rz(z))**2.0d+0

        end

C**************************************************************************************
C Funcao que fornece a idade do universo para um dado redshift
C Para Universo $\Lambda$CDM
C**************************************************************************************
        real*8 function age(z)
          real*8 z
          real*8 a,a3,fct
          real*8 omegalm,hsl
          real*8 omegab,omegam,omegal,h
          common/cparam/ omegab,omegam,omegal,h

c
Cf2py intent(in) :: z
Cf2py intent(out) ::tage
c

          omegalm=omegal/omegam
          hsl=h*sqrt(omegal)
          a=1.0/(1.0+z)
          a3=a**3.0
          fct=omegalm*a3
          age= 6.522916d+09*dlog(dsqrt(fct)+dsqrt(fct+1.0))/hsl
          return
        end



C**********************************************************************************
C  sigma: is the variance of the linear density field. As pointed out by Jenkis et al. (2001),
C              this definition of the mass function has the advantage that it does
C              not explicitly depend on redshift, power spectrum or cosmology.
C
C**********************************************************************************


        subroutine sigma(anorm,alfa1,beta1,gama1,ct2,kmass,km,sg,n)
            external dsigma2_dk
            real*8 dsigma2_dk
            integer i,n
            real*8 anorm,kmass(n),sg(n),km(n)
            real*8 tol
            real*8 sig2_1,sig2_2,sig2_3,sig2_4
            real*8 alfa1,beta1,gama1,ct2
            real*8 escala,alfa,beta,gama
            real*8 t0,t1,t2,t3,t4
            common/dados/ escala,alfa,beta,gama

cf2py real(DP) a
cf2py a = dsigma2_dk(a)
cf2py intent(in) :: anorm,ct2,alfa,beta,gama
cf2py intent(in) :: kmass
cf2py intent(out) :: sg,km
cf2py integer intent(hide), depend(kmass) :: n


            alfa=alfa1
            beta=beta1
            gama=gama1
            tol=1.45d-8

            do  i=1,n
                escala = (kmass(i)/ct2)**(1.0d+0/3.0d+0)
                km(i) = dlog10(kmass(i))

                t0=dlog10(1.0d-7/escala)
                t1=dlog10(1.0d-3/escala)
                t2=dlog10(1.0d+0/escala)
                t3=dlog10(10.0d+0/escala)
                t4=dlog10(100.0d+0/escala)

                call romberg(dsigma2_dk,t0,t1,tol,sig2_1)
                call romberg(dsigma2_dk,t1,t2,tol,sig2_2)
                call romberg(dsigma2_dk,t2,t3,tol,sig2_3)
                call romberg(dsigma2_dk,t3,t4,tol,sig2_4)

                sg(i)=dsqrt(anorm*(sig2_1+sig2_2+sig2_3+sig2_4))

              enddo
        end

C**********************************************************************************
C   dsigma2_dk : Derivative of the variance of the linear density field
C                         with respect to the scala factor
C**********************************************************************************
        real*8 function dsigma2_dk(kl)

c
Cf2py intent(in) :: kl
Cf2py intent(out) ::dsigma2_dk
c
            real*8 kl,k,x,pk1,pk2,pdmk
            real*8 escala,alfa,beta,gama
            common/dados/ escala,alfa,beta,gama
            k = dexp(kl)
            x = escala*k

            pk1 = 1.0d+0+(alfa*k+(beta*k)**1.5d+0
     &            +(gama*k)**2.0d+0)**1.13d+0
            pk2 = 1.0d+0/pk1
            pdmk = pk2*(k**3.0d+0)
            dsigma2_dk=pdmk*(3.0d+0*(dsin(x) - x*dcos(x))
     &                        /(x**3.0d+0))**2.0d+0
            return
        end


C**********************************************************************************
C grow : Growth function. In the more general case of a Universe with matter and
C           a cosmoogical constant, the exact solution for the growth function is given
C           by the work of Carrol, Press & Turner 1992.
C**********************************************************************************

        real*8 function grow(z)
            real*8 z
            real*8 a,a2,a3,a4
            real*8 ea,dz1,pi,s2pi
            real*8 omegamz,omegalz
            real*8 omegab,omegam,omegal,h
            common/cparam/ omegab,omegam,omegal,h


c
Cf2py intent(in) :: z
Cf2py intent(out) :: grow
c

            pi=3.14159265d+0
            s2pi=pi*dsqrt(2.0d+0)
            a=1.0d+0/(1.0d+0+z)
            a2=a*a
            a3=a2*a
            a4=a3*a
            ea = omegam*a+omegal*a4
            omegamz=omegam*a/ea
            omegalz=omegal*a4/ea
            dz1= (1.0d+0-omegalz)+(omegamz**(4.0d+0/7.0d+0))
     &              +(omegamz/2.0d+0)
            grow=(2.5d+0*omegamz*a/dz1)/s2pi

            return
        end

C**********************************************************************************
C rhodm : Evolution of the dark matter density
C**********************************************************************************


        subroutine  rhodm(rhodm0,z,rhodmz)
            real*8 z,rhodm0,rhodmz
            real*8 a,a2,a3

c
Cf2py intent(in) :: z
Cf2py real intent(in) :: rhodm0
Cf2py intent(out) ::rhodmz
c
            a=1.0d+0/(1.0d+0+z)
            a2=a*a
            a3=a2*a
            rhodmz=rhodm0/a3
            return
        end
C**********************************************************************************
C rhobr : Evolution of the barionic matter density.
C**********************************************************************************

        subroutine rhobr(rhob0,z,rhobrz)
         real*8 rhob0,z,rhobrz
         real*8 a,a3
c
Cf2py real intent(in) :: rhob0
Cf2py intent(in) :: z
Cf2py intent(out) :: rhobrz
c
          a=1.0d+0/(1.0d+0+z)
          a3=a*a*a
          rhobrz=rhob0/a3
          return
      end



C**********************************************************************************
C**********************************************************************************
C**********************************************************************************
C                       Numerical Method Section
C**********************************************************************************
C**********************************************************************************
C**********************************************************************************


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
