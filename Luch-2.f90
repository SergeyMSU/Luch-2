!  Luch2.f90 
!
!  FUNCTIONS:
!  Luch2 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Luch2
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

module PARAMETRS
	USE OMP_LIB
	integer(4) :: GD_par_n_points = 500 !400
	real(8) :: par_L = -1.0_8 !-5.0_8
	real(8) :: par_R = 0.0_8
	real(8) :: Mach_inf = 1.1_8  !2.0

	real(8) :: ddx = 0.0_8

	real(8), parameter :: par_Velosity_inf = -2.54351! -1.0 !-2.58_8
	real(8), parameter :: par_a_2 = 0.13043_8 ! 0.13043_8 
	real(8) :: par_n_H = 1.0_8 

	real(8), parameter :: par_pi_8 = acos(-1.0_8)         
	real(8), parameter :: par_pi = acos(-1.0_8)         
	real(8), parameter :: pi = acos(-1.0_8)            
	real(8), parameter :: par_sqrtpi = sqrt(par_pi_8)
	real(8), parameter :: sqrtpi = sqrt(par_pi_8)
    real(8), parameter :: cpi4 = 12.56637061435917295384
    real(8), parameter :: ggg = (5.0/3.0)
    real(8), parameter :: GD_ggg = (5.0/3.0)


	real(8), allocatable :: GD_mas_X(:)     ! Центры ячеек
	real(8), allocatable :: GD_mas_Q2(:)
	real(8), allocatable :: GD_mas_Q3(:)
	real(8), allocatable :: GD_mas_Q2_do(:)
	real(8), allocatable :: GD_mas_Q3_do(:)
	real(8), allocatable :: GD_mas_V(:)
	real(8), allocatable :: GD_mas_P(:)

	real(8), allocatable :: mas_n_H(:)
	real(8), allocatable :: mas_v_H(:)
	real(8), allocatable :: mas_vy_H(:)
	real(8), allocatable :: mas_vz_H(:)
	real(8), allocatable :: mas_T_H(:)

	integer (kind=omp_lock_kind), allocatable :: mas_lock(:)

end module PARAMETRS
	
module Monte_Karlo  
	USE PARAMETRS
	implicit none

    integer(4), parameter :: par_n_potok = 32!24
    integer(4), parameter :: par_n_claster = 1
	integer(4), allocatable :: sensor(:, :, :)
	
	real(8), parameter :: x0 = -0.5_8
	
	real(8), allocatable :: dist_f(:, :)
	real(8) :: vL = -6.0! -4.5                  ! Лучше симметричные концы, чтобы ноль был границей ячейки (или переделывать)
	real(8) :: vR = 4.0! 3.5
	integer(4) :: vN = 224! 224!224   !! 1 если не хотим считать

	real(8), allocatable :: dist_2D_f1(:, :)
	real(8) :: wL = 0.0                  ! Лучше симметричные концы, чтобы ноль был границей ячейки (или переделывать)
	real(8) :: wR = 2.0
	integer(4) :: cell_f1 = 432

	real(8), allocatable :: dist_2D_f2(:, :)
	integer(4) :: cell_f2 = 426

	real(8), allocatable :: dist_2D_f3(:, :)
	integer(4) :: cell_f3 = 410

	integer(4), parameter :: shock_n = 30
	integer(4) :: shock_cell
	real(8) :: shock_v(shock_n)
	real(8) :: shock_Q2(shock_n)
	real(8) :: shock_Q3(shock_n)

	contains
	
	subroutine Dist_funk(Vx, mu, cell)
		real(8), intent(in) :: Vx, mu
		real(8) :: dd
		integer :: n1, i, cell
		
		dd = (vR - vL)/vN
		i = INT(100.0 + (Vx - vL)/dd) - 99

		! print*, (vL + (i - 0.5) * (vR - vL)/vN), Vx
		! print*, (vL + (i - 1 - 0.5) * (vR - vL)/vN), (vL + (i - 0.5) * (vR - vL)/vN), (vL + (i + 1 - 0.5) * (vR - vL)/vN)
		! print*, "___"
		! pause

		if(i < 1) return !i = 1
		if(i > vN) return !i = vN
		dist_f(i, cell) = dist_f(i, cell) + (mu);
		!dist_f(i) = dist_f(i) + (mu/dabs(Vx));
	
	end subroutine Dist_funk

	subroutine Dist_funk2D(Vx, Vy, mu, nn)
		real(8), intent(in) :: Vx, mu, Vy
		real(8) :: dd
		integer :: n1, i, cell, j, nn
		
		dd = (vR - vL)/vN
		i = INT(100.0 + (Vx - vL)/dd) - 100

		dd = (wR - wL)/vN
		j = INT(100.0 + (Vy - wL)/dd) - 100

		if(i < 1) return !i = 1
		if(i > vN) return !i = vN
		if(j < 1) return !i = 1
		if(j > vN) return !i = vN

		if(nn == 1) dist_2D_f1(i, j) = dist_2D_f1(i, j) + (mu)
		if(nn == 2) dist_2D_f2(i, j) = dist_2D_f2(i, j) + (mu)
		if(nn == 3) dist_2D_f3(i, j) = dist_2D_f3(i, j) + (mu)
		
	
	end subroutine Dist_funk2D
	
	subroutine Get_sensor_sdvig(sdvig)
		! пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ пїЅпїЅ пїЅпїЅпїЅпїЅпїЅ
		! Variables
		integer, intent(in) :: sdvig
		logical :: exists
		integer(4) :: i, a, b, c, j
		
		
		inquire(file="rnd_my.txt", exist=exists)
		
		if (exists .eqv. .False.) then
			pause "net faila!!!  345434wertew21313edftr3e"
			STOP "net faila!!!"
		end if
		
		if (par_n_claster * par_n_potok * 2 > 1021) then
			print*, "NE XVATAET DATCHIKOV 31 miuhi8789pok9"
			pause
		end if
		
		open(1, file = "rnd_my.txt", status = 'old')
		
		do i = 1, sdvig
			read(1,*) a, b, c
			read(1,*) a, b, c
		end do
		
		do i = 1, par_n_potok
			read(1,*) a, b, c
			sensor(:, 1, i) = (/ a, b, c /)
			read(1,*) a, b, c
			sensor(:, 2, i) = (/ a, b, c /)
		end do
		
		close(1)
	
	end subroutine Get_sensor_sdvig
	
	subroutine M_K_Set()
		integer(4) :: i
		allocate(sensor(3, 2, par_n_potok))
		!allocate(dist_f(vN, int(GD_par_n_points/1) + 1))
		allocate(dist_f(vN, GD_par_n_points))
		allocate(dist_2D_f1(vN, vN))
		allocate(dist_2D_f2(vN, vN))
		allocate(dist_2D_f3(vN, vN))
		dist_f = 0.0_8
		sensor = 0.0_8
		dist_2D_f1 = 0.0_8
		dist_2D_f2 = 0.0_8
		dist_2D_f3 = 0.0_8

		ALLOCATE(mas_n_H(GD_par_n_points))
		ALLOCATE(mas_v_H(GD_par_n_points))
		ALLOCATE(mas_vy_H(GD_par_n_points))
		ALLOCATE(mas_vz_H(GD_par_n_points))
		ALLOCATE(mas_T_H(GD_par_n_points))

		mas_n_H = 0.0
		mas_v_H = 0.0
		mas_vy_H = 0.0
		mas_vz_H = 0.0
		mas_T_H = 0.0

		do i = 1, shock_n
			shock_v(i) = 0.55 + (i - 1) * (0.75 - 0.55)/(shock_n - 1)
		end do

		shock_Q2 = 0.0
		shock_Q3 = 0.0
		
	end subroutine M_K_Set
	
	subroutine M_K_rand(s1, s2, s3, b)
		integer(4), intent(in out) :: s1, s2, s3
		real(8), intent(out) :: b
		integer(4):: ic15 = 32768, ic10 = 1024
		integer(4):: mz = 710, my = 17784, mx = 11973
		real(8):: xi = 9.0949470177292824E-13, c = 1.073741824E9
		integer(4) :: i13, i12, i11, ii
		i13 = mz * s1 + my * s2 + mx * s3
		i12 = my * s1 + mx * s2
		i11 = mx * s1
		ii = i11 / ic15
		i12 = i12 + ii
		s1 = i11 - ic15 * ii
		ii = i12 / ic15
		i13 = i13 + ii
		s2 = i12 - ic15 * ii
		s3 = mod(i13,ic10)
		b = xi * (c * s3 + ic15 * s2 + s1)
	end subroutine M_K_rand
	
    subroutine MK_Velosity_initial2(potok, Vx, Vy, Vz, UU)
		
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		real(8), intent(in) :: UU
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(1.0 - ksi2))
		Vy = a * cos(2.0 * par_pi_8 * ksi1)
		Vz = a * sin(2.0 * par_pi_8 * ksi1)
		
		z = 0
		p1 = dabs(UU) * par_sqrtpi / (1.0 + dabs(UU) * par_sqrtpi);

		do
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi5)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi_8 * ksi5) * sqrt(-log(ksi4))
			else
				if(ksi4 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi4))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi4)))
				end if
			end if
			
			if((dabs(z + UU) / (dabs(UU) + dabs(z)) > ksi6 .and. z < -UU)) EXIT
		end do

		Vx = z + UU
		
		if (Vx >= 0.0) then
			print*, "Error iuygvbnmklo9890pljiouytrtyjhg"
		end if
		
		return
	
	end subroutine MK_Velosity_initial2

	subroutine MK_Velosity_initial3(potok, Vx, Vy, Vz, Vinf)
		! Запуск с левой границы
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		real(8), intent(in) :: Vinf  ! Скорость на левой границе (обезразмеренна на местное cp)
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
		call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(1.0 - ksi2))
		Vy = a * cos(2.0 * par_pi_8 * ksi1)
		Vz = a * sin(2.0 * par_pi_8 * ksi1)
		
		z = 0
		p1 = dabs(Vinf) * par_sqrtpi / (1.0 + dabs(Vinf) * par_sqrtpi);

		do
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi3)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi5)
			call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi_8 * ksi5) * sqrt(-log(ksi4))
			else
				if(ksi4 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi4))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi4)))
				end if
			end if
			
			if((dabs(z + Vinf) / (dabs(Vinf) + dabs(z)) > ksi6 .and. z > -Vinf)) EXIT
		end do

		Vx = z + Vinf
		
		if (Vx <= 0.0) then
			print*, "Error y6576m78onw4w4y45v4w45"
		end if
		
		return
	
	end subroutine MK_Velosity_initial3
	
	subroutine M_K_Change_Velosity4(potok, Ur, Uthe, Uphi, Vr, Vthe, Vphi, Wr_, Wthe_, Wphi_, cp)
		! пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ, пїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅ пїЅпїЅ-пїЅпїЅпїЅпїЅпїЅпїЅ
		! пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ I_ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ
		
		integer(4), intent(in) :: potok
		real(8), intent(in) ::  Ur, Uthe, Uphi, Vr, Vthe, Vphi, cp
		real(8), intent(out) ::   Wr_, Wthe_, Wphi_
		
		real(8) :: X
		real(8) :: ksi, gam1, gam2, Wr1, Wr2, ksi1, ksi2, W1, W2, Wa, pp1, pp2, c, p, u
		real(8) :: p4, om1, om2, om3, lo, y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h, ksi3, ksi4, ksi5, ksi6, D, ko, gg
		integer(4) :: ii, k
		
		X = sqrt( (Vr - Ur)**2 + (Vthe - Uthe)**2 + (Vphi - Uphi)**2 )

		! пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ
		p4 = 0.5 * par_sqrtpi * X / (1.0 + 0.5 * par_sqrtpi * X)
	
		do while(.True.)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi3)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi4)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi5)
			call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi6)
			
			if (p4 < ksi1) then
				om1 = 1.0 - 2.0 * ksi4
				om2 = sqrt(1.0 - (om1**2)) * cos(2.0 * par_pi_8 * ksi5)
				om3 = sqrt(1.0 - (om1**2)) * sin(2.0 * par_pi_8 * ksi5)

				lo = sqrt(-log(ksi2 * ksi3))
				y1 = lo * om1
				y2 = lo * om2
				y3 = lo * om3
			else
				y1 = sqrt(-log(ksi2)) * cos(par_pi_8 * ksi3)
				y2 = sqrt(-log(ksi4)) * cos(2.0 * par_pi_8 * ksi5)
				y3 = sqrt(-log(ksi4)) * sin(2.0 * par_pi_8 * ksi5)
			end if
			
			v1 = y1 + Ur
			v2 = y2 + Uthe
			v3 = y3 + Uphi
			u1 = Vr - v1
			u2 = Vthe - v2
			u3 = Vphi - v3
			uuu = sqrt(kvv(u1, u2, u3))
			yy = sqrt(kvv(y1, y2, y3))
			h = ((uuu * MK_sigma2(uuu, cp)) / (MK_sigma2(X, cp) * (X + yy)))
			
			if (h >= ksi6) EXIT
		end do


		Wr_ = v1
		Wthe_ = v2
		Wphi_ = v3

		return

	end subroutine M_K_Change_Velosity4
	
	real(8) pure function kvv(x, y, z)
	    real(8), intent (in) :: x, y, z
	    kvv = x**2 + y**2 + z**2
    end function kvv
	
	function MK_sigma(x)
		implicit none
		real(8), intent(in) :: x
		real(8) :: MK_sigma
	
		!MK_sigma = (1.0 - par_a_2 * log(x))**2
		!MK_sigma = (1.0 - 0.135838_8 * log(x/1.33301_8))**2
		MK_sigma = 1.0_8
	end function MK_sigma
	
	real(8) pure function MK_sigma2(x, y)
		real(8), intent (in) :: x, y
		
		!MK_sigma2 = (1.0 - par_a_2 * log(x * y))**2
		MK_sigma2 = (1.0 - 0.135838_8 * log(x * y/1.33301_8))**2
		!MK_sigma2 = 1.0_8
	end function MK_sigma2
	
	real(8) pure function MK_Velosity_1(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_1 = 2.0 * cp / par_sqrtpi + 2.0 * u * u / (3.0 * cp * par_sqrtpi) - &
				u * u * u * u / (15.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_1 =  exp(-u * u / cp**2) * cp / par_sqrtpi + (u + (cp**2) / (2.0 * u)) * erf(u / cp)
		end if
		
	end function MK_Velosity_1

	real(8) pure function MK_Velosity_2(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_2 = (8.0 / 3.0) * cp**4 * par_pi * u + (8.0 / 15.0) * cp**2 * par_pi * u * u * u - &
				(4.0 / 105.0) * par_pi * u**5
		else
			MK_Velosity_2 =  cp**3 * par_pi * (exp(-u * u / cp**2) * cp * u * 2.0 * (cp**2 + 2.0 * u**2) + &
			par_sqrtpi * (4.0 * u**4 + 4.0 * cp * cp * u**2 - cp**4) * erf(u / cp)) / (4.0 * u * u)
		end if
		
	end function MK_Velosity_2
	
	real(8) pure function MK_Velosity_3(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_3 = 8.0 * cp / (3.0 * par_sqrtpi) + 8.0 * u**2 / (9.0 * cp * par_sqrtpi) - &
				44.0 * u**4 / (135.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_3 =  exp(-u**2 / cp**2) * cp * (5.0 * cp**2 + 2.0 * u**2) / (par_sqrtpi * (3.0 * cp**2 + 2.0 * u**2)) + &
			(4.0 * u**4 + 12.0 * cp**2 * u**2 + 3.0 * cp**4) * erf(u / cp) / (2.0 * u * (3.0 * cp**2 + 2.0 * u**2))
		end if
		
	end function MK_Velosity_3
	
	real(8) pure function MK_int_1_f1(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_1_f1 =  6.283185155644284 + 0.000024846677279866114 * x + &
				2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x - &
				0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x + &
				0.006992274370591744 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f1 =  6.437524091973454 - 0.6331520099380095 * x + &
				3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x + &
				0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x - &
				0.000920750276197054 * x * x * x * x * x * x
		else if (x <= 5) then
			MK_int_1_f1 =  4.4920780630505135 + 2.5133093267020654 * x + &
				1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x + &
				0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x + &
				0.00005954097505746697 * x * x * x * x * x * x
		else if (x <= 7) then
			MK_int_1_f1 =  1.9138683588136232 + 5.350374732905213 * x - &
				0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x - &
				0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x - &
				0.00002691859374483661 * x * x * x * x * x * x
		else if (x <= 50.0) then
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		else
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		end if
				 
		return
	end function MK_int_1_f1
	
	real(8) pure function MK_int_1_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_1_f2 =  1.328216167939543 - 0.000004545681954848391 * x + &
				2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x - &
				0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x + &
				0.002767736179209713 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f2 =  1.2959616295629246 + 0.1533684067037866 * x + &
				2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x - &
				0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x - &
				0.0029978773724628604 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f2 =  1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x - &
				0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x - &
				0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f2 =  -4.484415105552316 + 5.3747429756690766 * x + &
				0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x - &
				0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x - &
				0.00006928845984076111 * x * x * x * x * x * x
		end if
				 
		return
	end function MK_int_1_f2
	
	real(8) pure function MK_int_1_f3(x)
		real(8), intent (in) :: x
 
		if (x <= 1.0) then
			MK_int_1_f3 = 1.2938345594193854 - 0.000031719847351174835 * x + &
				1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x + &
				0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x - &
				0.015347864048406958 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f3 = 0.9667460440956788 + 1.336271810704016 * x - &
				0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x - &
				0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x + &
				0.001780776080720323 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f3 = 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x + &
				0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x + &
				0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f3 = 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x - &
				0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x - &
				0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x
		end if
			 
		return
	end function MK_int_1_f3
	
	real(8) pure function MK_int_1(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_1 = (cp / (par_sqrtpi**3)) * (b * b * MK_int_1_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_1_f2(x / cp) + par_a_2**2 * MK_int_1_f3(x / cp))
	end function MK_int_1

	real(8) pure function MK_int_2_f1(x)
		real(8), intent (in) :: x
	
		if (x <= 1.0) then
				MK_int_2_f1 =   8.377571213788123 * x + 0.00047608508679086725 * x * x + &
					1.6710478320575737 * x * x * x + 0.016857530811432916 * x * x * x * x - &
					0.15132474125724812 * x * x * x * x * x + 0.030723378194358945 * x * x * x * x * x * x
				return
		else if (x <= 3.0) then
				MK_int_2_f1 =   -0.11788367995598747 + 8.937936705157014 * x - &
						1.119886471634001 * x * x + 2.8831031948885917 * x * x * x - &
						0.735250146386514 * x * x * x * x + 0.10356311378423572 * x * x * x * x * x - &
						0.006231417172309398 * x * x * x * x * x * x
				return
		else if (x <= 5.0) then
				MK_int_2_f1 =   2.9044497739429858 + 2.757712415967557 * x + 4.239161941189675 * x * x + &
						0.36198838294786784 * x * x * x - 0.05737777787138304 * x * x * x * x + &
						0.004956250079677106 * x * x * x * x * x - 0.0001809238236975877 * x * x * x * x * x * x
				return
		else if (x <= 7.0) then
				MK_int_2_f1 =  41.6323689028892 - 38.118317864344135 * x + 22.211189528076645 * x * x -&
					3.8547348524931246 * x * x * x + 0.5000517174807501 * x * x * x * x -&
					0.03446294709493891 * x * x * x * x * x + 0.0009860204070962582 * x * x * x * x * x * x
				return
		end if
	
		MK_int_2_f1 =  0.0
	end  function MK_int_2_f1
	
	real(8) pure function MK_int_2_f2(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_2_f2 =  3.8653461103376667 * x + 0.0001975300512691014 * x * x + &
				2.4468141895384012 * x * x * x + 0.005987984681429616 * x * x * x * x - &
				0.06453987836713967 * x * x * x * x * x + 0.0066920981111229004 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f2 =  -0.10983889480661446 + 4.321087890898017 * x - &
				0.7707850845797619 * x * x + 3.1237901158486583 * x * x * x - &
				0.31485222316123385 * x * x * x * x + 0.010270249760261791 * x * x * x * x * x + &
				0.0008259803934338584 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f2 =  1.8468011847509729 + 1.1986396743254275 * x + &
				1.1421489029509448 * x * x + 2.606316149569781 * x * x * x - &
				0.2788783468089509 * x * x * x * x + 0.019815317035281846 * x * x * x * x * x - &
				0.0006379970557448899 * x * x * x * x * x * x
			return
		else if (x <= 7.0) then
			MK_int_2_f2 =  9.480707804348185 - 8.022988228952784 * x + 5.823555900242488 * x * x + &
				1.3277220473440972 * x * x * x - 0.08074921612981413 * x * x * x * x + &
				0.0033058587723954185 * x * x * x * x * x - 0.00006041810279926061 * x * x * x * x * x * x
			return
		end if
			 
		MK_int_2_f2 = 0.0
	end function MK_int_2_f2
		
	real(8) pure function MK_int_2_f3(x)
        real(8), intent (in) :: x
		if (x <= 1.0) then
			MK_int_2_f3 =  2.6106039258326 * x - 0.0008357997793049243 * x * x + &
				2.0764571907368174 * x * x * x - 0.03182275644841273 * x * x * x * x + &
				0.26310521962808975 * x * x * x * x * x - 0.06034325471643871 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f3 =   0.20784760901369737 + 1.5920325291316857 * x + &
				2.0985329535259014 * x * x - 0.26286255221171206 * x * x * x + &
				1.4610329096120886 * x * x * x * x - 0.25626862029131897 * x * x * x * x * x + &
				0.01684969647300594 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f3 =   -6.284115352064703 + 15.665162343948523 * x &
				- 10.766105772158252 * x**2 + 6.0821074614870465 * x**3 - &
				0.3181501403196319 * x**4 + 0.01232319194701587 * x**5 &
				- 0.00017890661550876597 * x**6
			return
		else if (x <= 7.0) then
			MK_int_2_f3 =   -4.355962170454177 + 13.835332665069274 * x &
				- 10.12766071646978 * x**2 + 5.999392227482686 * x**3 - &
				0.32171318647989955 * x**4 + 0.014181987261856027 * x**5 &
				- 0.00030579035551447497 * x**6
			return
		end if
				
		MK_int_2_f3 = 0.0
		
	end function MK_int_2_f3
	
	real(8) pure function MK_int_2(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_2 = -(cp * cp / (par_sqrtpi * par_sqrtpi * par_sqrtpi)) * (b * b * MK_int_2_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_2_f2(x / cp) + (par_a_2**2) * MK_int_2_f3(x / cp))
	end function MK_int_2
	
	real(8) pure function MK_int_3(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_3 = (cp**3 / (par_sqrtpi**3)) * (b * b * MK_int_3_f1(x / cp) - 2.0 * par_a_2 * b * MK_int_3_f2(x / cp) + &
					 (par_a_2**2) * MK_int_3_f3(x / cp))
	end function MK_int_3

	real(8) pure function MK_int_3_f1(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f1 = 12.566370586001975 - 0.00001944816384202852 * x + 12.567558607381049 * x**2 - &
				0.010507444068349692 * x**3 + 1.2911398125420694 * x**4 - 0.05048482363937502 * x**5 - &
				0.029947937607835744 * x**6
		else if (x <= 3.0) then
			MK_int_3_f1 = 12.451555448799724 + 0.40252674353016715 * x + 12.081033298182223 * x**2 + 0.12939193776331415 * x**3 + &
				1.478876561367302 * x**4 - 0.2237491583356496 * x**5 + 0.014474521138620033 * x**6
		else if (x <= 5.0) then
			MK_int_3_f1 = 8.281962852844913 + 9.783032429527339 * x + 3.165848344887614 * x**2 + 4.711462666119614 * x**3 + &
				0.13739453130827392 * x**4 - 0.012096015315889195 * x**5 + 0.0004514225555943018 * x**6
		else if (x <= 7.0) then
			MK_int_3_f1 = 17.29966669035025 + 2.1457227895928668 * x + 5.572082327920818 * x**2 + 4.4083748449004645 * x**3 + &
				0.13640200155890422 * x**4 - 0.00854302147917508 * x**5 + 0.00022205921430504255 * x**6
		end if
		
	end function MK_int_3_f1

	real(8) pure function MK_int_3_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f2 = 5.798024979296493 - 0.00001772671478406096 * x + 9.98769025405073 * x**2 - 0.0073593516014156535 * x**3 + &
				2.27820901023418 * x**4 - 0.03135086958655956 * x**5 - 0.030403716978821237 * x**6
		else if (x <= 3.0) then
			MK_int_3_f2 = 5.864728705834779 - 0.34799875480550213 * x + 10.760249318358127 * x**2 - 0.9493738978205943 * x**3 + &
				2.948810798239708 * x**4 - 0.29690823082625284 * x**5 + 0.015284639719532296 * x**6
		else if (x <= 5.0) then
			MK_int_3_f2 = -3.21810405152587 + 17.08054504204813 * x - 3.43427364926184 * x**2 + 5.342946243936558 * x**3 + &
				1.3459970589832742 * x**4 - 0.07448103623442687 * x**5 + 0.0021616888853670958 * x**6
		else if (x <= 7.0) then
			MK_int_3_f2 = -59.42860988375287 + 79.24709914283142 * x - 32.304295129821256 * x**2 + 12.558114389999929 * x**3 + &
				0.32138208180756495 * x**4 + 0.003976327853989966 * x**5 - 0.0003703205838879336 * x**6
		end if
		
	end function MK_int_3_f2

	real(8) pure function MK_int_3_f3(x)
		real(8), intent(in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f3 = 3.915885322797866 + 0.000044284982651632276 * x + 7.778946280540595 * x**2 + 0.01990833982643192 * x**3 + &
				2.710531013102172 * x**4 + 0.09480498076917274 * x**5 + 0.026454483470905545 * x**6
		else if (x <= 3.0) then
			MK_int_3_f3 = 4.371710139648201 - 1.787430730965795 * x + 10.530386754306065 * x**2 - 1.9945949141813966 * x**3 + &
				3.310710912254515 * x**4 + 0.1356460090430464 * x**5 - 0.019853464614841342 * x**6
		else if (x <= 5.0) then
			MK_int_3_f3 = 4.6940208896644435 - 6.423081982183987 * x + 18.137713177954524 * x**2 - 7.298291873534797 * x**3 + &
				5.202095367169497 * x**4 - 0.20620743501245353 * x**5 + 0.005094092519200813 * x**6
		else if (x <= 7.0) then
			MK_int_3_f3 = 194.2826492092263 - 195.77327124311523 * x + 95.83253215419927 * x**2 - 23.99281397751645 * x**3 + &
				7.170549820159753 * x**4 - 0.32568903981020303 * x**5 + 0.007955090180048264 * x**6
		end if
		
	end function MK_int_3_f3
	
end module Monte_Karlo

module GD
	USE PARAMETRS
	USE Monte_Karlo
	implicit none
	
	contains
	
	subroutine Set_GD()
		real(8) :: a1, a2
		integer(4) :: i
		
		allocate(GD_mas_X(GD_par_n_points))
		allocate(GD_mas_Q2(GD_par_n_points))
		allocate(GD_mas_Q3(GD_par_n_points))
		allocate(GD_mas_Q2_do(GD_par_n_points))
		allocate(GD_mas_Q3_do(GD_par_n_points))
		allocate(GD_mas_V(GD_par_n_points))
		allocate(GD_mas_P(GD_par_n_points))
		
		allocate(mas_lock(GD_par_n_points))

		do i = 1, GD_par_n_points
			call omp_init_lock(mas_lock(i))
		end do

		ddx = (par_R - par_L)/(GD_par_n_points)

		do i = 1, GD_par_n_points
			GD_mas_X(i) = par_L + (i - 1) * (par_R - par_L)/(GD_par_n_points) + ddx/2.0
			if(GD_mas_X(i) < 0.0) shock_cell = i
		end do


		a1 = (ggg + 1) * Mach_inf**2 / (  (ggg - 1) * Mach_inf**2 + 2.0 )
		a2 = 2 * ggg * Mach_inf**2 / (ggg + 1) - (ggg - 1)/(ggg + 1)

		do i = 1, GD_par_n_points
			! if (GD_mas_X(i) < 0.0) then !(.False.) then!
			! 	GD_mas_P(i) = 1.0/(Mach_inf**2 * ggg) * a2
			! 	GD_mas_V(i) = par_Velosity_inf/ a1
			! else
			! 	GD_mas_V(i) = par_Velosity_inf
			! 	GD_mas_P(i) = 1.0/(Mach_inf**2 * ggg)
			! end if
			GD_mas_P(i) = 0.3
			GD_mas_V(i) = -0.33333333333
		end do

		GD_mas_Q2 = 0.0
		GD_mas_Q3 = 0.0
		GD_mas_Q2_do = 0.0
		GD_mas_Q3_do = 0.0
	end subroutine Set_GD

	subroutine Print_GD(step_)
		integer(4) :: i, step, j
		integer(4), intent(in), optional :: step_
		character(len=3) :: name, name2
		real(8) :: M1, ro1, dL, dR, p1, V1, sred, ro, ro_start, w
		real(8) :: rho_eff, p, p_eff, M, roro, roro1, pp1, pp, u, Ux, cp, uu, cp2, uu2, Ux2

		step = 0
		if(PRESENT(step_)) step = step_

		write(unit=name,fmt='(i3.3)') step

		rho_eff = 1.0 + par_n_H
		p = 1.0/(ggg * Mach_inf**2)
		p_eff = p * (1.0 + par_n_H/2.0)
		M = 1.0/sqrt(ggg * p_eff/rho_eff)

		roro1 = 1.0 * (ggg + 1) * M**2 / ((ggg - 1) * M**2 + 2)
		pp1 = p * ( 2 * ggg * M**2/(ggg + 1) - (ggg - 1)/(ggg + 1) )

		open(1, file = 'Print_GD_' // name //'.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = X, rho, V, p, T, Q2, Q3, nH, VH, TH, ro_eff, pp_eff"

		do i = 1, GD_par_n_points
			if(GD_mas_X(i) >= 0.0) then
				roro = 1.0
				pp = p
			else
				roro = roro1
				pp = pp1
			end if

			write(1,*) GD_mas_X(i), -1.0/GD_mas_V(i), GD_mas_V(i), GD_mas_P(i), -GD_mas_P(i) * GD_mas_V(i)/2.0, GD_mas_Q2(i), GD_mas_Q3(i), & 
				mas_n_H(i), mas_v_H(i), mas_T_H(i), roro, pp
		end do

		print*, "par_n_H = ", par_n_H, " delta = ", -(roro1 + 1.0/GD_mas_V(1)), roro1
		print*, "Max, Max_eff = ", Mach_inf, M

		close(1)

		Ux = par_Velosity_inf
		cp = sqrt(-GD_mas_P(GD_par_n_points) * GD_mas_V(GD_par_n_points))
		Ux2 = GD_mas_V(1)
		cp2 = sqrt(-GD_mas_P(1) * GD_mas_V(1))

		open(2, file = 'info_points.txt')
			do j = 1, size(dist_f(1, :))
				write(2,*) j, GD_mas_X(j)
			end do
		close(2)
		
		if(.False.) then
			do j = 1, size(dist_f(1, :))
				write(unit=name2,fmt='(i3.3)') j
				open(1, file = 'dist_f_' // name2 //'.txt')
				write(1,*) "TITLE = 'HP'  VARIABLES = u, f, f1, f2"
				
				! print*, "par_L = ", par_L
				! print*, "par_R = ", par_R

				do i = 1, vN
					u = (vL + (i + 0.5) * (vR - vL)/vN)
					uu = (sqrt(par_pi) * cp)**(-1.0) * exp(-(u - Ux)**2/cp**2)
					uu2 = (sqrt(par_pi) * cp2)**(-1.0) * exp(-(u - Ux2)**2/cp2**2)
					WRITE (1, *) u, dist_f(i, j), uu, uu2 * (-1.0/GD_mas_V(1))
				end do
				close(1)
			end do
		end if

		open(1, file = 'dist_f_2D_1.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = Vx, Vy, ff"
		do j = 1, vN
			do i = 1, vN
				u = (vL + (i + 0.5) * (vR - vL)/vN)
				w = (wL + (j + 0.5) * (wR - wL)/vN)
				WRITE (1, *) u, w, dist_2D_f1(i, j)/w
			end do
		end do
		close(1)

		open(1, file = 'dist_f_2D_2.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = Vx, Vy, ff"
		do j = 1, vN
			do i = 1, vN
				u = (vL + (i + 0.5) * (vR - vL)/vN)
				w = (wL + (j + 0.5) * (wR - wL)/vN)
				WRITE (1, *) u, w, dist_2D_f2(i, j)/w
			end do
		end do
		close(1)

		open(1, file = 'dist_f_2D_3.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = Vx, Vy, ff"
		do j = 1, vN
			do i = 1, vN
				u = (vL + (i + 0.5) * (vR - vL)/vN)
				w = (wL + (j + 0.5) * (wR - wL)/vN)
				WRITE (1, *) u, w, dist_2D_f3(i, j)/w
			end do
		end do
		close(1)

		open(1, file = 'characteristic_GD_' // name //'.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = n_H, M1, rho1, p1, V1, dL, dR"

		do i = 1, GD_par_n_points
			if (GD_mas_X(i) >= 0.0) then 
				ro1 = -1.0/GD_mas_V(i)
				V1 = GD_mas_V(i)
				p1 = GD_mas_P(i)
				M1 = dabs(V1)/sqrt(ggg * p1/ro1)
				EXIT
			end if
		end do

		sred = (-1.0/GD_mas_V(1) -1.0/GD_mas_V(GD_par_n_points))/2.0
		ro_start = -1.0/GD_mas_V(1)
		do i = 1, GD_par_n_points
			ro = -1.0/GD_mas_V(i)
			if(dabs(ro - ro_start) * 100.0/ro_start >= 3.0) then
				dL = dabs(GD_mas_X(i))
				EXIT
			end if
		end do

		ro_start = -1.0/GD_mas_V(GD_par_n_points)
		do i = GD_par_n_points, 1, -1
			ro = -1.0/GD_mas_V(i)
			if(dabs(ro - ro_start) * 100.0/ro_start >= 3.0) then
				dR = dabs(GD_mas_X(i))
				EXIT
			end if
		end do



		write(1,*) par_n_H, M1, ro1, p1, V1, dL, dR

		close(1)



		open(1, file = 'shock_GD_' // name //'.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = v, Q2, Q3, QQ"

		do i = 1, shock_n
			write(1,*) shock_v(i), shock_Q2(i), shock_Q3(i), ggg * shock_v(i) * shock_Q2(i) - shock_Q3(i) * (ggg - 1)
		end do
		close(1)

	end subroutine Print_GD

	subroutine Save_all(step_)
		integer(4) :: step
		integer(4), intent(in), optional :: step_
		character(len=3) :: name

		step = 0
		if(PRESENT(step_)) step = step_

		write(unit=name,fmt='(i3.3)') step

		open(1, file = "Save_all_" // name // ".bin", FORM = 'BINARY')

		write(1) GD_par_n_points, par_L, par_R, Mach_inf, ddx

		write(1) GD_mas_X, GD_mas_Q2, GD_mas_Q3, GD_mas_V, GD_mas_P, mas_n_H, mas_v_H, mas_T_H

		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
		write(1) 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1

		close(1)

	end subroutine Save_all

	subroutine Read_all(step_)
		integer(4) :: step, N, i
		integer(4), intent(in), optional :: step_
		character(len=3) :: name
		logical :: exists

		step = 0
		if(PRESENT(step_)) step = step_

		write(unit=name,fmt='(i3.3)') step

		inquire(file= "Save_all_" // name // ".bin", exist=exists)
		if (exists == .False.) then
            print*, "net faila 1898 tfgdhfwy4rfetrgfw4rwter!!!"
            pause
            STOP "net faila!!!"
        end if

		open(11, file = "Save_all_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")

		read(11) GD_par_n_points, par_L, par_R, Mach_inf, ddx

		print*, "GD_par_n_points", GD_par_n_points
		print*, "par_L", par_L
		print*, "par_R", par_R
		print*, "Mach_inf", Mach_inf

		call Set_GD()
		call M_K_Set()
		call Get_sensor_sdvig(0)

		read(11) GD_mas_X, GD_mas_Q2, GD_mas_Q3, GD_mas_V, GD_mas_P, mas_n_H

		close(11)

		GD_mas_Q2_do = GD_mas_Q2
		GD_mas_Q3_do = GD_mas_Q3

	end subroutine Read_all

	subroutine Print_MK()
		integer(4) :: i

		open(1, file = 'Print_MK.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = X, n_H"

		do i = 1, GD_par_n_points
			write(1,*) GD_mas_X(i), mas_n_H(i)
		end do

		close(1)

	end subroutine Print_MK

	subroutine Start_right_GD()
		integer :: i, N1
		real(8) :: dx, Mach, a1, a2

		GD_mas_V(GD_par_n_points) = par_Velosity_inf
		GD_mas_P(GD_par_n_points) = 1.0/(Mach_inf**2 * ggg)

		GD_mas_V(GD_par_n_points - 1) = par_Velosity_inf
		GD_mas_P(GD_par_n_points - 1) = 1.0/(Mach_inf**2 * ggg)
		
		do i = GD_par_n_points - 1, 1, -1
			N1 = i - 1
			if(GD_mas_X(i-1) < 0.0) EXIT
			Mach = sqrt( dabs(GD_mas_V(i))/GD_mas_P(i)/GD_ggg)
			!print*, Mach
			if(Mach < 1.00001) then
				N1 = i + 1
				print*, "Shockless tranzition"
				print*, Mach, GD_ggg * GD_mas_V(i)* GD_mas_Q2(i) - GD_mas_Q3(i) * (GD_ggg - 1.0)
				STOP "STOPPPPPPPPPPPPPPPPPPPPPPPPPP"
				EXIT
			end if
			
			dx = GD_mas_X(i) - GD_mas_X(i-1)
			GD_mas_V(i-1) = GD_mas_V(i) + dx * (GD_ggg * GD_mas_V(i)* GD_mas_Q2(i) - GD_mas_Q3(i) * (GD_ggg - 1.0)) /(GD_mas_V(i) + GD_ggg * GD_mas_P(i))
			
			GD_mas_P(i-1) = GD_mas_P(i) - dx * (GD_mas_Q2(i) - (GD_ggg * GD_mas_V(i)* GD_mas_Q2(i) - GD_mas_Q3(i) * (GD_ggg - 1.0)) /(GD_mas_V(i) + GD_ggg * GD_mas_P(i)))

			! print*, GD_mas_V(i-1), GD_mas_P(i-1)
		end do
		
		Mach = sqrt( dabs(GD_mas_V(N1 + 1))/GD_mas_P(N1 + 1)/GD_ggg)
		
		a1 = (GD_ggg + 1) * Mach**2 / (  (GD_ggg - 1) * Mach**2 + 2.0 )
		a2 = 2 * GD_ggg * Mach**2 / (GD_ggg + 1) - (GD_ggg - 1)/(GD_ggg + 1)
		
		GD_mas_P(N1) = GD_mas_P(N1 + 1) * a2
		GD_mas_V(N1) = GD_mas_V(N1 + 1) / a1
		
		print*, "Mach do = ", Mach 
		Mach = sqrt( dabs(GD_mas_V(N1))/GD_mas_P(N1)/GD_ggg)
	
		print*, "Mach = ", Mach 
		
		do i = N1, 2, -1
			dx = GD_mas_X(i) - GD_mas_X(i-1)
			GD_mas_V(i-1) = GD_mas_V(i) + dx * (GD_ggg * GD_mas_V(i)* GD_mas_Q2(i) - GD_mas_Q3(i) * (GD_ggg - 1.0)) & 
				/(GD_mas_V(i) + GD_ggg * GD_mas_P(i))
			
			GD_mas_P(i-1) = GD_mas_P(i) - dx * (GD_mas_Q2(i) - (GD_ggg * GD_mas_V(i)* GD_mas_Q2(i) - GD_mas_Q3(i) * (GD_ggg - 1.0)) & 
				/(GD_mas_V(i) + GD_ggg * GD_mas_P(i)))
		end do
		
		GD_mas_P(2) = GD_mas_P(3)
		GD_mas_V(2) = GD_mas_V(3)
		
		GD_mas_P(1) = GD_mas_P(2)
		GD_mas_V(1) = GD_mas_V(2)

		
	end subroutine Start_right_GD

	subroutine MK_ADD_MOMENT2(potok, cell, time, mu, VV, Source)
		integer(4), intent(in) :: cell, potok
		real(8), intent(in) :: time, VV(3), mu
		logical, intent(in) :: Source
		real(8) :: cp, vx, vy, vz, u, u1, u2, u3, skalar, uz, uz_M, uz_E, k1, k2, k3
		integer(4) :: i

		! Моменты функции распределения
		call omp_set_lock(mas_lock(cell))
		mas_n_H(cell) = mas_n_H(cell) + time * mu
		mas_v_H(cell) = mas_v_H(cell) + time * mu * VV(1)
		mas_T_H(cell) = mas_T_H(cell) + time * mu * kvv(VV(1), VV(2), VV(3))

		if(.True.) then
			!call Dist_funk(VV(1), time * mu, int((cell - 1)/1) + 1)
			call Dist_funk(VV(1), time * mu, cell)
		end if

		if(.False.) then
			if(cell == cell_f1) then
				call Dist_funk2D(VV(1), sqrt(VV(2)**2 + VV(3)**2), time * mu, 1)
			end if

			if(cell == cell_f2) then
				call Dist_funk2D(VV(1), sqrt(VV(2)**2 + VV(3)**2), time * mu, 2)
			end if

			if(cell == cell_f3) then
				call Dist_funk2D(VV(1), sqrt(VV(2)**2 + VV(3)**2), time * mu, 3)
			end if
		end if

		call omp_unset_lock(mas_lock(cell))

		! Источники
		if(Source == .True.) then
			cp = sqrt(-GD_mas_P(cell) * GD_mas_V(cell))
			vx = GD_mas_V(cell)
			vy = 0.0
			vz = 0.0

			u = sqrt(kvv(VV(1) - vx, VV(2) - vy, VV(3) - vz))
			u1 =  vx - VV(1)
			u2 =  vy - VV(2)
			u3 =  vz - VV(3)
			skalar = VV(1) * u1 + VV(2) * u2 + VV(3) * u3

			if (u / cp > 7.0) then
				uz = MK_Velosity_1(u, cp)
				uz_M = MK_Velosity_2(u, cp)/ (uz * (cp**2) * cp * par_pi * par_sqrtpi)
				uz_E = MK_Velosity_3(u, cp)
				
				call omp_set_lock(mas_lock(cell))
				GD_mas_Q2(cell) = GD_mas_Q2(cell) - mu * uz_M * u1 / u
				GD_mas_Q3(cell) = GD_mas_Q3(cell) + &
					mu * (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E / uz) - uz_M * skalar / u)
				call omp_unset_lock(mas_lock(cell))
			else
				k1 = MK_int_1(u, cp)
				k2 = MK_int_2(u, cp)
				k3 = MK_int_3(u, cp)
				
				call omp_set_lock(mas_lock(cell))
				GD_mas_Q2(cell) = GD_mas_Q2(cell) + mu * (k2/k1) * u1 / u
				GD_mas_Q3(cell) = GD_mas_Q3(cell) + mu * (-0.5 * k3/k1 + k2/k1 * skalar / u)
				call omp_unset_lock(mas_lock(cell))
			end if


			if (.False.) then!(cell == shock_cell .or. cell == shock_cell + 1) then  !(.False.) then!
				do i = 1, shock_n
					cp = sqrt(shock_v(i)**2 / ggg)
					vx = -shock_v(i)
					vy = 0.0
					vz = 0.0

					u = sqrt(kvv(VV(1) - vx, VV(2) - vy, VV(3) - vz))
					u1 =  vx - VV(1)
					u2 =  vy - VV(2)
					u3 =  vz - VV(3)
					skalar = VV(1) * u1 + VV(2) * u2 + VV(3) * u3

					if (u / cp > 7.0) then
						uz = MK_Velosity_1(u, cp)
						uz_M = MK_Velosity_2(u, cp)/ (uz * (cp**2) * cp * par_pi * par_sqrtpi)
						uz_E = MK_Velosity_3(u, cp)
						
						call omp_set_lock(mas_lock(i))
						shock_Q2(i) = shock_Q2(i) - mu * uz_M * u1 / u
						shock_Q3(i) = shock_Q3(i) + &
							mu * (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E / uz) - uz_M * skalar / u)
						call omp_unset_lock(mas_lock(i))
					else
						k1 = MK_int_1(u, cp)
						k2 = MK_int_2(u, cp)
						k3 = MK_int_3(u, cp)
						
						call omp_set_lock(mas_lock(i))
						shock_Q2(i) = shock_Q2(i) + mu * (k2/k1) * u1 / u
						shock_Q3(i) = shock_Q3(i) + mu * (-0.5 * k3/k1 + k2/k1 * skalar / u)
						call omp_unset_lock(mas_lock(i))
					end if
				end do
			end if
		end if
	end subroutine MK_ADD_MOMENT2

	subroutine Proverka_source()
		real(8) :: cp, vx, vy, vz, VV(3), u, u1, u2, u3, skalar, uz, uz_M, uz_E, k1, k2, k3, kkk1, kkk2, kkk
		integer(4) :: i
		cp = 0.3
		vx = 0.001
		vy = 0.0
		vz = 0.0

		VV(2) = 0.0
		VV(3) = -0.2


		open(1, file = 'Proverka_source.txt')
		write(1,*) "TITLE = 'HP'  VARIABLES = u, u1, u2, u3, cp, skalar, M2, MM2, M3, I2, I3, MM3"


		do i = 1, 5000
			VV(1) = i * 0.01

			u = sqrt(kvv(VV(1) - vx, VV(2) - vy, VV(3) - vz))
			u1 =  vx - VV(1)
			u2 =  vy - VV(2)
			u3 =  vz - VV(3)
			skalar = VV(1) * u1 + VV(2) * u2 + VV(3) * u3

			uz = MK_Velosity_1(u, cp)
			uz_M = MK_Velosity_2(u, cp)/ (uz * (cp**2) * cp * par_pi * par_sqrtpi)
			uz_E = MK_Velosity_3(u, cp)
			!print*, "uz_M, uz_E", uz, uz_M, uz_E
			!pause
		
			k1 = MK_int_1(u, cp)
			k2 = MK_int_2(u, cp)
			k3 = MK_int_3(u, cp)

			kkk1 = 1.0
			kkk2 = 1.0

			kkk = MK_sigma(uz_E)/MK_sigma(uz)
			if(u/cp > 0.8) kkk = (1.0 + kkk)/2.0

			if(u/cp < 0.4) then
				kkk1 =  MK_sigma(uz_E)/MK_sigma(uz)
			else
				kkk2 =  MK_sigma(uz_M)/MK_sigma(uz)
			end if
			
			write(1,*) u/cp, u1, u2, u3, cp, skalar, -uz_M * u1 / u , -uz_M * u1 / u * kkk, (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E / uz) - uz_M * skalar / u), &
				(k2/k1) * u1 / u, (-0.5 * k3/k1 + k2/k1 * skalar / u),&
					 (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E * kkk1 / (uz)) - uz_M * skalar / u * kkk2)
			
			
		end do

		close(1)
	end subroutine Proverka_source

	subroutine MK_FLY(potok, VVx0, VVy0, VVz0, x0, cell_0, mu)
		! cell - текущая ячейка частицы
		real(8), intent(in) :: VVx0, VVy0, VVz0, x0, mu
		integer(4), intent(in) :: cell_0, potok

		real(8) :: ro, p, vx, cp, vy, vz, time, uz, Wr_, Wthe_, Wphi_, I_do, x, VVx, VVy, VVz, KSI
		real(8) :: L, R, time1, time2, nu_ex, lenght
		real(8) :: skalar, u, u1, u2, u3, sig, II, t_ex, x_ex, ksi1, VV(3)
		integer(4) :: next, cell

		I_do = 0.0
		cell = cell_0
		x = x0
		VVx = VVx0
		VVy = VVy0
		VVz = VVz0
		call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
		KSI = -log(1.0 - ksi1)

		do while (.True.)
			L = GD_mas_X(cell) - ddx/2.0
			R = GD_mas_X(cell) + ddx/2.0
			
			time1 = (L - x)/VVx
			time2 = (R - x)/VVx


			if(time1 > 0) then
				next = cell - 1
				time = time1
			else
				next = cell + 1
				time = time2
			end if
			
			time = time * 1.00001

			if(time <= 0.0) then
				print*, time
				print*, time1
				print*, time2
				print*, "__________________"
				print*, L, x, R
				print*, "__________________"
				print*, cell, next
				print*, "__________________"
				print*, "__________________"
				x = (L + R)/2.0
				CYCLE
			end if

			vx = -0.963449!! GD_mas_V(cell)
			vy = 0.0
			vz = 0.0
			ro = 1.0!! -1.0/vx
			p = 3.06593096_8!!GD_mas_P(cell)
			cp = sqrt(p/ro)
			u = sqrt(kvv(VVx - vx, VVy - vy, VVz - vz))
			u1 =  vx - VVx
			u2 =  vy - VVy
			u3 =  vz - VVz
			skalar = VVx * u1 + VVy * u2 + VVz * u3
			II = I_do
			lenght = sqrt(kvv(VVx * time, VVy * time, VVz * time))

			if (.True.) then!(u / cp > 7.0) then
				uz = MK_Velosity_1(u, cp)
				nu_ex = (ro * uz * MK_sigma(uz))
			else
				nu_ex = (ro * MK_int_1(u, cp))
			end if

			sig = sqrt(kvv(VVx, VVy, VVz))/nu_ex
			II = II + lenght / sig

			if (II < KSI)  then !(.True.) then!
				I_do = II

				! Моменты
				!mas_n_H(cell) = mas_n_H(cell) + time * mu
				VV(1) = VVx
				VV(2) = VVy
				VV(3) = VVz
				call MK_ADD_MOMENT2(potok, cell, time, mu, VV, Source = .False.)
				!print*, time, mu
				!pause

				!if(cell == 1 .and. next < 1 .and. VVx <= 0.0) call Dist_funk(VVx, mu)

				cell = next
				x = x + time * VVx
				if(next < 1 .or. next > GD_par_n_points) EXIT
				CYCLE
			else
				!t_ex = (II - KSI) * sig / sqrt(kvv(VVx, VVy, VVz))
				t_ex = (KSI - I_do) * sig / sqrt(kvv(VVx, VVy, VVz))
				if(time <= t_ex) then
					print*, "ERROR oi09u39y7808yqupwtu8q93vt"
					print*, time
					print*, t_ex
					print*, "_______________________"
					print*, II
					print*, KSI
					print*, "_______________________"
					pause
				end if

				I_do = 0.0
				x_ex = x + t_ex * VVx

				!!if( (x_ex - (-2.9)) * (x - (-2.9)) < 0.0 ) call Dist_funk(VVx, mu)

				if(x_ex > R .or. x_ex < L) then
					t_ex = t_ex * 0.999
					x_ex = x + t_ex * VVx
				end if

				! Моменты
				!mas_n_H(cell) = mas_n_H(cell) + t_ex * mu
				VV(1) = VVx
				VV(2) = VVy
				VV(3) = VVz
				call MK_ADD_MOMENT2(potok, cell, t_ex, mu, VV, Source = .True.)


				call M_K_Change_Velosity4(potok, vx/cp, vy/cp, vz/cp, &
					VVx/cp, VVy/cp, VVz/cp, Wr_, Wthe_, Wphi_, cp)
				VVx = Wr_ * cp
				VVy = Wthe_ * cp
				VVz = Wphi_ * cp
				
				!print*, VVx, VVy, VVz
				!pause

				call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
				KSI = -log(1.0 - ksi1)
				x = x_ex
				CYCLE
			end if

		end do

	end subroutine MK_FLY

	subroutine Start_MK()
		! Запус Монте-Карло
		USE Monte_Karlo
		real(8) :: Vx, Vy, Vz, S, cp, S1, S2, Ux, mu1, mu2, swap
		real(8) :: X, u, u1, u2, u3, skalar, ksi1, uz, L, nu_ex, time, lenght, Wr_, Wthe_, Wphi_, ro
		integer(4) :: potok, i, Num, no, Num1, Num2, step
		
		potok = 1
		Num1 = 10000000! 1000000 * 30! 30! 20! 25! * 32 * 4! * 32! * 32 * 2
		Num2 = 0! 1000000 * 10! 10! * 16 * 4! * 16! * 40
		Num = Num1 + Num2

		mas_n_H = 0.0
		mas_v_H = 0.0
		mas_T_H = 0.0
		GD_mas_Q2 = 0.0
		GD_mas_Q3 = 0.0
		dist_f = 0.0
		dist_2D_f1 = 0.0
		dist_2D_f2 = 0.0
		dist_2D_f3 = 0.0

		! Сначала для вылета справа
		Ux = -3.0!! par_Velosity_inf
		cp = 1.0! sqrt(-GD_mas_P(GD_par_n_points) * GD_mas_V(GD_par_n_points))
		S1 = 0.5 * (cp * exp(-Ux**2/cp**2)/sqrtpi - Ux + Ux * erf(Ux/cp))
		!print*, cp, Ux, S1 

		! Для вылета слева
		Ux = GD_mas_V(1)
		cp = sqrt(-GD_mas_P(1) * GD_mas_V(1))
		S2 = 0.0!!0.5 * (cp * exp(-Ux**2/cp**2)/sqrtpi + Ux + Ux * erf(Ux/cp)) * (-1.0/Ux)

		!print*, cp, Ux , S2
		
		S = S1 + S2
		mu1 = (S1/S) * (1.0 * Num/Num1)
		mu2 = (S2/S) * (1.0 * Num/Num2)

		print*, "Mu = ", mu1, mu2
		
		Ux = -3.0!! par_Velosity_inf
		cp = 1.0!
		step = 1

		call omp_set_num_threads(par_n_potok)
		! !$omp parallel
		! !$omp do private(X, Vx, Vy, Vz, potok)
		do i = 1, Num1
			if(mod(i, 100000) == 0) print*, "I = ", i
			potok = (omp_get_thread_num() + 1) 

			! !$omp critical
			! 	if(mod(step, 300000) == 0) print*, "START = ", step
			! 	step = step + 1
			! !$omp end critical
			

			X = par_R - 0.000001


			call MK_Velosity_initial2(potok, Vx, Vy, Vz, Ux/cp)
			Vx = Vx * cp
			Vy = Vy * cp
			Vz = Vz * cp
			!!call MK_FLY(potok, Vx, Vy, Vz, X, GD_par_n_points, mu1)

			!! Сразу накопим функцию распределения
			call Dist_funk(Vx, -1.0/Vx, 1)
		end do
		! !$omp end do

		! !$omp barrier

		! !$omp single
		!print*, "Ux = ", Ux, cp
		Ux = GD_mas_V(1)
		cp = sqrt(-GD_mas_P(1) * GD_mas_V(1))
		step = 1
		! !$omp end single

		! !$omp barrier

		! !$omp do private(X, Vx, Vy, Vz, potok)
		do i = 1, Num2
			potok = (omp_get_thread_num() + 1) 
			! !$omp critical
			! 	if(mod(step, 100000) == 0) print*, "START = ", step
			! 	step = step + 1
			! !$omp end critical
			X = par_L + 0.000001
			call MK_Velosity_initial3(potok, Vx, Vy, Vz, Ux/cp)
			Vx = Vx * cp
			Vy = Vy * cp
			Vz = Vz * cp
			!call Dist_funk(Vx, mu2)
			call MK_FLY(potok, Vx, Vy, Vz, X, 1, mu2)
		end do
		! !$omp end do


		! !$omp end parallel

		!print*, "Ux = ", Ux, cp

		!! ТУТ НИЖЕ ДОБАВИТЬ
		no = Num!! * (GD_mas_X(2) - GD_mas_X(1))  ! ДОБАВИТЬ
		dist_f = dist_f * S/(no)/((vR - vL)/vN)
		dist_2D_f1 = dist_2D_f1 * S/(no)/((vR - vL)/vN)/((wR - wL)/vN)
		dist_2D_f2 = dist_2D_f2 * S/(no)/((vR - vL)/vN)/((wR - wL)/vN)
		dist_2D_f3 = dist_2D_f3 * S/(no)/((vR - vL)/vN)/((wR - wL)/vN)

		mas_n_H = S * mas_n_H/(no)
		mas_v_H = S * mas_v_H/(no)
		mas_T_H = S * mas_T_H/(no)

		do i = 1, size(mas_n_H)
			mas_v_H(i) = mas_v_H(i)/mas_n_H(i)
			mas_T_H(i) = (1.0/3.0) * (mas_T_H(i)/mas_n_H(i) - mas_v_H(i)**2)
		end do

		mas_n_H = par_n_H * mas_n_H
		GD_mas_Q2 = par_n_H * S * GD_mas_Q2/(no)
		GD_mas_Q3 = par_n_H * S * GD_mas_Q3/(no)

		shock_Q2 = par_n_H * S * shock_Q2/(2.0 * no)
		shock_Q3 = par_n_H * S * shock_Q3/(2.0 * no)

		do i = 1, size(GD_mas_Q2)
			! swap = GD_mas_Q2_do(i)
			! GD_mas_Q2_do(i) = GD_mas_Q2(i)
			! GD_mas_Q2(i) = (GD_mas_Q2(i) + swap)/2.0

			GD_mas_Q2_do(i) = (GD_mas_Q2_do(i) + GD_mas_Q2(i))/2.0   !! БЫЛИ ЭТИ 
			GD_mas_Q2(i) = GD_mas_Q2_do(i)


			! swap = GD_mas_Q3_do(i)
			! GD_mas_Q3_do(i) = GD_mas_Q3(i)
			! GD_mas_Q3(i) = (GD_mas_Q3(i) + swap)/2.0

			GD_mas_Q3_do(i) = (GD_mas_Q3_do(i) + GD_mas_Q3(i))/2.0
			GD_mas_Q3(i) = GD_mas_Q3_do(i)

			if(GD_mas_X(i) < -4.0) then
				GD_mas_Q3(i) = 0.0
				GD_mas_Q2(i) = 0.0
				GD_mas_Q3_do(i) = 0.0
				GD_mas_Q2_do(i) = 0.0
			end if
		end do

		!! dist_f(i) * S/Num/((vR - vL)/vN)

		!print*, "END 1 = ", S, no, Num, (GD_mas_X(2) - GD_mas_X(1)), GD_mas_Q2(100)
		!print*, "END 2 = ", mas_n_H(1), mas_n_H(10)

		open(1, file = "Dist_funk.txt")
	
		write(1, *)  "TITLE = 'HP'  VARIABLES = Vx, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, fp, ftest"
		do i = 1, vN
			Ux = -0.963449!GD_mas_V(1)
			cp = 1.75098!sqrt(-GD_mas_P(1) * GD_mas_V(1))
			ro = 1.0! -1.0/Ux
			u = (vL + (i - 0.5) * (vR - vL)/vN)
			if(u > 0) then
				WRITE (1, *) u, dist_f(i, 1), dist_f(i, 50), dist_f(i, 100), dist_f(i, 150), dist_f(i, 200), dist_f(i, 250), dist_f(i, 300),&
				dist_f(i, 350), dist_f(i, 400), dist_f(i, 450), ro * 1.0/(cp * par_sqrtpi) * exp(-(u - Ux)**2 / (cp)**2), 0.0
			else
				WRITE (1, *) u, dist_f(i, 1), dist_f(i, 50), dist_f(i, 100), dist_f(i, 150), dist_f(i, 200), dist_f(i, 250), dist_f(i, 300),&
				dist_f(i, 350), dist_f(i, 400), dist_f(i, 450), ro * 1.0/(cp * par_sqrtpi) * exp(-(u - Ux)**2 / (cp)**2), exp(-(u + 3.0)**2)/1.77243
			end if
		end do
		

	end subroutine Start_MK
	
end module GD
	
subroutine Luch2_do()
	USE Monte_Karlo
    implicit none
	
	real(8) :: Vx, Vy, Vz, S, cp
	real(8) :: X, u, u1, u2, u3, skalar, ksi1, uz, L, nu_ex, time, lenght, Wr_, Wthe_, Wphi_
	integer(4) :: potok, i, Num
	
	! potok = 1
	! Num = 10000
	! S = 0.5 * (exp(-par_Velosity_inf**2)/par_sqrtpi - par_Velosity_inf * erfc(par_Velosity_inf))
	! !S = 0.5 * (exp(-1.0**2)/par_sqrtpi - (-1.0) * erfc(-1.0))
	
	! call M_K_Set()
	! call Get_sensor_sdvig(0)
	
	! do i = 1, Num
	! 	X = par_R
	! 	call MK_Velosity_initial2(potok, Vx, Vy, Vz)
	! 	!print*, "V = ", i, Vx, Vy, Vz
		
	! 	do while (.True.)
	! 		! пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ пїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ пїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ  ****************************************************************************************
	! 		u = sqrt(kvv(Vx - par_vx, Vy - par_vy, Vz - par_vz))
	! 		u1 =  par_vx - Vx
	! 		u2 =  par_vy - Vy
	! 		u3 =  par_vz - Vz
	! 		skalar = Vx * u1 + Vy * u2 + Vz * u3
			
	! 		if (.True.) then!(u / par_cp > 7.0) then
	! 			!print*, "tut"
	! 			uz = MK_Velosity_1(u, par_cp)
	! 			nu_ex = (par_ro * uz * MK_sigma(uz))
	! 		else
	! 			nu_ex = (par_ro * MK_int_1(u, par_cp))        ! пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	! 			!print*, "2 =", nu_ex, par_ro, u, par_cp
	! 		end if
		
	! 		L = sqrt(Vx**2 + Vy**2 + Vz**2)/nu_ex
	! 		call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi1)
	! 		lenght = -L * log(ksi1)
	! 		time = lenght/sqrt(Vx**2 + Vy**2 + Vz**2)
			
	! 		if ((x0 - X)/Vx > 0.0 .and. (x0 - X)/Vx < time) then
	! 			call Dist_funk(Vx)
	! 		end if
		
	! 		X = X + time * Vx
			
			
	! 		!print*, "x = ", X
	! 		!print*, "V = ", i, Vx, Vy, Vz, time, lenght, ksi1, nu_ex
	! 		!pause
	! 		if(X < par_L .or. X > par_R) EXIT
		
	! 		call M_K_Change_Velosity4(potok, par_vx/par_cp, par_vy/par_cp, par_vz/par_cp, &
	! 			Vx/par_cp, Vy/par_cp, Vz/par_cp, Wr_, Wthe_, Wphi_, par_cp)
	! 		Vx = Wr_ * par_cp
	! 		Vy = Wthe_ * par_cp
	! 		Vz = Wphi_ * par_cp
	! 	end do
		
	! end do
	

    ! Variables

    ! Body of Luch2
    ! print *, 'Hello World'
	
	! open(1, file = "Dist_funk.txt")
	! open(2, file = "Dist_funk_0.txt")
	! open(3, file = "Dist_funk_1.txt")
	! open(4, file = "sravnenie secheniy.txt")
	
	! do i = 1, vN
	! 	u = (vL + (i + 0.5) * (vR - vL)/vN)
	! 	cp = 0.386209
	! 	WRITE (1, *) u, dist_f(i) * S/Num/((vR - vL)/vN)
	! 	!WRITE (2, *) u, 1.0/(par_sqrtpi) * exp(-(u - par_Velosity_inf)**2)
	! 	cp = 0.55672027
	! 	WRITE (3, *) u, par_ro * 1.0/(cp * par_sqrtpi) * exp(-(u - par_vx/par_rat_v)**2 / (cp)**2)
	! end do
	
	! do i = 1, 1000
	! 	u = (vL + (i + 0.5) * (vR - vL)/1000)
	! 	WRITE (2, *) u, 1.0/(par_sqrtpi) * exp(-(u - par_Velosity_inf)**2)
	! end do
	
	
	! do i = 1, vN
	! 	u = i * 0.1
	! 	uz = MK_Velosity_1(u, par_cp)
	! 	WRITE (4, *) u, uz * MK_sigma(uz), MK_int_1(u, par_cp)
	! end do
	
		
	
	! close(1)
	! close(2)
	! close(3)
	! close(4)
	
	
	
	pause

end subroutine Luch2_do


program Luch2
	USE GD
	integer(4) :: step, kk

	call Set_GD()
	call M_K_Set()
	call Get_sensor_sdvig(0)

	! call Read_all(203)

	call Print_GD()

	! pause
	! Mach_inf = 1.3_8
	! par_n_H = 1.95_8
	kk = 0
	do step = 1, 1
		kk = kk + 1
		print*, "M-K ", step, "nH ", par_n_H
		call Start_MK()
		print*, "G-D", step
		! call Start_right_GD()

		call Save_all(step)
		call Print_GD(step)

		if(mod(kk, 5) == 0) then
			par_n_H = par_n_H + 0.025
			print*, "MENYAEM nH", par_n_H
		end if

		! if(par_n_H < 4.99_8) then
		! 	par_n_H = par_n_H + 0.07
		! else
		! 	par_n_H = 5.0_8
		! end if

	end do

	print*, "END"

end program Luch2