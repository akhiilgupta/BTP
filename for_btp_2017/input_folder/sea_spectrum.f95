! program to generate sea spectrum of irregular waves

program sea_spec
    implicit none
	integer, parameter :: max_n = 100
    real, parameter :: g = 9.81, h_sig = 3, U = 10, pi = 3.14159265358979323846
	real :: we(max_n), w(max_n), heave(max_n), pitch(max_n), surge(max_n), area_curve_h, area_curve_p
    real :: sway(max_n), roll(max_n), yaw(max_n), sw(max_n), sz(max_n), junk, heave_res_spec(max_n), beta, V, pitch_res_spec(max_n)
	integer :: i, j, io
    
    write(*,*) 'Enter the heading angle and speed at which encounter spectrum is to be calculated'
    read(*,*) beta, V
    
    open(unit = 1, file = "out16", status='old', action='read')
    do i = 1,8
        read(1,*)
    end do
    
    i = 1
    
    do
        read(1,*,IOSTAT = io) we(i),w(i),junk,junk,surge(i),junk,sway(i),junk,heave(i),junk,roll(i),junk,pitch(i),junk,yaw(i)
        if (io .lt. 0.) exit
        ! write(*,*) surge(i), sway(i)
        i = i+1
    end do
    
    do j = 1, i-1
        we(j) = we(j)*sqrt(9.8060/100.0)
        w(j) = w(j)*sqrt(9.8060/100.0)
    end do
    
    open(unit = 2, file = 'sea_spec')          
    
    write(2,*) 'w       sw      we      sz      heave   heave_res   pitch   pitch_res'
    write(2,*)
    
    do j = 1, i-1
        sw(j) = .78*exp(-0.13/(w(j)**4))/(w(j)**5)
        sz(j) = sw(j)/(1-2*w(j)*V*.514*cos(beta*pi/180)/g)
        heave_res_spec(j) = ((heave(j))**2)*sz(j)
        pitch_res_spec(j) = (pitch(j)**2)*sz(j)
        write(2,55) w(j),sw(j),we(j),sz(j),heave(j),heave_res_spec(j),pitch(j),pitch_res_spec(j)
    end do
    
    55	format(8(f6.4, 3X))
    
    do j = 2, i-1
        area_curve_h = area_curve_h + abs(we(j)-we(j-1))*heave_res_spec(j)
        area_curve_p = area_curve_p + abs(we(j)-we(j-1))*pitch_res_spec(j)
    end do
    
    write(*,*) 2*sqrt(area_curve_h), 2*sqrt(area_curve_p)
    
end program sea_spec