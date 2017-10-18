! program to change length and create new input geometry

program new_length 
	implicit none
	real :: length_overall, stn_aft, stn_fwd, new_lengthoa, new_breadth, jt
	real, allocatable, dimension(:) :: stn_num
	integer, allocatable, dimension(:) :: num_points
	integer, parameter :: max_n_points = 100, max_n_stn = 50
	real :: stn_dis, yjt(max_n_stn,max_n_points), zjt(max_n_stn, max_n_points)
	real :: ini_breadth = 10.0, shr_factor1, shr_factor2, ini_depth = 6.25, new_depth
	integer :: i, j, n_stn
	
	open(unit = 1, file = "input_geometry")
	
	read(1,*) n_stn,length_overall,stn_aft,stn_fwd
	stn_dis = length_overall/(stn_fwd-stn_aft)
	
	allocate(stn_num(n_stn), num_points(n_stn))
	
	do i = 1,n_stn
		read(1,*) stn_num(i),num_points(i)
		do j = 1,num_points(i)
			read(1,*) jt,yjt(i,j),zjt(i,j)
		enddo
	enddo
	
	close(1)
	
	print *, "enter new desired length, breadth and depth of the ship"
	read *, new_lengthoa, new_breadth, new_depth
	
	shr_factor1 = new_breadth/ini_breadth
	shr_factor2 = new_depth/ini_depth
	
	open(unit = 2, file = "input_folder/new_input_geometry")
	
	write(2,*) n_stn, new_lengthoa, stn_aft, stn_fwd
	
	do i = 1,n_stn
		write(2,*) stn_num(i), num_points(i)
		do j = 1,num_points(i)
			write(2,*) j, yjt(i,j)*shr_factor1, zjt(i,j)*shr_factor2
		enddo
	enddo
	
	
end program new_length