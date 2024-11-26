module Math_Libraries
contains
recursive subroutine Fast_Sort(A,sort_type)
        real :: A(:)
        character(len=*),optional :: sort_type
        integer :: head,tail,i,s
        real :: temp
        real :: pivot
        if (.not. present(sort_type)) then
                sort_type='asc'
        end if
        head=0
        tail=size(A)+1
        pivot=A(1)
loop1:  do while (.true.)
                tail=tail-1
                head=head+1
tail_loop:      do while (.true.)
                        if (A(tail) <= pivot) then
                                exit tail_loop
                        end if
                        tail=tail-1
                end do tail_loop

head_loop:      do while (.true.) 
                        if (A(head) > pivot) then
                                exit head_loop
                        end if
                        head=head+1
                        if (head > size(A)) then
                                exit head_loop
                        end if
                end do head_loop

                if (head < tail) then
                        temp=A(head)
                        A(head)=A(tail)
                        A(tail)=temp

                else if (head == tail) then
                        temp=A(head)
                        A(head)=A(1)
                        A(1)=temp

                        exit loop1
                else if (head > tail) then
                        temp=A(tail)
                        A(tail)=A(1)
                        A(1)=temp

                        exit loop1
                end if
        end do loop1
                if (size(A(:tail-1)) > 1) then
                        call Fast_Sort(A(:tail-1))
                end if

                if (size(A(tail+1:)) > 1) then
                        call Fast_Sort(A(tail+1:))
                else 
                        return
                end if

        if (trim(adjustl(sort_type)) == 'desc') then
                call reverse_array(A)
        end if
        

contains
        
subroutine reverse_array(array)
        implicit none
        real :: array(:)
        integer :: s,i,temp

        s=size(array)

        do i=1,s/2
                temp=array(i)
                array(i)=array(s+1-i)
                array(s+1-i)=temp
        end do
end subroutine reverse_array

end subroutine Fast_Sort

function differentiate(x,y) result(d)
implicit none
real :: x(:),y(:)
real,dimension(size(x)-1) :: d
integer :: n,i,j,k,l,m

if (size(x)/=size(y)) then
    write(*,*) 'x size is: ',size(x),'y size is: ',size(y),' mismatch!'
    return
end if

n=size(x)

do i=1,n-1
    d(i)=(y(i+1)-y(i))/(x(i+1)-x(i))
end do

end function differentiate
end module Math_Libraries

program Amplitude
!$ use mpi_lib
use Math_Libraries
implicit none
type :: atom
        integer :: atom_type
        real :: x(3)
end type atom

integer :: i,j,k,l,m,n
integer :: alive
integer :: step
integer,dimension(5) :: atom_number
integer,allocatable :: indecies(:),t_x(:),t_y(:),t_z(:)
real,parameter :: cutoff_r=2.0
real,dimension(6) :: x_b,x_b_origin
real :: a,b,c,x0,y0,z0,time
real,dimension(3) :: temp
real,allocatable :: dr(:),dx(:,:),dy(:,:),dz(:,:),distance(:),d_x(:),d_y(:),d_z(:),input_step(:)
real,allocatable :: V_x(:),V_y(:),V_z(:),V(:),F_x(:),F_y(:),F_z(:),F(:),px(:),py(:),pz(:)!This means atoms vibration.
type(atom),allocatable :: R(:,:)
type(atom),allocatable :: R_origin(:)
step=0
atom_number(:)=0
!-----------------------------------------------------------------------|
!                                                                       |
!                               read file                               |
!                                                                       |
!-----------------------------------------------------------------------|
open(unit=100,file='dump.atom',status='old',action='read')
open(unit=200,file='dump1.atom',status='old',action='read')
open(unit=300,file='data/amplitude/amplitude-900k.dat',status='replace',action='write')
open(unit=400,file='data/frequency/frequency_raw-900k.dat',status='replace',action='write')
!-----------------------------------------------------------------------|
!                                                                       |
!                 confirm atom number and crystal constant              |
!                                                                       |
!-----------------------------------------------------------------------|                                                                      
do i=1,3
        read(100,*)
end do
read(100,*) atom_number(1)
read(100,*)
read(100,*) x_b(1),x_b(2)
read(100,*) x_b(3),x_b(4)
read(100,*) x_b(5),x_b(6)
a=x_b(2)-x_b(1)
b=x_b(4)-x_b(3)
c=x_b(6)-x_b(5)
rewind(100)
!-----------------------------------------------------------------------|
!                                                                       |
!                               confirm step                            |
!                                                                       |
!-----------------------------------------------------------------------|                                                                           
m=1
do while (.true.)
        do i=1,9
                read(100,*,iostat=alive)
                if (alive/=0) then
                        goto 55
                end if
        end do
        
        do i=1,atom_number(1)
                read(100,*) 
        end do
        step=step+1
 end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               read file                               |
!                                                                       |
!-----------------------------------------------------------------------|
55 continue
rewind(100)
allocate(R(atom_number(1),step))
allocate(R_origin(atom_number(1)))
do i=1,step
        do j=1,9
                read(100,*)
        end do

        do j=1,atom_number(1)
                read(100,*) n,R(n,i)%atom_type,R(n,i)%x(1),R(n,i)%x(2),R(n,i)%x(3)
                R(n,i)%x(1)=R(n,i)%x(1)-x_b(1)
                R(n,i)%x(2)=R(n,i)%x(2)-x_b(3)
                R(n,i)%x(3)=R(n,i)%x(3)-x_b(5)
                if (i==1) then
                        atom_number(1+R(n,i)%atom_type)=atom_number(1+R(n,i)%atom_type)+1
                end if
        end do
end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               read origin Li file                     |
!                                                                       |
!-----------------------------------------------------------------------|
do i=1,5
        read(200,*)
end do
read(200,*) x_b_origin(1),x_b_origin(2)
read(200,*) x_b_origin(3),x_b_origin(4)
read(200,*) x_b_origin(5),x_b_origin(6)
read(200,*) 
do i=1,atom_number(1)
        read(200,*) l,R_origin(l)%atom_type,R_origin(l)%x(1),R_origin(l)%x(2),R_origin(l)%x(3)
        R_origin(l)%x(1)=R_origin(l)%x(1)-x_b_origin(1)
        R_origin(l)%x(2)=R_origin(l)%x(2)-x_b_origin(3)
        R_origin(l)%x(3)=R_origin(l)%x(3)-x_b_origin(5)
end do
close(200)
!-----------------------------------------------------------------------|
!                                                                       |
!                     calclulate Amplitude of vibration                 |
!                                                                       |
!-----------------------------------------------------------------------|
allocate(dx(atom_number(2),step))
allocate(dy(atom_number(2),step))
allocate(dz(atom_number(2),step))
allocate(dr(atom_number(2)))
allocate(distance(atom_number(2)))
allocate(d_x(step-1))
allocate(d_y(step-1))
allocate(d_z(step-1))
allocate(V_x(atom_number(2)))
allocate(V_y(atom_number(2)))
allocate(V_z(atom_number(2)))
allocate(V(atom_number(2)))
allocate(F_x(atom_number(2)))
allocate(F_y(atom_number(2)))
allocate(F_z(atom_number(2)))
allocate(F(atom_number(2)))
allocate(t_x(step))
allocate(t_y(step))
allocate(t_z(step))
allocate(input_step(step))
allocate(px(atom_number(2)))
allocate(py(atom_number(2)))
allocate(pz(atom_number(2)))
V_x(:)=0.0
V_y(:)=0.0
V_z(:)=0.0
V(:)=0.0
F_x(:)=0.0
F_y(:)=0.0
F_z(:)=0.0
F(:)=0.0
t_x(:)=0
t_y(:)=0
t_z(:)=0
px(:)=0.0
py(:)=0.0
pz(:)=0.0
do i=1,step
    input_step(i)=i
end do

do i=1,atom_number(2)
    t_x(:)=0
    t_y(:)=0
    t_z(:)=0
    x0=R(i,1)%x(1)
    y0=R(i,1)%x(2)
    z0=R(i,1)%x(3)
    dx(i,1)=0.0
    dy(i,1)=0.0
    dz(i,1)=0.0
    do j=2,step
        temp(1)=R(i,j)%x(1)-x0
        temp(2)=R(i,j)%x(2)-y0
        temp(3)=R(i,j)%x(3)-z0

        !x:
        if (temp(1)>=0.5*a) then
            temp(1)=temp(1)-a
        else if (temp(1)<=-0.5*a) then
            temp(1)=a+temp(1)
        else
            continue
        end if

        !y:
        if (temp(2)>=0.5*b) then
            temp(2)=temp(2)-b
        else if (temp(2)<=-0.5*b) then
            temp(2)=b+temp(2)
        else
            continue
        end if

        !z:
        if (temp(3)>=0.5*c) then
            temp(3)=temp(3)-c
        else if (temp(3)<=-0.5*c) then
            temp(3)=c+temp(3)
        else
            continue
        end if

        if (abs(temp(1))>cutoff_r**(1.0/4.0) .or. &
            abs(temp(2))>cutoff_r**(1.0/4.0) .or. &
            abs(temp(3))>cutoff_r**(1.0/4.0)) then

            do k=1,atom_number(2)
                distance(k)=sqrt((R(i,j)%x(1)-R_origin(k)%x(1))**2+&
                (R(i,j)%x(2)-R_origin(k)%x(2))**2+&
                (R(i,j)%x(3)-R_origin(k)%x(3))**2)
            end do
            
            call Fast_Sort(distance)

            if (distance(1)<=cutoff_r) then ! diffusion condition
                x0=R(i,j)%x(1)
                y0=R(i,j)%x(2)
                z0=R(i,j)%x(3)
                dx(i,j)=dx(i,j-1)
                dy(i,j)=dy(i,j-1)
                dz(i,j)=dz(i,j-1)
            else if (distance(1)>cutoff_r) then ! i condition
                dx(i,j)=dx(i,j-1)
                dy(i,j)=dy(i,j-1)
                dz(i,j)=dz(i,j-1)
            end if
        else
            dx(i,j)=temp(1)
            dy(i,j)=temp(2)
            dz(i,j)=temp(3)
        end if
        !write(300,*) i,j,dx(i,j),dy(i,j),dz(i,j)
    end do
    d_x=differentiate(input_step,dx(i,:))
    d_y=differentiate(input_step,dy(i,:))
    d_z=differentiate(input_step,dz(i,:))
    !x:
    m=0
    do j=1,step-1
        if (d_x(j)*d_x(j+1)<0) then
            V_x(i)=V_x(i)+dx(i,j+1)
            m=m+1
            t_x(m)=j+1
        end if
    end do
    V_x(i)=V_x(i)/m
    do j=1,step-1
        if (t_x(j+1)==0) then
            cycle
        end if
        px(i)=(t_x(j+1)-t_x(j))*1.0+px(i)
    end do
    px(i)=2*px(i)/(m*1.0)
    px(i)=(1*2e13)/px(i)

    !y:
    m=0
    do j=1,step-1
        if (d_y(j)*d_y(j+1)<0) then
            V_y(i)=V_y(i)+dy(i,j+1)
            m=m+1
            t_y(m)=j+1
        end if
    end do   
    V_y(i)=V_y(i)/m
    do j=1,step-1
        if (t_y(j+1)==0) then
            cycle
        end if
        py(i)=(t_y(j+1)-t_y(j))*1.0+py(i)
    end do
    py(i)=2*py(i)/(m*1.0)
    py(i)=(1*2e13)/py(i)

    !z:
    m=0
    do j=1,step-1
        if (d_z(j)*d_z(j+1)<0) then
            V_z(i)=V_z(i)+dz(i,j+1)
            m=m+1
            t_z(m)=j+1
        end if
    end do   
    V_z(i)=V_z(i)/m
    do j=1,step-1
        if (t_z(j+1)==0) then
            cycle
        end if
        pz(i)=(t_z(j+1)-t_z(j))*1.0+pz(i)
    end do
    pz(i)=2*pz(i)/(m*1.0)
    pz(i)=(1*2e13)/pz(i)

    V(i)=sqrt(V_x(i)**2+V_y(i)**2+V_z(i)**2)
    write(300,*) i,V_x(i),V_y(i),V_z(i),V(i)
    write(400,*) i,px(i),py(i),pz(i)
    call cpu_time(time)
    if (mod(i,10)==0) then
        write(*,*) i,'is done. ',' cpu time: ',time
    end if
end do
stop
end program