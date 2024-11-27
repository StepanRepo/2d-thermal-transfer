module multidim
        use :: precision
        use :: tridiagonal
        use :: nle_solve
        implicit none

        real(prec) :: dx, dy, t, temp_b
        real(prec) :: k1, s1, k2, s2
        real(prec), dimension(:), allocatable :: x, y
        real(prec), dimension(:, :), allocatable :: u
        integer :: n, k

        private

        public :: init, step, save_state, t

contains



subroutine init(xx, yy, zz, tb, kk1, ss1, kk2, ss2)
        real(prec), dimension(0:) :: xx, yy
        real(prec), dimension(0:, 0:) :: zz
        real(prec) :: tb, kk1, kk2, ss1, ss2

        t = 0.0_prec

        dx = xx(1) - xx(0)
        dy = yy(1) - yy(0)

        allocate(x(0:size(xx)-1), y(0:size(yy)-1))
        allocate(u(0:(size(xx)-1), 0:(size(yy)-1)))

        x = xx
        y = yy
        u = zz

        temp_b = tb
        k1 = kk1
        s1 = ss1
        k2 = kk2
        s2 = ss2

        n = size(u, 1)
        k = size(u, 2)
end subroutine


subroutine step (courant)
        real(prec) :: dt, courant

        real(prec), dimension(0:n-1, 0:k-1) :: unew, temp
        integer :: iter

        ! find maximal D
        dt = maxval(k1*u**s1)
        dt = max(maxval(k2*u**s2), dt)

        dt = courant* (dx*dx + dy*dy)/4.0_prec / dt

        t = t+dt

        !-----------
        ! make a step along x-axis
        !-----------
        temp = 0.0_prec
        unew = u
        iter = 0

        ! iterate over solution until it converges
        do while(sum(abs(unew - temp))/n/k > 10.0_prec*ZERO)

                temp = unew

                ! find new solution
                unew = fx(temp)


                ! throw an exception if th solution 
                ! doesn't converge 
                iter = iter + 1

                if (iter > 100) then
                        write (*, "('Cant solve x after ', i6, ' iterations. t = ', f16.8, ', error = ', es16.8)")&
                                iter, t, sum(abs(unew - temp))/n/k
                        stop
                endif

        enddo

        u = unew

        !-----------
        ! make a step along y-axis
        !-----------
        temp = 0.0_prec
        unew = u
        iter = 0

        ! iterate over solution until it converges
        do while(sum(abs(unew - temp))/n/k > ZERO)

                temp = unew

                ! find new solution
                unew = fy(temp)



                ! throw an exception if th solution 
                ! doesn't converge 
                iter = iter + 1

                if (iter > 100) then
                        write (*, "('Cant solve y after ', i6, ' iterations. t = ', f16.8, ', error = ', es16.8)")&
                                iter, t, sum(abs(unew - temp))/n/k
                        stop
                endif

        enddo


        u = unew


contains 

function fx(unew) result (res)
        real(prec), dimension(0:n-1, 0:k-1) :: unew, res

        real(prec), dimension(0:n-1) :: Bx
        real(prec), dimension(0:n-1, 0:2) :: Ax
        real(prec), dimension(0:n-1, 0:k-1) :: D

        real(prec) :: c
        integer :: i

        c = dt/dx/dx/2.0_prec

        Ax(0, 1) = 1.0_prec
        Ax(0, 2) = 0.0_prec
        Ax(n-1, 0) = 0.0_prec
        Ax(n-1, 1) = 1.0_prec


        D = k1*(unew**s1)

        do i = 1, k-2
                Ax(1:n-2, 0) = -c * (D(0:n-3, i) + D(1:n-2, i))/2.0_prec
                Ax(1:n-2, 2) = -c * (D(1:n-2, i) + D(2:n-1, i))/2.0_prec

                Ax(1:n-2, 1) = c*(D(1:n-2, i) + &
                        .5_prec*D(0:n-3, i) + .5_prec*D(2:n-1, i)) + 1.0_prec

                Bx = (u(:, i)**s2 * (u(:, i+1) - u(:, i)) -&
                        u(:, i-1)**s2 * ((u(:, i) - u(:, i-1)))) * k2*dt/dy/dy/2.0_prec

                Bx = Bx + u(:, i)

                res(:, i) = tri_solve(Ax, Bx)
        enddo
        

        ! boundary conditions
        res(0, :)   = temp_b
        res(n-1, :) = temp_b
        res(:, 0)   = res(:, 1)
        res(:, k-1) = res(:, k-2)
end function fx

function fy(unew) result (res)
        real(prec), dimension(0:n-1, 0:k-1) :: unew, res

        real(prec), dimension(0:k-1) :: By
        real(prec), dimension(0:k-1, 0:2) :: Ay
        real(prec), dimension(0:n-1, 0:k-1) :: D

        real(prec) :: c
        integer :: i

        c = dt/dy/dy/2.0_prec


        ! this is from boundary conditions
        Ay(0, 1) = 1.0_prec
        Ay(0, 2) = -1.0_prec
        Ay(k-1, 0) = -1.0_prec
        Ay(k-1, 1) = 1.0_prec


        D = k2*(unew**s2)

        do i = 1, n-2
                Ay(1:k-2, 0) = -c * (D(i, 0:k-3) + D(i, 1:k-2))/2.0_prec
                Ay(1:k-2, 2) = -c * (D(i, 1:k-2) + D(i, 2:k-1))/2.0_prec

                Ay(1:k-2, 1) = c*(D(i, 1:k-2) + &
                        .5_prec*D(i, 0:k-3) + .5_prec*D(i, 2:k-1)) + 1.0_prec

                By = (u(i, :)**s1 * (u(i+1, :) - u(i, :)) -&
                        u(i-1, :)**s1 * ((u(i, :) - u(i-1, :)))) * k1*dt/dx/dx/2.0_prec

                By = By + u(i, :)


                By(0) = 0.0_prec
                By(k-1) = 0.0_prec


                res(i, :) = tri_solve(Ay, By)
        enddo

        res(0, :)   = temp_b
        res(n-1, :) = temp_b
        !res(:, 0)   = temp_b 
        !res(:, k-1) = temp_b 

        res(:, 0)   = res(:, 1)
        res(:, k-1) = res(:, k-2)

end function fy




end subroutine step

subroutine save_state(filename)
        character (len = *) :: filename

        integer :: fu, io, i
        logical :: if_exists

        fu = 1

        open (action = 'write', file = filename, iostat = io, newunit = fu)

        ! check whether file exists
        inquire (file = filename, exist = if_exists)

        if (.not. if_exists) then
                write (*, '("Cannot open file ", a)') filename
                close (fu)

                STOP 1
        end if

        write(fu, *) t
        write(fu, *) x
        write(fu, *) y

        do i = 0, n-1
                write(fu, *) u(i, :)
        enddo

        if (io /= 0) then
                write (*, '("Error occured while writing to the file ", a)')&
                        filename
                close (fu)
        end if
        
        close (fu)

end subroutine save_state

end module multidim
