module tridiagonal
        use :: precision
implicit none

interface tri_solve
        module procedure ::&
                tri_solve_vec, tri_solve_mat
end interface tri_solve
contains


function tri_det (A) result (det)
! Calculate determinant of tridiagonal matirx
! input:
! A - 0 b1 c1
!     a2 b2 c2
!       ...
! output: x - real number
        real(prec), dimension(1:, 1:), intent(in) :: A
        real(prec) :: det

        real(prec), dimension(-2:0) :: f
        integer :: n, i

        n = size(A(1, :))

        !f(-2) = 0
        f(-1) = 1
        f(0) = A(2, 1)

        do i = 2, n
                f(-2) = f(-1)
                f(-1) = f(0)
                f(0) = A(2, i)*f(-1) - A(1, i)*A(3, i-1)*f(-2)
        enddo

        det = f(0)

end function tri_det


function tri_solve_vec (A_in, b) result (x)
! Find roots of tridiagonal SLE
! input:
! A - 0 b1 c1
!     a2 b2 c2
!       ...
!     an bn 0
! b - y1
!     y2
!     ...
!     yn
! output: x - roots of SLE
! x - x1
!     x2
!     ...
!     xn
        real(prec), dimension(1:, 1:), intent(in) :: A_in
        real(prec), dimension(1:), intent(in) :: b

        real(prec), dimension(3, 1:size(b)) :: A
        real(prec), dimension(1:size(b)) :: x

        real(prec), dimension(1:size(b)) :: alpha, beta
        real(prec) :: det

        integer :: n, i

        n = size(b)
        if (size(A, 2) .ne. 3) then
                A = transpose(A_in)
        else
                A = A_in
        endif

        det = tri_det(A)

        if (abs(det) < ZERO) then
                write (*, '("Error: the matrix is singular, det = ", e8.2)') det
                stop 
        endif


        alpha(1) = -A(3, 1) / A(2, 1)
        beta(1) = b(1) / A(2, 1)

        do i = 2, n
                alpha(i) = -A(3, i) / (A(1, i) * alpha(i-1) + A(2, i))
                beta(i) = (b(i) - A(1, i) * beta(i-1)) / (A(1, i) * alpha(i-1) + A(2, i))
        enddo


        x(n) = beta(n)

        do i = n-1, 1, -1
                x(i) = alpha(i) * x(i+1) + beta(i)
        enddo


end function tri_solve_vec

function tri_solve_mat (A_in, b) result (x)
! Find roots of tridiagonal SLE
! input:
! A - 0 b1 c1
!     a2 b2 c2
!       ...
!     an bn 0
! b - y1 z1 ...
!     y2 z2 ...
!     ...
!     yn zn ...
! output: x - roots of SLE
! x - x1
!     x2
!     ...
!     xn
        real(prec), dimension(1:, 1:), intent(in) :: A_in
        real(prec), dimension(1:, 1:), intent(in) :: b
        real(prec), dimension(1:size(b, 1), 1:size(b, 2)) :: x

        integer :: k, i

        k = size(b, 2)

        do i = 1, k
                x(:, i) = tri_solve_vec(A_in, b(:, i))
        enddo
end function tri_solve_mat


end module tridiagonal
