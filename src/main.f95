program rieman
        use :: precision
        use :: multidim
        use :: tridiagonal
        use :: config
implicit none

real(prec), dimension(:), allocatable :: x, y
real(prec), dimension(:, :), allocatable :: vec
real(prec) :: dx, dy
integer :: i

character(len = 10) :: filename

call configure()

allocate(x(0:xn-1), y(0:yn-1), vec(0:xn-1, 0:yn-1))

dx = (x1-x0) / real(xn-1, prec)
dy = (y1-y0) / real(yn-1, prec)


do i = 0, xn-1
        x(i) = real(i, prec)*dx + x0
enddo

do i = 0, yn-1
        y(i) = real(i, prec)*dy + y0
enddo

vec = temp_0
vec(0, :)    = temp_b
vec(xn-1, :) = temp_b
vec(:, 0)    = temp_b
vec(:, yn-1) = temp_b

call init(x, y, vec, temp_b, k1, s1, k2, s2)
call save_state(trim(output_dir)//"/000.dat")

i = 0

do while (t < t1)

        call step(C) 

        write(*, "('Culculating... t = ', es16.8, a)", advance="no")  t, achar(13)

        i = i + 1

        if (mod(i, every) == 0) then
                write(filename, "(i3.3, '.dat')") i/every
                call save_state(trim(output_dir)//"/"//trim(filename))
        endif
enddo

write(*, *) ""

deallocate(x, y, vec)












contains

subroutine configure()
        implicit none

        character(len = 100):: arg, config_file = ""
        integer :: idx

        idx = 1
        do
                if (idx > command_argument_count()) exit
                call get_command_argument(idx, arg)

                select case (arg)

                case ("-h", "--help")
                        write(*, "(a)") "Command-line options:"
                        write(*, "(a)") "  -h, --help             print this help message"
                        write(*, "(a)") "  -c, --cfg {file_name}  set configuration file"

                        STOP 0

                case ("-c", "--cfg")
                        idx = idx+1
                        call get_command_argument(idx, arg)

                        config_file = arg

                case default
                        write(*, "(2a)") "Unknown command line option: ", trim(arg)
                end select

                idx = idx+1
        end do

        if (trim(config_file) /= "") then
                call read_config(trim(config_file)) ! read custom config
        else
                call read_config("default.cfg") ! read default config
        endif

end subroutine configure



end program 
