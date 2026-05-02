

module STORAGE 
    use cudafor
    real(8), parameter :: par_pi_8 = acos(-1.0_8)         
    real(8), parameter :: par_pi_4 = acos(-1.0_4)      
	real(8), parameter :: par_sqrtpi = sqrt(par_pi_8)
    real(8), parameter :: cpi4 = 12.56637061435917295384_8
    real(8), parameter :: ggg = (5.0_8/3.0_8)
    real(8), parameter :: par_R_character = 1.0_8


    integer(4), parameter  :: host_num_param = 8     ! Число параметров в ячейке


    real(8) :: host_time_all                       ! Текущее время расчёта
    INTEGER(4) :: host_N_cell                       ! Число ячеек в сетке
    INTEGER(4) :: host_N_gran                       ! Число граней в сетке

	real(8), allocatable :: host_Cell_par(:, :)     ! (8, :)

	real(8), allocatable :: host_Cell_center(:, :)      ! (3, :)
	real(8), allocatable :: host_Cell_Volume(:)         ! (:)
    integer(4), allocatable :: host_Cell_gran(:,:)            ! (6, :) Набор из 6 граней для каждой ячейки (если номер = 0, то грани нет в этом направлении)
	real(8), allocatable :: host_Cell_dist(:)                 ! Минимальное расстояние до грани в каждой ячейки 


	real(8), allocatable :: host_Gran_normal(:,:)       ! (3, :) Нормаль грани   
    real(8), allocatable :: host_Gran_square(:)         ! (:) Площадь грани
	real(8), allocatable :: host_Gran_center(:,:)			! (3, :)

	integer(4), allocatable :: host_Gran_neighbour(:,:) ! Соседи каждой грани (2,:) имеют по 2 соседа, нормаль ведёт от первого ко второму
	integer(4), allocatable :: host_Gran_neighbour_TVD(:,:) ! TVD-Соседи каждой грани (2,:) имеют по 2 соседа

	integer(4), allocatable :: host_Gran_type(:)                 ! Показывает тип грани

	real(8), allocatable :: host_Gran_POTOK(:, :)       ! (10, :) поток грани    последний - дивергенция магнитного поля для очистки

    real(8), allocatable, public :: host_cooling_T(:)
    real(8), allocatable, public :: host_cooling_Lambda(:)
    real(8), allocatable, public :: host_cooling_lnT(:)
    real(8), allocatable, public :: host_cooling_lnLambda(:)

    real(8), allocatable, public :: host_heating_T(:)
    real(8), allocatable, public :: host_heating_Lambda(:)
    real(8), allocatable, public :: host_heating_lnT(:)
    real(8), allocatable, public :: host_heating_lnLambda(:)

    contains

    subroutine Set_Storage()
        real(8) :: vv
        integer(4) :: i
        logical :: file_exists
        integer :: ierr
        character(len=256) :: filename

        filename = 'FCMHD_1.bin'

        ierr = 1
        inquire(file=filename, exist=file_exists)
        if (file_exists) then
            open(1, file = filename, FORM = 'BINARY', ACTION = "READ", iostat=ierr)
        else
            inquire(file='../'//filename, exist=file_exists)
            if (file_exists) then
                open(1, file = '../'//filename, FORM = 'BINARY', ACTION = "READ", iostat=ierr)
            end if
        end if

        if (ierr /= 0) then
            write(*,*) 'Error: cannot open FCMHD_1.bin (current dir, CUDA_FORT/ or ../)'
            stop
        end if


        read(1) host_time_all
        read(1) host_N_cell
        read(1) host_N_gran

        print*, "host_time_all, host_N_cell, host_N_gran = ", host_time_all, host_N_cell, host_N_gran

        allocate(host_Cell_par(host_num_param, host_N_cell))
        allocate(host_Cell_center(3, host_N_cell))
        allocate(host_Cell_Volume(host_N_cell))
        allocate(host_Cell_gran(6, host_N_cell))


        allocate(host_Gran_normal(3, host_N_gran))
        allocate(host_Gran_square(host_N_gran))
        allocate(host_Gran_center(3, host_N_gran))
        allocate(host_Gran_neighbour(2, host_N_gran))
        allocate(host_Gran_neighbour_TVD(2, host_N_gran))
        allocate(host_Gran_type(host_N_gran))
        allocate(host_Gran_POTOK(host_num_param + 1, host_N_gran))

        read(1) host_Cell_par
        !print*, "S1 = ", host_Cell_par(1, 1), host_Cell_par(2, 1), host_Cell_par(3, 1), host_Cell_par(4, 1), host_Cell_par(5, 1), host_Cell_par(6, 1), host_Cell_par(7, 1), host_Cell_par(8, 1)
        read(1) host_Cell_center
        read(1) host_Cell_Volume
        read(1) host_Cell_gran

        read(1) host_Gran_normal
        read(1) host_Gran_square
        read(1) host_Gran_center
        read(1) host_Gran_neighbour
        read(1) host_Gran_neighbour_TVD
        !print*, "S2 = ", host_Gran_neighbour_TVD(1, 2), host_Gran_neighbour_TVD(2, 2)
        read(1) vv
        !print*, "VV = ", vv
        read(1) host_Gran_type
        read(1) vv
        !print*, "VV = ", vv
        read(1) host_Gran_POTOK

        read(1) vv
        print*, "VV (123)  =  ", vv

        close(1)

        do i = 1, size(host_Cell_par(1, :)) 
            if(host_Cell_par(1, i) < 0.000000001) then
                print*, "rho < 0  =  ", host_Cell_par(1, i)
                STOP
            end if
        end do

        do i = 1, size(host_Gran_neighbour_TVD(1, :)) 
            if(host_Gran_neighbour_TVD(1, i) > host_N_cell) then
                print*, "ERROR host_Gran_neighbour_TVD"
                STOP
            end if

            if(host_Gran_neighbour_TVD(1, i) < 0) then
                print*, "ERROR host_Gran_neighbour_TVD", host_Gran_neighbour_TVD(1, i)
                STOP
            end if

            if(host_Gran_neighbour_TVD(2, i) < 0) then
                print*, "ERROR host_Gran_neighbour_TVD", host_Gran_neighbour_TVD(2, i)
                STOP
            end if
        end do

        call flush(6)


        call read_cooling_function('combined_cooling_function.txt', &
            host_cooling_T, host_cooling_Lambda, host_cooling_lnT, host_cooling_lnLambda)
        call read_cooling_function('combined_heating_function.txt', &
            host_heating_T, host_heating_Lambda, host_heating_lnT, host_heating_lnLambda)

    end subroutine Set_Storage


    subroutine Save_Storage()
        integer :: unit, ierr
        real(8) :: cf

        !3.0  0 градусов
        !3.1  10 градусов
        !3.2  30 градусов


        ! Открываем файл для записи в бинарном формате
        open(newunit=unit, file="FCMHD_1.3_out.bin", form='unformatted', access='stream', &
            action='write', status='replace', iostat=ierr)
        
        if (ierr /= 0) then
            print *, "Error opening file for writing: ", "FCMHD_3.0_out.bin"
            return
        endif
        
        ! Записываем данные в ТОМ ЖЕ порядке, что и при чтении
        write(unit) host_time_all
        write(unit) host_Cell_par 

        cf = 321.0_8
        write(unit) cf

        close(unit)
    end subroutine Save_Storage

    subroutine Fill_data()
        integer :: unit, ierr
        real(8) :: cf

        ! Открываем файл для записи в бинарном формате
        open(3, file = "FCMHD_1.2_out.bin", FORM = 'BINARY', ACTION = "READ")
        
        ! Записываем данные в ТОМ ЖЕ порядке, что и при чтении
        read(3) host_time_all
        read(3) host_Cell_par 

        read(3) cf
        print*, "Proverka (321)  ", cf

        close(3)
        call flush(6)
    end subroutine Fill_data

    !---------------------------------------------------------------------
    ! Чтение данных из файла, заполнение массивов T, Lambda и их логарифмов
    !---------------------------------------------------------------------
    subroutine read_cooling_function(filename, T, Lambda, lnT, lnLambda)
        character(len=*), intent(in)                          :: filename
        real(8), allocatable, intent(out)                :: T(:), Lambda(:), lnT(:), lnLambda(:)

        integer               :: unit, i, n, ierr
        character(len=256)    :: line
        real(8)          :: tval, lval
        logical               :: sorted

        ! Открываем файл
        open (newunit=unit, file=filename, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            write (*,*) 'Error: cannot open file ', trim(filename)
            stop
        end if

        ! Пропускаем строку заголовка
        read (unit, '(A)', iostat=ierr) line

        ! Временные динамические массивы
        allocate (T(0), Lambda(0))

        do
            read (unit, '(A)', iostat=ierr) line
            if (ierr < 0) exit            ! конец файла
            if (ierr > 0) exit            ! ошибка чтения
            if (len_trim(line) == 0) cycle
            read (line, *, iostat=ierr) tval, lval
            if (ierr /= 0) then
                write (*,*) 'Warning: skipped invalid line: ', trim(line)
                cycle
            end if
            ! Расширяем массивы на 1
            T = [T, tval]
            Lambda = [Lambda, lval]
        end do
        close (unit)

        n = size(T)
        if (n == 0) then
            write (*,*) 'Error: no data read from file.'
            stop
        end if

        ! Вычисляем логарифмы
        allocate (lnT(n), lnLambda(n))
        do i = 1, n
            lnT(i)      = log(T(i))
            lnLambda(i) = log(Lambda(i))
        end do

        ! Проверка монотонности по lnT
        sorted = .true.
        do i = 2, n
            if (lnT(i) <= lnT(i-1)) then
                sorted = .false.
                exit
            end if
        end do
        if (.not. sorted) then
        write (*,*) 'Warning: T values are not in increasing order.'
        end if

        ! Создаём проверочный файл
        call write_check_file(filename, T, Lambda, lnT, lnLambda)
    end subroutine read_cooling_function

    !---------------------------------------------------------------------
    ! Запись проверочного файла с интерполяцией на равномерной лог-сетке
    !---------------------------------------------------------------------
    subroutine write_check_file(filename, T, Lambda, lnT, lnLambda)
        character(len=*), intent(in)   :: filename
        real(8), intent(in)       :: T(:), Lambda(:), lnT(:), lnLambda(:)

        integer, parameter             :: nsteps = 10000
        real(8)                   :: T_min, T_max, logT_min, logT_max, log_step
        real(8)                   :: Tcur, Lcur
        integer                        :: i, unit, ierr

        T_min = T(1)
        T_max = T(size(T))
        logT_min = log(T_min)
        logT_max = log(T_max)
        log_step = (logT_max - logT_min) / real(nsteps-1, 8)

        open (newunit=unit, file=trim(filename)//'_check.txt', &
            status='replace', action='write', iostat=ierr)
        if (ierr /= 0) then
            write (*,*) 'Warning: could not create check file.'
            return
        end if

        write (unit, '(A)') '# T [K]    Lambda_interpolated [erg cm^3/s]'
        do i = 0, nsteps-1
            Tcur = exp(logT_min + i * log_step)
            Lcur = interpolate_cooling(Tcur, T, Lambda, lnT, lnLambda)
            write (unit, '(ES18.9,2X,ES18.9)') Tcur, Lcur
        end do
        close (unit)

        write (*,*) 'Check file written: ', trim(filename)//'_check.txt'
    end subroutine write_check_file


    !---------------------------------------------------------------------
    ! Интерполяция Lambda/n_H^2 по T (линейная в log-log)
    ! Для использования на GPU – pure функция (может вызываться из ядра)
    !---------------------------------------------------------------------
    !@cuf attributes(host, device) & 
    pure function interpolate_cooling(T_in, T_arr, Lambda_arr, lnT_arr, lnLambda_arr) result(lambda_out)
        real(8), intent(in) :: T_in
        real(8), intent(in) :: T_arr(:), Lambda_arr(:), lnT_arr(:), lnLambda_arr(:)
        real(8)             :: lambda_out

        real(8) :: logT, lnL
        integer      :: n, idx

        n = size(T_arr)
        ! Если температура вне диапазона данных – возвращаем 0
        if (T_in <= T_arr(1) .or. T_in >= T_arr(n)) then
            lambda_out = 0.0_8
            return
        end if

        logT = log(T_in)

        ! Поиск индекса первого lnT >= logT (аналог std::lower_bound)
        idx = 2
        do while (idx <= n .and. lnT_arr(idx) < logT)
            idx = idx + 1
        end do

        ! Точное совпадение с узлом (с учётом погрешности)
        if (idx <= n .and. abs(lnT_arr(idx) - logT) < 1.0e-12_8) then
            lambda_out = Lambda_arr(idx)
            return
        end if

        ! Линейная интерполяция в log-log пространстве
        lnL = lnLambda_arr(idx-1) + &
            (logT - lnT_arr(idx-1)) * (lnLambda_arr(idx) - lnLambda_arr(idx-1)) &
            / (lnT_arr(idx) - lnT_arr(idx-1))

        lambda_out = exp(lnL)
    end function interpolate_cooling


end module STORAGE