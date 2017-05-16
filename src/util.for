      function find_data_file(basename)
        ! Find a required data file by searching several locations.
        ! Returns the path if found; exits the program if not.
        ! Can be used directly in an open() statement, e.g.:
        !
        !   open(<unit>, file=find_data_file(<basename>), status='old')
        !
        implicit none
        character(len=*), intent(in) :: basename
        character(len=1000) :: find_data_file, data_home
        character(len=100)  :: suffix
        logical :: file_exists

        ! check current directory first
        find_data_file = basename
        inquire(file=find_data_file, exist=file_exists)
        if (file_exists) return

        ! path suffix for the remaining candidates
        suffix = '/osu-hydro/' // basename

        ! check user directory: XDG_DATA_HOME or ~/.local/share if unset
        call get_environment_variable('XDG_DATA_HOME', data_home)
        if (len_trim(data_home) == 0) then
          call get_environment_variable('HOME', data_home)
          data_home = trim(data_home) // '/.local/share'
        end if

        find_data_file = trim(data_home) // suffix
        inquire(file=find_data_file, exist=file_exists)
        if (file_exists) return

        ! check system directories
        find_data_file = '/usr/local/share' // suffix
        inquire(file=find_data_file, exist=file_exists)
        if (file_exists) return

        find_data_file = '/usr/share' // suffix
        inquire(file=find_data_file, exist=file_exists)
        if (file_exists) return

        print *, basename, ' not found'
        call exit(1)
      end function
