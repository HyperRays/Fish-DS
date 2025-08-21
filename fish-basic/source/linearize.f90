module linearize

    use comm_module
    use mhd_parameter_module, only: ns, nv
    use prof_module

    implicit none

contains

    subroutine index_to_3d(index, indices3D, nx, ny, nz)

        integer, intent(in) :: index
        integer, dimension(3), intent(out) :: indices3D
        type(comm_wld), intent(in) :: nx, ny, nz
        integer :: ix, iy, iz
        integer :: idx

        idx = index - 1

        ! Convert linear index to 3D indices using row-major ordering
        ! For array(nx%l, ny%l, nz%l), linear index = i + j*nx%l + k*nx%l*ny%l
        iz = idx / (nx%l * ny%l) + 1
        iy = mod(idx, nx%l * ny%l) / nx%l + 1
        ix = mod(idx, nx%l) + 1

        indices3D = [ix, iy, iz]
        
        
    end subroutine index_to_3d

    subroutine get_neighbor_index(index, neighbor, nx,ny,nz, indexout, error)
        integer, intent(in) :: index
        integer, dimension(3), intent(in) :: neighbor
        integer, intent(out) :: indexout
        type(comm_wld), intent(in) :: nx,ny,nz
        logical, intent(out) :: error
        integer :: indices3D(3)
        integer :: new_indices(3)
        integer :: i

        call index_to_3d(index, indices3D, nx, ny, nz)
        do i = 1, 3
            new_indices(i) = indices3D(i) + neighbor(i)
        end do

        ! Check bounds
        if (new_indices(1) < 1 .or. new_indices(1) > nx%l .or. &
            new_indices(2) < 1 .or. new_indices(2) > ny%l .or. &
            new_indices(3) < 1 .or. new_indices(3) > nz%l) then
            error = .true.
            indexout = -1
            return
        else
            error = .false.
        end if

        ! Convert new_indices back to linear index using row-major ordering
        ! linear index = (i-1) + (j-1)*nx%l + (k-1)*nx%l*ny%l
        indexout = (new_indices(1)-1) + (new_indices(2)-1)*nx%l + (new_indices(3)-1)*nx%l*ny%l + 1

    end subroutine get_neighbor_index

    subroutine is_out_of_bounds(indices, nx, ny, nz, out_of_bounds)
        integer, dimension(3), intent(in) :: indices
        type(comm_wld), intent(in) :: nx, ny, nz
        logical, intent(out) :: out_of_bounds

        out_of_bounds = (indices(1) < 1 .or. indices(1) > nx%l-1 .or. &
                         indices(2) < 1 .or. indices(2) > ny%l-1 .or. &
                         indices(3) < 1 .or. indices(3) > nz%l-1)
    end subroutine is_out_of_bounds

    subroutine scalar_array(array, index, nx, ny, nz, values)

        type(comm_wld), intent(in) :: nx, ny, nz
        integer, intent(in) :: index
        integer, dimension(3) :: indices
        real(4), dimension(ns,nx%l,ny%l,nz%l), intent(in) :: array
        real(4), dimension(ns), intent(out) :: values
        
        call index_to_3d(index, indices, nx, ny, nz)
        values = array(:, indices(1), indices(2), indices(3))

    end subroutine scalar_array

    subroutine vector_array(array, index, nx, ny, nz, values)

        type(comm_wld), intent(in) :: nx, ny, nz
        integer, intent(in) :: index
        integer, dimension(3) :: indices
        real(4), dimension(3,nv,nx%l,ny%l,nz%l), intent(in) :: array
        real(4), dimension(3, nv), intent(out) :: values
        
        call index_to_3d(index, indices, nx, ny, nz)
        values = array(:, :, indices(1), indices(2), indices(3))

    end subroutine vector_array

    subroutine linear_count(nx, ny, nz, count)

        type(comm_wld), intent(in) :: nx, ny, nz
        integer, intent(out) :: count

        count = nx%l * ny%l * nz%l
    
    end subroutine linear_count
end module linearize