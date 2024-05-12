
    module octree

        type, public :: Node
            type(Node), dimension(8), allocatable:: children(:)
            integer, dimension(3) :: position
        end type Node

        contains
                subroutine hello(this)
                    type(Node), intent(in) :: this
                    write (*,*) this%position
                end subroutine hello
    end module octree