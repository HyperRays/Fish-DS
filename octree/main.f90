    
    program ocprog
        use octree
        type(Node) :: value
        integer, dimension(3) :: pos
        pos = (/1,1,1/)
        value%position = pos
        
        call hello(value)

    end program ocprog