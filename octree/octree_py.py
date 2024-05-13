# pylint: disable=redefined-outer-name
## A python implementation of the octree
# 1 indexed to ease reimplementation in fortran
from dataclasses import dataclass
import pprint

@dataclass(slots=True, frozen=True)
class Vec3:
    x: int | float
    y: int | float
    z: int | float

    def __add__(self, other):
        return Vec3(self.x + other.x, self.y + other.y, self.z+other.z)

    def __neg__(self):
        return Vec3(-self.x, -self.y, -self.z)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Vec3(self.x * other,
            self.y * other,
            self.z * other)

    # useful for better accuracy method
    def __lshift__(self, other):
        return Vec3(
            self.x << other,
            self.y << other,
            self.z << other
        )

    def __rshift__(self, other):
        return Vec3(
            self.x >> other,
            self.y >> other,
            self.z >> other
        )
    
    def __eq__(self, other):
        return self.x == other.x and \
                self.y == other.y and \
                self.z == other.z
    
    def __and__(self, other):
        return Vec3(self.x & other.x,
                self.y & other.y,
                self.z & other.z)

@dataclass
class Node:
    children: list
    data: any
    
map_pos = {
    1: Vec3(0,0,0),
    2: Vec3(1,0,0),
    3: Vec3(0,1,0),
    4: Vec3(1,1,0),

    5: Vec3(0,0,1),
    6: Vec3(1,0,1),
    7: Vec3(0,1,1),
    8: Vec3(1,1,1)
}

map_child = {
    Vec3(0,0,0): 1,
    Vec3(1,0,0): 2,
    Vec3(0,1,0): 3,
    Vec3(1,1,0): 4,

    Vec3(0,0,1): 5,
    Vec3(1,0,1): 6,
    Vec3(0,1,1): 7,
    Vec3(1,1,1): 8
}

def create_full_tree(depth, position = Vec3(0,0,0), counter = 1):
    if counter == depth:
        return Node([None for _ in range(8)],(position, counter))
    
    return Node([create_full_tree(depth, position + map_pos[i+1] << (counter-1), counter+1) for i in range(8)],(position, counter))

## naive code
# def bfs_find_node(root: Node, target: Vec3, depth: int):
#     if (root.position == target) and (root.depth == depth):
#         return root

#     for node in root.children:
#         value = bfs_find_node(node, target, depth)
#         if value != None:
#             return value

# use postion as map
def find_position(root: Node, target: Vec3, depth: int):
    for i in range(depth-1):
        idx = map_child[(target >> i) & Vec3(1,1,1)]
        node = root.children[idx-1]
        if node == None:
            return None
        else:
            root = node
        
    return root

def addr_to_abs(position: Vec3, max_depth):
    return sum((((position >> i) & Vec3(1,1,1))*(1/2**(i+1)) for i in range(max_depth-1)), start=Vec3(0,0,0))


root = create_full_tree(2)
pprint.pp(root)

print(addr_to_abs(Vec3(6,2,4), 5))

search = find_position(root, Vec3(2,2,2), 3)
if search != None:
    print(search.data)
else:
    print("Found Nothing")