# pylint: disable=redefined-outer-name
## A python implementation of the octree
# 1 indexed to ease reimplementation in fortran
from dataclasses import dataclass
import pprint
from typing import Any

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

    def __iter__(self):
        for v in (self.x,self.y,self.z):
            yield v

@dataclass
class Node:
    children: list
    data: Any
    
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

#                           |Debug only          |            
def create_full_tree(depth, position = Vec3(0,0,0), counter = 1):
    if counter == depth:
        return Node([None for _ in range(8)],(position, counter))
    
    return Node([create_full_tree(depth, position + (map_pos[i+1] << (counter-1)), counter+1) for i in range(8)],(position, counter))

## naive code
# def bfs_find_node(root: Node, target: Vec3, depth: int):
#     if (root.position == target) and (root.depth == depth):
#         return root

#     for node in root.children:
#         value = bfs_find_node(node, target, depth)
#         if value != None:
#             return value


# use postion as map
def find_position(root: Node, target: Vec3, depth: int, last_found = False):
    for i in range(depth-1):
        
        idx = map_child[(target >> i) & Vec3(1,1,1)]
        node = root.children[idx-1]

        if node is None and not last_found:
            return None
        if last_found:
            return root
        root = node
        
    return root

def inv(v):
    if v == 0:
        return 0 
    return 1/v

def addr_to_abs(position: Vec3):
    position = position << 1
    return Vec3(inv(position.x),inv(position.y),inv(position.z))

def abs_to_addr(position: Vec3):
    return to_int(Vec3(inv(position.x),inv(position.y),inv(position.z))) >> 1
 
def to_int(vec):
    return Vec3(
        int(vec.x),
        int(vec.y),
        int(vec.z)
    )


max_height = 2

root = create_full_tree(max_height)
# pprint.pp(root)

search = find_position(root, Vec3(0,0,0), 2)
if search != None:
    print("Search Result:",search.data)
else:
    print("Found Nothing")

v = search.data[0] # type: ignore
abs_addr = addr_to_abs(v) + Vec3(0,0,1/2**(2-1))
rel_addr = to_int(abs_to_addr(abs_addr))

print("Absolute Address:",abs_addr)
print("Relative Address", rel_addr)
print("New Node",find_position(root, rel_addr, max_height).data) # type: ignore
