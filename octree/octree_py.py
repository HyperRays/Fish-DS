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
    position: Vec3
    depth: int
    
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
        return Node([],position,counter)
    
    return Node([create_full_tree(depth, position + map_pos[i+1] << (counter-1), counter+1) for i in range(8)], position, counter)

def bfs_find_node(root: Node, target: Vec3, depth: int):
    # naive code
    if (root.position == target) and (root.depth == depth):
        return root

    for node in root.children:
        value = bfs_find_node(node, target, depth)
        if value != None:
            return value

# use postion as map
def find_position(root: Node, target: Vec3, depth: int, count = 0):
    idx = map_child[target & Vec3(1,1,1) << count]
    if root.depth == depth:
        return root
    
    return find_position(root.children[idx-1], target, depth, count+1)


root = create_full_tree(3)
pprint.pp(root)

search = find_position(root, Vec3(1,0,0), 2)
if search != None:
    print(search.position, search.depth)
else:
    print("Found Nothing")

# pprint.pp(root)
