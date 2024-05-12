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

@dataclass
class Node:
    children: list
    postion: Vec3
    depth: int
    
map_pos = {
    1: Vec3(0,0,0),
    2: Vec3(1,0,0),
    3: Vec3(0,1,0),
    4: Vec3(1,1,0),

    5: Vec3(0,0,0),
    6: Vec3(1,0,0),
    7: Vec3(0,1,0),
    8: Vec3(1,1,0)
}

num_calls = 0

def create_full_tree(depth, position = Vec3(0,0,0), counter = 1):
    global num_calls
    num_calls += 1
    print(num_calls)
    
    if counter-1 == depth:
        return Node([],position * (1/counter),counter)
    
    return Node([create_full_tree(depth, map_pos[i+1], counter+1) for i in range(8)],position * (1/counter), counter)

root = create_full_tree(3)



# pprint.pp(root)
