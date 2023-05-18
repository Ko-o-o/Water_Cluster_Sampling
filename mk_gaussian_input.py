#usage: In wsl ubuntu, python3 mk_gaussian input.py 100

import sys
import random
import numpy as np
import quaternion
from multiprocessing import Pool

DISTANCE_CUTOFF = 1.5
def check_collision(water1, water2):
    for atom1 in water1:
        for atom2 in water2:
            distance = np.linalg.norm(atom1 - atom2)
            if distance < DISTANCE_CUTOFF:
                return True
    return False


def make_cluster():
    #H2Oが維持できていない構造をどれくらい取り込むか
    #r = random.choice([0,0.1,0.2])
    side = random.uniform((10*n)**(1/3), (10*n)**(1/3)+2)
    #n1 = int((1-r)*n) #H2Oの個数
    #n2 = int(r*n) #H,H,Oの個数

    #n1個のpackmolインプット用water.xyzファイル作成
    a_name = ['O','H','H']

    #n個の水を配置
    coors = []
    #1個目
    r1 =  0.96#random.triangular(0.83, 1.2, 0.96) #random.gauss(0.96,0.13)
    r2 =  0.96#random.triangular(0.83, 1.2, 0.96)
    theta = np.radians(109.5)#random.triangular(60, 180, 109.5)) #gauss(109.5,15))
    coor = np.array([
        [0,0,0],
        [r1,0,0],
        [r2*np.cos(theta),r2*np.sin(theta),0]
    ])
    coor_i = np.zeros((3,3))
    #回転
    vartheta_x = np.pi*random.random() #x軸周り回転角
    vartheta_y = np.pi*random.random() #y回転角
    vartheta_z = np.pi*random.random() #z回転角
    ux = np.quaternion(np.cos(vartheta_x), np.sin(vartheta_x), 0, 0) #軸
    uy = np.quaternion(np.cos(vartheta_y), 0, np.sin(vartheta_y), 0)
    uz = np.quaternion(np.cos(vartheta_z), 0, 0, np.sin(vartheta_z))
    for j in range(1,3):
        b = uz * uy * ux * np.quaternion(0,coor[j,0],coor[j,1],coor[j,2]) * ux.inverse() * uy.inverse() * uz.inverse()
        coor_i[j] = np.array([b.imag[0], b.imag[1], b.imag[2]])
    #並進
    R = np.array([random.uniform(0, side),random.uniform(0, side),random.uniform(0, side)])#np.array([random.uniform(side-1, side+1),random.uniform(side-1, side+1),random.uniform(side-1, side+1)])#np.random.rand(3)*random.uniform(side-1, side+1)
    coor_i += R[np.newaxis, :]
    coors.append(coor_i)

    #2個目以降
    for i in range(1,n):
        #ある確率で変な構造発生
        random1 = random.random()
        random2 = random.random()
        random3 = random.random()
        r1 = random.triangular(0.83, 1.2, 0.96) if random1<0.2 else 0.96 #?割りの確率で変な構造
        r2 = random.triangular(0.83, 1.2, 0.96) if random2<0.2 else 0.96
        theta = np.radians(random.triangular(60, 180, 109.5)) if random3<0.2 else 109.5
        coor = np.array([
            [0,0,0],
            [r1,0,0],
            [r2*np.cos(theta),r2*np.sin(theta),0]
        ])
        coor_i = np.zeros((3,3))
        #回転
        vartheta_x = np.pi*random.random() #x軸周り回転角
        vartheta_y = np.pi*random.random() #y回転角
        vartheta_z = np.pi*random.random() #z回転角
        ux = np.quaternion(np.cos(vartheta_x), np.sin(vartheta_x), 0, 0) #軸
        uy = np.quaternion(np.cos(vartheta_y), 0, np.sin(vartheta_y), 0)
        uz = np.quaternion(np.cos(vartheta_z), 0, 0, np.sin(vartheta_z))
        for j in range(1,3):
            b = uz * uy * ux * np.quaternion(0,coor[j,0],coor[j,1],coor[j,2]) * ux.inverse() * uy.inverse() * uz.inverse()
            coor_i[j] = np.array([b.imag[0], b.imag[1], b.imag[2]])
        collision = True
        while collision:
            print(f'i={i+1}, Collision True.')
            #並進
            R = np.array([random.uniform(0, side),random.uniform(0, side),random.uniform(0, side)])#np.random.rand(3)*random.uniform(side-1, side+1)
            coor_check_coli = coor_i + R[np.newaxis, :]
            for existing_water in coors:
                collision = check_collision(coor_check_coli,existing_water)
                if collision:break
        coors.append(coor_check_coli)
    return coors

n = int(sys.argv[1])

def mk_file(mn):
    Gtxt = f'# force AM1\n\nWater{n}\n\n0 1\n'
    cluster = make_cluster()
    for i in range(n):
        Gtxt += f' O        {cluster[i][0,0]:.6f}    {cluster[i][0,1]:.6f}    {cluster[i][0,2]:.6f}\n'
        Gtxt += f' H        {cluster[i][1,0]:.6f}    {cluster[i][1,1]:.6f}    {cluster[i][1,2]:.6f}\n'
        Gtxt += f' H        {cluster[i][2,0]:.6f}    {cluster[i][2,1]:.6f}    {cluster[i][2,2]:.6f}\n'
    Gtxt += '\n'
    with open(f'gjfs/asw_{mn+201}.gjf', 'w') as f:
        f.write(Gtxt)

with Pool(20) as p:
    p.map(mk_file, range(300))
