# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:58:48 2016

@author: Visscher
"""

# LET'S DRAW A PURPLE SPHERE !
from vapory import *
import NucleosomeMC as nuc
from helixmc.pose import HelixPose
import numpy as np
import FiberMC as FMC

#config.POVRAY_BINARY = "path/to/the/povray/binary"


def MakeDNAobject(dna, dyads):
    NRL = dyads[1]-dyads[0]
    coords = dna.coord

    #subtract COM of fiber here
    # coords-=np.average(coords)

    #for dinucleosome, subtrackt frame of 1st nuc / center of linker DNA

    objects = []
    z = np.array([[0,0,1],[0,0,1],[0,0,1],[0,0,1]])
    n=nuc.NucPose()
    n.from_file('1KX5.3DNA')
    P = nuc.of2coords(nucl.coords[n.d_index+n.fixed_i[7]], np.transpose(n.frames[n.d_index+n.fixed_i[7]]))
    n_coords = []
    j = 0
    for d in dyads:
        fixed_frame = np.transpose(dna.frames[d+n.fixed_i[7]-NRL%2-2])
        fixed_origin = dna.coord[d+n.fixed_i[7]-NRL%2-2]
        Q = nuc.of2coords(fixed_origin, fixed_frame)
        P = nuc.of2coords(n.coords[n.d_index+n.fixed_i[7]], \
                            np.transpose(n.frames[n.d_index+n.fixed_i[7]]))
        tf = nuc.get_transformation(P,Q)
        n_coords.append(nuc.apply_transf_coords(n.n_coords,tf))

    #create nucleosome
    h = 25
    rc = 35
    # for i in range(len(dyads)):
    #    objects.append(Sphere(n_coords[i][0,:],rc, Texture( Pigment( 'color', [.9,0,0] ))))
#    
    vec1 = np.array([0,-7,-3])
    vec2 = np.array([0,7,3])
    vec3 = np.array([0,1.75,0])
    vec6 = np.array([0,-1.75,0])
    vec5 = np.array([1.8,5.6,0])
    vec4 = np.array([1.8,-5.6,0])    
#    vec1 = np.array([0,-1.75,0])
#    vec2 = np.array([1.8,-5.6,0])
#    vec3 = np.array([1.8,-7,3])
#    vec6 = np.array([0,-1.75,0])
#    vec5 = np.array([-1.8,5.6,0])
#    vec4 = np.array([-1.8,7,-3])
    frames = dna.frames
    
    r = 4.5
    for i in range(len(coords)):
        color = [.9, .9, .9]
        objects.append(Sphere(np.dot(frames[i],vec1)+coords[i], r, Texture( Pigment( 'color', color ))))
        objects.append(Sphere(np.dot(frames[i],vec2)+coords[i], r, Texture( Pigment( 'color', color ))))
        objects.append(Sphere(np.dot(frames[i],vec3)+coords[i], r, Texture( Pigment( 'color', color ))))
        objects.append(Sphere(np.dot(frames[i],vec4)+coords[i], r, Texture( Pigment( 'color', color ))))
        objects.append(Sphere(np.dot(frames[i],vec5)+coords[i], r, Texture( Pigment( 'color', color ))))
        objects.append(Sphere(np.dot(frames[i],vec6)+coords[i], r, Texture( Pigment( 'color', color ))))
    return objects


def CreateScene(objects, title):
    # bert params
    # camloc = [500,500,750]
    # camdir = [0,0,750]
    cam_loc = [250,175, 0] # bovenaanzicht
    # cam_loc = [250, 175, 0]  # wat voor aanzicht?
    cam_dir = [0, 0, 0]
    light_loc = [600, 600, 600]
    camera = Camera( 'location', cam_loc, 'look_at', cam_dir )
    light = LightSource(light_loc, 'color', [1, 1, 1] )
    objects = [light]+objects
    objects.append(Background( "color", [1,1,1] ))

    scene = Scene( camera, objects= objects)
    scene.render(title, width=resolution_factor*600, height=resolution_factor*400)

def GetCam(dna, r):
    COM = np.mean(dna.coord, axis = 0)
    camdir = COM
    camloc = COM + np.array([r,r,0])
    return camloc, camdir

# initiate program

resolution_factor=2  # double

dna, dyads, nucl = FMC.create_fiber2(0, 2, 197, dna_file = 'simulations\\20180215_197_-10_unwrap.npz', nuc_file = '1KX5.3DNA', CL = False)

objects = MakeDNAobject(dna, dyads)


CreateScene(objects, "20180215_197_-10_unwrap.png")

#n_nuc = 2 
#NRL = 197
#n_bp = 147 + NRL*2
#dna, dyads, nucl = FMC.create_dna(n_bp, n_nuc, NRL)
#dna_file = '180208_197_2000_2.npz'
#dna = HelixPose(params = np.load(dna_file)['params'] , frame0 = np.load(dna_file)['frame0'])
#
##dna, dyads, nucl = FMC.create_fiber(200, 12, 170, dna_file = 'DinucleosomePoses\\1KX5NRL170.npz')
##dna = HelixPose(np.load('NRL1701KX5render.npz')['params'])
#
#objects = MakeDNAobject(dna, dyads)
#camloc = [500,500,750]
#camdir = [0,0,750]
#CreateScene(camloc, camdir, objects, 'picture.png')
#dyad = 111
#dnas = []
#for i in range(0,6):
#    temp = np.load('dna%d.npz' %i)
#    dnas.append( HelixPose(params = temp['params'], frame0 = temp['frame0']) )
#    
#    
#
#    
#def makepicture(i):
#    objects = MakeDNAobject(dnas[i], dyad)
#    camloc, camdir = GetCam(dnas[i], dyad, 200)
#    CreateScene(camloc, camdir, objects, 'snapshot%d' %i)
#    return i+1
#
#makepicture(5)







