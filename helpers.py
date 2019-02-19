import numpy

die_faces = "⚀⚁⚂⚃⚄⚅"

def roll_die():
    return die_faces[numpy.random.randint(len(die_faces))]

