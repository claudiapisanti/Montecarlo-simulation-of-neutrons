# FUNCTIONS
#%% FUNZIONI da mettere in un altro file
def source(face, x_rect,y_rect,z_rect):
    x_source = random() * x_rect
    y_source = random() * y_rect
    z_source = random() * z_rect
    pos = [0,0,0]

    if (face == 1): pos = np.array([x_source, y_source, z_rect])
    if (face == 2): pos = np.array([x_rect, y_source, z_source])
    if (face == 3): pos = np.array([x_source, 0.0, z_source])
    if (face == 4): pos = np.array([0.0, y_source, z_source])
    if (face == 5): pos = np.array([x_source, y_rect, z_source])
    if (face == 6): pos = np.array([x_source, y_source, 0.0])

    return pos