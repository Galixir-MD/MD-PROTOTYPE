


def integrate(box, coords, velocities, forces, masses, time_step, lambda_T):
    ekinetic = 0
    for i in range(len(coords)):
        velocity2 = 0
        for m in range(3):
            # v = T * (v + F/m*t), 
            v_new = lambda_T * (velocities[i][m] + forces[i][m]*time_step/masses[i])
            coords[i][m] += (v_new * time_step)
            velocities[i][m] = v_new
            velocity2 += v_new**2
        ekinetic += 0.5 * masses[i] * velocity2
    
    return [ ekinetic, coords, velocities ]