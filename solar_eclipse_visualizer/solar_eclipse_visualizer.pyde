import math
time = 358.88687782909255 #Time of New Moon of eclipse to be mapped

showGrid = True #Show longitude latitude grid lines (recommended)
longitudeCorrect = True #Tries to center the eclipse in the image

def setup():
    size(1500, 750) #Resolution of the image generated, 2nd number *must* be 1/2 the first.
    
R = 673320.0 #Radius of the sun in km
r = 5512.0 #Radius of the planet in km
rm = 2260.113 #Radius of the moon in km
D = 126096000.0 #Distance to the sun in km
d = 294000.241 #Distance to the moon in km
S = 298.0 #Planetary orbital period in days
moonPeriod = 22.4435146444 #Lunar sidereal orbital period in days
precession = radians(15.42840432988) #Lunar nodal precession rate in radians per year
tilt = math.radians(24.7) #Planetary axial tilt in radians
lunarInclination = radians(4.92) #Inclination of the Moon's orbit in radians
moon_offset = -78 #Copy from eclipse predictor program
node_offset = -78 #Copy from eclipse predictor program
earth_hours_in_day = 24 #Length of the synodic day in Earth hours

###DO NOT TOUCH FROM HERE BELOW###
def ecliptic_equat(coords):
    ascension = atan2(cos(tilt) * sin(coords[0]) - sin(tilt) * tan(coords[1]), cos(coords[0]))
    declension = asin(cos(tilt) * sin(coords[1]) + sin(tilt) * cos(coords[1]) * sin(coords[0]))
    if (ascension < 0):
        ascension += 2 * PI
    return [ascension, declension, coords[2]]

def calculate_sun_equat(time):
    sun_ecliptic = [((time % S) / S * 2 * PI) % (2 * PI), 0, D]
    return ecliptic_equat(sun_ecliptic)

def calculate_moon_equat(time):
    node = (((-(time - node_offset) % precessionPeriod) / precessionPeriod + (S + node_offset) / S) * 2 * PI) % (2 * PI)
    moon_long = ((((time - moon_offset) % moonPeriod) / moonPeriod + (S + moon_offset) / S) * 2 * PI) % (2 * PI)
    diff = moon_long - node
    moon_ecliptic = [moon_long, asin(sin(lunarInclination) * sin(diff)), d]
    return ecliptic_equat(moon_ecliptic)

def calculate_shadow_axis_equat(sun_equat, moon_equat):
    #a
    shadow_axis_ascension = sun_equat[0] - distanceRatio / (1 - distanceRatio) * cos(moon_equat[1]) / cos(sun_equat[1]) * (moon_equat[0] - sun_equat[0])
    #d
    shadow_axis_declension = sun_equat[1] - distanceRatio / (1 - distanceRatio) * (moon_equat[1] - sun_equat[1])
    return [shadow_axis_ascension, shadow_axis_declension, shadow_axis_distance]

def calculate_moon_FundCoords(shadow_axis_equat, moon_equat):
    moon_FundCoords_x = lunar_distance * cos(moon_equat[1]) * sin(moon_equat[0] - shadow_axis_equat[0])
    moon_FundCoords_y = lunar_distance * (sin(moon_equat[1] - shadow_axis_equat[1]) * pow(cos(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2) + sin(moon_equat[1] + shadow_axis_equat[1]) * pow(sin(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2))
    moon_FundCoords_z = lunar_distance * (cos(moon_equat[1] - shadow_axis_equat[1]) * pow(cos(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2) - cos(moon_equat[1] + shadow_axis_equat[1]) * pow(sin(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2))
    return [moon_FundCoords_x, moon_FundCoords_y, moon_FundCoords_z]

def calculate_sidereal(time):
    time_1 = time + 0.5 * (S / (S + 1))
    return ((time_1 * (S + 1) / S) - floor(time_1 * (S + 1) / S)) % 1 * 2 * PI

sunAngular = atan(R/D)
moonAngular = atan(rm/d)
quotient = S/moonPeriod
distanceRatio = d / D #b
lunarParallax_mean = asin(r / d)
solarParallax_mean = asin(r / D)
lunarInclination_prime = atan(quotient / (quotient - 1) * math.tan(lunarInclination))

precessionPeriod = 2 * PI / precession * S

#g
shadow_axis_distance = 1 - distanceRatio
#k
radius_ratio = rm / r

lunar_distance = d / r
newMoon = 0


def draw():
    global time
    global newMoon
    newMoon = time
    
    sinf_exterior = (sin(sunAngular) + radius_ratio * sin(solarParallax_mean)) / shadow_axis_distance #sine of the angle of the cone, penumbra
    sinf_interior = (sin(sunAngular) - radius_ratio * sin(solarParallax_mean)) / shadow_axis_distance #", umbra
    i_exterior = tan(asin(sinf_exterior)) #tangent of the angle of the cone, penumbra
    i_interior = tan(asin(sinf_interior)) #", umbra
    
    time_0 = newMoon
    sun_equat_nm = calculate_sun_equat(time_0)
    moon_equat_nm = calculate_moon_equat(time_0)
    shadow_axis_equat_nm = calculate_shadow_axis_equat(sun_equat_nm, moon_equat_nm)
    moon_FundCoords_nm = calculate_moon_FundCoords(shadow_axis_equat_nm, moon_equat_nm)
    c_exterior_nm = moon_FundCoords_nm[2] + radius_ratio / sinf_exterior 
    c_interior_nm = moon_FundCoords_nm[2] - radius_ratio / sinf_interior
    penumbra_radius_nm = i_exterior * c_exterior_nm #l
    umbra_radius_nm = i_interior * c_interior_nm #l
        
    time_prime = newMoon + 1.0 / earth_hours_in_day
    sun_equat_2 = calculate_sun_equat(time_prime)
    moon_equat_2 = calculate_moon_equat(time_prime)
    shadow_axis_equat_2 = calculate_shadow_axis_equat(sun_equat_2, moon_equat_2)
    moon_FundCoords_2 = calculate_moon_FundCoords(shadow_axis_equat_2, moon_equat_2)
    moon_FundCoords_2_x = moon_FundCoords_2[0]
    moon_FundCoords_2_y = moon_FundCoords_2[1]
    time_prime = newMoon - 1.0 / earth_hours_in_day
    sun_equat_2 = calculate_sun_equat(time_prime)
    moon_equat_2 = calculate_moon_equat(time_prime)
    shadow_axis_equat_2 = calculate_shadow_axis_equat(sun_equat_2, moon_equat_2)
    moon_FundCoords_2 = calculate_moon_FundCoords(shadow_axis_equat_2, moon_equat_2)
    x_prime = (moon_FundCoords_2_x - moon_FundCoords_2[0]) / (1.0 / (earth_hours_in_day/2))
    y_prime = (moon_FundCoords_2_y - moon_FundCoords_2[1]) / (1.0 / (earth_hours_in_day/2))
    
    if(moon_FundCoords_nm[1] == 0):
        moon_FundCoords_nm[1] = 0.0000001
    M_0 = atan2(moon_FundCoords_nm[0] , moon_FundCoords_nm[1])
    if(M_0 == 0):
        M_0 = 0.0000001
    m_0 = moon_FundCoords_nm[0] / sin(M_0)
    N = atan2(x_prime , y_prime)
    n = x_prime / sin(N)
    
    p = 1
    time_0_correction = time_0 - m_0 / n * cos(M_0 - N)
    psi_exterior = asin((m_0 * sin(M_0 - N)) / (p + penumbra_radius_nm))
    offset_exterior = abs((p + penumbra_radius_nm) / n * cos(psi_exterior))
    eclipse_begin_exterior = time_0_correction - offset_exterior
    eclipse_end_exterior = time_0_correction + offset_exterior
    
    psi_interior = asin((m_0 * sin(M_0 - N)) / (p - umbra_radius_nm))
    offset_interior = abs((p - umbra_radius_nm) / n * cos(psi_interior))
    eclipse_begin_interior = time_0_correction - offset_interior
    eclipse_end_interior = time_0_correction + offset_interior
    
    time = eclipse_begin_exterior
    moon_x_memory = 0
    moon_y_memory = 0
    sun_equat = calculate_sun_equat(time)
    moon_equat = calculate_moon_equat(time)
    shadow_axis_equat = calculate_shadow_axis_equat(sun_equat, moon_equat)
    moon_FundCoords = calculate_moon_FundCoords(shadow_axis_equat, moon_equat)
    moon_dist_memory = moon_x_memory * moon_x_memory + moon_y_memory * moon_y_memory + 10
    count = 0
    num = 0
    max_duration = 0
    max_time = 0
    max_lat = -999
    max_long = -999
    print("New Moon: " + str(newMoon))
    background(255)
    if (showGrid):
        drawGrid(newMoon)
    
    while(time <= eclipse_end_exterior):
        moon_x_memory = moon_FundCoords[0]
        moon_y_memory = moon_FundCoords[1]
        
        sun_equat = calculate_sun_equat(time)
        moon_equat = calculate_moon_equat(time)
        shadow_axis_equat = calculate_shadow_axis_equat(sun_equat, moon_equat)
        moon_FundCoords = calculate_moon_FundCoords(shadow_axis_equat, moon_equat)
        
        c_exterior = moon_FundCoords[2] + radius_ratio / sinf_exterior 
        c_interior = moon_FundCoords[2] - radius_ratio / sinf_interior
        penumbra_radius = i_exterior * c_exterior #l
        umbra_radius = i_interior * c_interior #l
        
        sidereal_time = calculate_sidereal(time)
        mu_1 = sidereal_time - shadow_axis_equat[0]
        
        time_2 = time + 0.01
        sun_equat_2 = calculate_sun_equat(time_2)
        moon_equat_2 = calculate_moon_equat(time_2)
        shadow_axis_equat_2 = calculate_shadow_axis_equat(sun_equat_2, moon_equat_2)
        moon_FundCoords_2 = calculate_moon_FundCoords(shadow_axis_equat_2, moon_equat_2)
        derivative_x = (moon_FundCoords_2[0] - moon_FundCoords[0]) / 0.01
        derivative_y = (moon_FundCoords_2[1] - moon_FundCoords[1]) / 0.01
        c_exterior_2 = moon_FundCoords_2[2] + radius_ratio / sinf_exterior 
        c_interior_2 = moon_FundCoords_2[2] - radius_ratio / sinf_interior
        penumbra_radius_2 = i_exterior * c_exterior_2
        umbra_radius_2 = i_interior * c_interior_2
        derivative_l_exterior = (penumbra_radius_2 - penumbra_radius) / 0.01
        derivative_l_interior = (umbra_radius_2 - umbra_radius) / 0.01
        derivative_d = (shadow_axis_equat_2[1] - shadow_axis_equat[1]) / 0.01
        sidereal_time_2 = calculate_sidereal(time_2)
        derivative_mu = (sidereal_time_2 - shadow_axis_equat_2[0] - mu_1) / 0.01
        
        a_prime = - derivative_l_exterior - derivative_mu * i_exterior * moon_FundCoords[0] * cos(shadow_axis_equat[1])
        b_prime = - derivative_y + derivative_mu * moon_FundCoords[0] * sin(shadow_axis_equat[1])
        c_prime = derivative_x + derivative_mu * moon_FundCoords[1] * sin(shadow_axis_equat[1]) + derivative_mu * i_exterior * penumbra_radius * cos(shadow_axis_equat[1])
        
        E = atan2(b_prime , c_prime)
        e = b_prime / sin(E)
        F = atan2(derivative_d , (derivative_mu * cos(shadow_axis_equat[1])))
        f = derivative_d / sin(F)
        
        if(moon_FundCoords[1] == 0):
            moon_FundCoords[1] = 0.0000001
        M = atan2(moon_FundCoords[0] , moon_FundCoords[1])
        if(M == 0):
            M = 0.0000001
        m = moon_FundCoords[0] / sin(M)
        if(m == 0):
            m += 0.000001
            
        #PENUMBRA
        lam = sqrt((penumbra_radius + m - p) * (penumbra_radius - m + p) / (4 * m))
        lamb = 2 * asin(lam)
        if(lamb > PI/2):
            lamb = PI - lamb
            
        #CURVES OF PENUMBRA CONTACT AT HORIZON
        gamma = M + lamb
        col = color(0, 0, 255)
        if(m * sin(M - E) < sin(gamma - E)):
            col = color(255, 0, 0)
        calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col)
    
        gamma = M - lamb
        col = color(0, 0, 255)
        if(m * sin(M - E) < sin(gamma - E)):
            col = color(255, 0, 0)
        calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col)
        
        #UMBRA
        lam = sqrt((umbra_radius + m - p) * (umbra_radius - m + p) / (4 * m))
        lamb = 2 * asin(lam)
        if(lamb > PI/2):
            lamb = PI - lamb
            
        #CURVES OF UMBRA CONTACT AT HORIZON
        gamma = M + lamb
        col = color(0, 0, 255)
        if(m * sin(M - E) < sin(gamma - E)):
            col = color(255, 0, 0)
        calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col)
    
        gamma = M - lamb
        col = color(0, 0, 255)
        if(m * sin(M - E) < sin(gamma - E)):
            col = color(255, 0, 0)
        calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col)
        
        #CURVE OF MAXIMUM AT HORIZON
        psi = asin(m * sin(M - E))
        Delta = abs(m * cos(M - E) - cos(psi))
        psi_2 = PI - asin(m * sin(M - E))
        Delta_2 = abs(m * cos(M - E) - cos(psi_2))
        if(Delta < penumbra_radius):
            gamma = E + psi
            col = color(255, 0, 255)
            calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col)    
        if(Delta_2 < penumbra_radius):
            gamma = E + psi_2
            col = color(255, 0, 255)
            calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col)    
        
        #SOUTHERN LIMIT OF PARTIALITY
        nu = atan2(f, e)
        psi = atan(tan(PI/4 + nu) * tan(E/2))
        Q_1 = E
        Q_2 = E/2 + psi
        Q = (Q_1 + Q_2)/2
        Q = correctQ(Q, moon_FundCoords, penumbra_radius, f, e, E, False)
        calculateCoordFromQ(Q, moon_FundCoords, i_exterior, shadow_axis_equat, mu_1, penumbra_radius, color(0, 0, 0), 2)
        
        #NORTHERN LIMIT OF PARTIALITY
        Q_1 = E + PI
        Q_2 = E/2 + psi + PI
        Q = (Q_1 + Q_2)/2
        Q = correctQ(Q, moon_FundCoords, penumbra_radius, f, e, E, True)
        calculateCoordFromQ(Q, moon_FundCoords, i_exterior, shadow_axis_equat, mu_1, penumbra_radius, color(0, 0, 0), 2)
        
        #PATH OF CENTRALITY
        a = moon_FundCoords[0]
        b = moon_FundCoords[1]
        gamma = atan2(a, b)
        beta = asin(a / sin(gamma))
        C = atan2(moon_FundCoords[1], cos(beta))
        c = moon_FundCoords[1] / sin(C)
        latitude = asin(c * sin(C + shadow_axis_equat[1]))
        sin_theta = moon_FundCoords[0] / cos(latitude)
        cos_theta = c * cos(C + shadow_axis_equat[1]) / cos(latitude)
        theta = atan2(sin_theta, cos_theta)
        longitude = mu_1 - theta
        col = color(0, 150, 0)
        if(umbra_radius > 0):
            col = color(150, 0, 0)
        drawPoint(longitude, latitude, col, 2)
        L = abs(umbra_radius - i_interior * cos(beta))
        a = c_prime - f * cos(beta)
        Q = atan2(a, b_prime)
        t = earth_hours_in_day * 60 * L * sin(Q) / a
        #MAXIMUM ECLIPSE AT NOON
        if(sign(moon_x_memory) != sign(moon_FundCoords[0])):
            drawPoint(mu_1, asin(moon_FundCoords[1]) + shadow_axis_equat[1], col, 10)
            print("Maximum Eclipse at Noon: " + str(-mu_1 * 180 / PI) + "E, " + str((asin(moon_FundCoords[1]) + shadow_axis_equat[1]) * 180 / PI) + "N")
        if(t > max_duration):
            max_duration = t
        if(moon_FundCoords[0] * moon_FundCoords[0] + moon_FundCoords[1] * moon_FundCoords[1] <= moon_dist_memory):
            max_time = time - 0.00005
            max_lat = latitude
            max_long = longitude
            moon_dist_memory = moon_FundCoords[0] * moon_FundCoords[0] + moon_FundCoords[1] * moon_FundCoords[1]
        
        #SOUTHERN LIMIT OF TOTALITY
        nu = atan2(f, e)
        psi = atan(tan(PI/4 + nu) * tan(E/2))
        Q_1 = E
        Q_2 = E/2 + psi
        Q = (Q_1 + Q_2)/2
        Q = correctQ(Q, moon_FundCoords, umbra_radius, f, e, E, False)
        col = color(0, 150, 0)
        if(umbra_radius > 0):
            col = color(150, 0, 0)
        calculateCoordFromQ(Q, moon_FundCoords, i_interior, shadow_axis_equat, mu_1, umbra_radius, col, 2)
        
        #NORTHERN LIMIT OF TOTALITY
        Q_1 = E + PI
        Q_2 = E/2 + psi + PI
        Q = (Q_1 + Q_2)/2
        Q = correctQ(Q, moon_FundCoords, umbra_radius, f, e, E, True)
        col = color(0, 150, 0)
        if(umbra_radius > 0):
            col = color(150, 0, 0)
        calculateCoordFromQ(Q, moon_FundCoords, i_interior, shadow_axis_equat, mu_1, umbra_radius, col, 2)
        
        if(abs((count - 0.01) - num * (eclipse_end_interior - eclipse_begin_interior) / 6) < 0.001 and time < eclipse_end_interior and time > eclipse_begin_interior):
            for i in list(range(200)):
                Q = 2 * PI / 200 * i
                col = color(0, 200, 0)
                if(umbra_radius > 0):
                    col = color(200, 0, 0)
                calculateCoordFromQ(Q, moon_FundCoords, i_interior, shadow_axis_equat, mu_1, umbra_radius, col, 1.5)
            print("Moon Shadow " + str(num + 1) + ": " + str(time))
            num += 1
        
        if (time < eclipse_end_exterior):
            time += 0.0001
            if(time < eclipse_end_interior and time > eclipse_begin_interior):
                count += 0.0001

    col = color(0, 200, 0)
    if(umbra_radius > 0):
        col = color(200, 0, 0)
    if(max_long != -999):
        drawPoint(max_long, max_lat, col, 10)
    print("First Exterior Contact: " + str(eclipse_begin_exterior))
    print("First Interior Contact: " + str(eclipse_begin_interior))
    print("Maximum Eclipse: " + str(max_time))
    print("Maximum Eclipse Location: " + str(-max_long * 180 / PI) + "E, " +  str(max_lat * 180 / PI) + "N")
    print("Last Interior Contact: " + str(eclipse_end_interior))
    print("Last Exterior Contact: " + str(eclipse_end_exterior))
    print("Max Duration: " + str(max_duration) + " Minutes")
    noLoop()
    saveFrame(str(newMoon) + ".png")
    
def drawPoint(longitude, latitude, col, thick):
    if(longitudeCorrect):
        correct = 0.5 - (newMoon % 1)
        longitude += correct * (2 * PI)
    longitude = longitude % (2 * PI)
    latitude = (latitude + PI/2) % PI - PI/2
    stroke(col)
    strokeWeight(thick)
    x = width / (2 * PI) * -longitude + width / 2
    if(x < 0):
        x += width
    y = height - height / PI * (latitude + PI/2)
    point(x, y)

def drawGrid(newMoon):
    for i in list(range(36)):
        stroke(0)
        strokeWeight(0.25)
        if(i == 18):
            stroke(255, 0, 0)
            strokeWeight(0.75)
        elif(i == 0):
            strokeWeight(0.75)
        if(longitudeCorrect):
            correct = 0.5 - (newMoon % 1)
            x = (i * width/36 - correct * width) % width
            line(x, 0, x, height)
        else:
            line(i * width/36, 0, i * width/36, height)
    for i in list(range(18)):
        strokeWeight(0.25)
        if(i == 9):
            strokeWeight(0.75)
        stroke(0)
        line(0, i * height/18, width, i * height/18)
    
def calculateCoordFromGamma(gamma, shadow_axis_equat, mu_1, col):
    latitude = asin(cos(gamma) * cos(shadow_axis_equat[1]))
    sin_theta = sin(gamma) / cos(latitude)
    cos_theta = -cos(gamma) * sin(shadow_axis_equat[1]) / cos(latitude)
    theta = atan2(sin_theta,cos_theta)
    longitude = mu_1 - theta
    drawPoint(longitude, latitude, col, 2)
    
def calculateCoordFromQ(Q, moon_FundCoords, i, shadow_axis_equat, mu_1, radius, col, n):
    a = moon_FundCoords[0] - radius * sin(Q)
    b = moon_FundCoords[1] - radius * cos(Q)
    gamma = atan2(a, b)
    beta = asin(a / sin(gamma))
    epsilon = i * cos(Q - gamma)
    zeta = cos(beta + epsilon)
    xi = a + i * zeta * sin(Q)
    eta = b + i * zeta * cos(Q)
    C = atan2(eta, zeta)
    c = eta / sin(C)
    latitude = asin(c * sin(C + shadow_axis_equat[1]))
    sin_theta = xi / cos(latitude)
    cos_theta = c * cos(C + shadow_axis_equat[1]) / cos(latitude)
    theta = atan2(sin_theta, cos_theta)
    longitude = mu_1 - theta
    drawPoint(longitude, latitude, col, n)
    
def correctQ(Q, moon_FundCoords, radius, f, e, E, north):
    a = moon_FundCoords[0] - radius * sin(Q)
    b = moon_FundCoords[1] - radius * cos(Q)
    gamma = atan2(a, b)
    beta = a / sin(gamma)
    nu_prime = atan(f / e * cos(beta))
    Q = atan(tan(PI/4 + nu_prime) * tan(E/2)) + E/2
    if(north):
        return Q + PI
    return Q

def sign(x):
    if(x == 0):
        return 0
    return x / abs(x)
