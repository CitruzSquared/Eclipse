import math
R = 673320.0 #Radius of the sun in km
r = 5512.0 #Radius of the planet in km
rm = 2260.113 #Radius of the moon in km
D = 126096000.0 #Distance to the sun in km
d = 294000.241 #Distance to the moon in km
S = 298.0 #Planetary orbital period in days
M = 22.4435146444 #Lunar sidereal orbital period in days
precession = math.radians(15.42840432988) #Lunar nodal precession rate in radians per year
tilt = math.radians(24.7) #Planetary axial tilt in radians
lunarInclination = math.radians(4.92) #Inclination of the Moon's orbit in radians
moon_offset = -78 #Can be any number
node_offset = -78 #Can be any number, doesn't have to be equal to moon_offset

###DO NOT TOUCH FROM HERE BELOW###
sunAngular = math.atan(R/D)
moonAngular = math.atan(rm/d)
earthUmbra = math.atan((R - (R - r) * (D + d)/D)/d)
earthPenumbra = math.atan(((R + r) * (D + d)/D - R)/d)
quotient = S/M
precessionPeriod = 2 * math.pi / precession * S
lunarParallax_mean = math.asin(r / d)
solarParallax_mean = math.asin(r / D)
lunarInclination_prime = math.atan(quotient / (quotient - 1) * math.tan(lunarInclination))
solarEclipseLimit = (lunarParallax_mean - solarParallax_mean + moonAngular + sunAngular) * 1 / math.cos(lunarInclination_prime)
U = abs((earthUmbra - moonAngular))
P = (moonAngular + earthUmbra)
R = (moonAngular + earthPenumbra) 
precessionPerDay = precession / S
draconic = 1/((2 * math.pi / M + precessionPerDay)/(2 * math.pi))
month = S * M/abs(S - M)
n = 0
i = 0
node1 = 0
node2 = node1 + math.pi
k = 0
error = 0.001
diff = 5
j = 0
PI = math.pi
firstNewMoon = moon_offset
firstFullMoon = moon_offset + month/2
while(diff > error):
    j += 1
    n = month * j / draconic
    diff = abs(n - round(n))
saros = j * month

sarosmode = float(input("Saros mode (1 for yes)? "))
eclipsetype_sarosmode = 0
if(sarosmode == 1):
    firstEclipse = float(input("Enter date of member eclipse: ")) - moon_offset
    x = firstEclipse / month
    if(abs(x - round(x)) >= 0.1):
        eclipsetype_sarosmode = 1
    firstEclipse += moon_offset
    firstEclipse -= saros * 150
    k = 300
else:    
    beginday = int(input("Begin predicting from what date?: "))
    if(firstNewMoon < beginday):
        while(firstNewMoon < beginday):
            firstNewMoon += month
    else:
        while(firstNewMoon > beginday):
            firstNewMoon -= month
    if(firstNewMoon - month/2 < beginday):
        firstFullMoon = firstNewMoon + month/2
    else:
        firstFullMoon = firstNewMoon - month/2
    k = int(input("Predict how many eclipses?: "))
day = 0
n = 0

def ecliptic_equat(coords):
    ascension = math.atan2(math.cos(tilt) * math.sin(coords[0]) - math.sin(tilt) * math.tan(coords[1]), math.cos(coords[0]))
    declension = math.asin(math.cos(tilt) * math.sin(coords[1]) + math.sin(tilt) * math.cos(coords[1]) * math.sin(coords[0]))
    if (ascension < 0):
        ascension += 2 * math.pi
    return [ascension, declension, coords[2]]

def moon_ecliptic(day):
    node = (((-(day - node_offset) % precessionPeriod) / precessionPeriod + (S + node_offset) / S) * 2 * PI) % (2 * PI)
    moon_long = ((((day - moon_offset) % M) / M + (S + moon_offset) / S) * 2 * PI) % (2 * PI)
    diff = moon_long - node
    return [moon_long, math.asin(math.sin(lunarInclination) * math.sin(diff)), d]

def sun_ecliptic(day):
    return [((day % S) / S * 2 * PI) % (2 * PI), 0, D]

def calculateLunarEclipseDistance(day):
    moonEclipticPosition = moon_ecliptic(day)
    shadowPosition = sun_ecliptic(day)[0] + PI
    deltaLongitude = (abs(moonEclipticPosition[0] - shadowPosition) % (2 * math.pi))
    distance = abs((math.acos(math.cos(moonEclipticPosition[1]) * math.cos(deltaLongitude)) + math.pi) % (2 * math.pi) - math.pi)
    return distance

def calculateSolarEclipseDistance(day):
    moonEclipticPosition = moon_ecliptic(day)
    sunPosition = sun_ecliptic(day)[0]
    deltaLongitude = (abs(moonEclipticPosition[0] - sunPosition) % (2 * math.pi))
    distance = abs((math.acos(math.cos(moonEclipticPosition[1]) * math.cos(deltaLongitude)) + math.pi) % (2 * math.pi) - math.pi)
    return distance

distanceRatio = d / D
lunar_distance = d / r
shadow_axis_distance = 1 - distanceRatio
def calculate_shadow_axis_equat(sun_equat, moon_equat):
    #a
    shadow_axis_ascension = sun_equat[0] - distanceRatio / (1 - distanceRatio) * math.cos(moon_equat[1]) / math.cos(sun_equat[1]) * (moon_equat[0] - sun_equat[0])
    #d
    shadow_axis_declension = sun_equat[1] - distanceRatio / (1 - distanceRatio) * (moon_equat[1] - sun_equat[1])
    return [shadow_axis_ascension, shadow_axis_declension, shadow_axis_distance]
    
lunar_distance = d / r
radius_ratio = rm / r
def calculate_moon_FundCoords(shadow_axis_equat, moon_equat):
    moon_FundCoords_x = lunar_distance * math.cos(moon_equat[1]) * math.sin(moon_equat[0] - shadow_axis_equat[0])
    moon_FundCoords_y = lunar_distance * (math.sin(moon_equat[1] - shadow_axis_equat[1]) * pow(math.cos(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2) + math.sin(moon_equat[1] + shadow_axis_equat[1]) * pow(math.sin(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2))
    moon_FundCoords_z = lunar_distance * (math.cos(moon_equat[1] - shadow_axis_equat[1]) * pow(math.cos(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2) - math.cos(moon_equat[1] + shadow_axis_equat[1]) * pow(math.sin(0.5 * (moon_equat[0] - shadow_axis_equat[0])), 2))
    return [moon_FundCoords_x, moon_FundCoords_y, moon_FundCoords_z]

while i < k:
    if(sarosmode == 1):
        day = firstEclipse + saros * n
        count = i
    else:
        day = n * month + firstNewMoon
    if((sarosmode == 1 and eclipsetype_sarosmode == 0) or sarosmode != 1):
        distance = calculateSolarEclipseDistance(day)
        dayMemory = day
        if(distance < solarEclipseLimit):
            moon_equat = ecliptic_equat(moon_ecliptic(day))
            shadow_axis_equat = calculate_shadow_axis_equat(ecliptic_equat(sun_ecliptic(day)), moon_equat)
            moon_FundCoords = calculate_moon_FundCoords(shadow_axis_equat, moon_equat)
            sinf_exterior = (math.sin(sunAngular) + radius_ratio * math.sin(solarParallax_mean)) / shadow_axis_distance
            sinf_interior = (math.sin(sunAngular) - radius_ratio * math.sin(solarParallax_mean)) / shadow_axis_distance
            i_exterior = math.tan(math.asin(sinf_exterior))
            i_interior = math.tan(math.asin(sinf_interior))
            penumbra_radius = i_interior * (moon_FundCoords[2] + radius_ratio / sinf_exterior )
            umbra_radius = i_interior * (moon_FundCoords[2] - radius_ratio / sinf_interior)
            day_minus = day - 0.05
            moon_equat_minus = ecliptic_equat(moon_ecliptic(day_minus))
            shadow_axis_equat_minus = calculate_shadow_axis_equat(ecliptic_equat(sun_ecliptic(day_minus)), moon_equat_minus)
            moon_FundCoords_minus = calculate_moon_FundCoords(shadow_axis_equat_minus, moon_equat_minus)
            day_plus = day + 0.05
            moon_equat_plus = ecliptic_equat(moon_ecliptic(day_plus))
            shadow_axis_equat_plus = calculate_shadow_axis_equat(ecliptic_equat(sun_ecliptic(day_plus)), moon_equat_plus)
            moon_FundCoords_plus = calculate_moon_FundCoords(shadow_axis_equat_plus, moon_equat_plus)
            x_prime = (moon_FundCoords_plus[0] - moon_FundCoords_minus[0]) / 0.1
            y_prime = (moon_FundCoords_plus[1] - moon_FundCoords_minus[1]) / 0.1
            if(moon_FundCoords[1] == 0):
                moon_FundCoords[1] = 0.0000001
            M_0 = math.atan2(moon_FundCoords[0] , moon_FundCoords[1])
            if(M_0 == 0):
                M_0 = 0.0000001
            m_0 = moon_FundCoords[0] / math.sin(M_0)
            N = math.atan2(x_prime , y_prime)
            n_calc = x_prime / math.sin(N)
            time_0_correction = day - m_0 / n_calc * math.cos(M_0 - N)
            sin_psi_exterior = (m_0 * math.sin(M_0 - N)) / (1 + penumbra_radius)
            sin_psi_interior = (m_0 * math.sin(M_0 - N)) / (1 - umbra_radius)
            if(abs(sin_psi_interior) > 1 and abs(sin_psi_exterior) <= 1):
                psi_exterior = math.asin(sin_psi_exterior)
                offset_exterior = abs((1 + penumbra_radius) / n_calc * math.cos(psi_exterior))
                contacts_list = [round((time_0_correction - offset_exterior) * 10000) / 10000, "-", round(time_0_correction * 10000) / 10000, "-", round((time_0_correction + offset_exterior) * 10000) / 10000]
                print("PS: " + str(contacts_list) + ", NM: " + str(day))
            elif(abs(sin_psi_interior) <= 1):
                psi_exterior = math.asin(sin_psi_exterior)
                offset_exterior = abs((1 + penumbra_radius) / n_calc * math.cos(psi_exterior))
                psi_interior = math.asin(sin_psi_interior)
                offset_interior = abs((1 - umbra_radius) / n_calc * math.cos(psi_interior))
                contacts_list = [round((time_0_correction - offset_exterior) * 10000) / 10000, round((time_0_correction - offset_interior) * 10000) / 10000, round(time_0_correction * 10000) / 10000, round((time_0_correction + offset_interior) * 10000) / 10000, round((time_0_correction + offset_exterior) * 10000) / 10000]
                print("TS: " + str(contacts_list) + ", NM: " + str(day))
            i += 1
        if i == k:
            break
    if((sarosmode == 1 and eclipsetype_sarosmode == 1) or sarosmode != 1):
        if(sarosmode != 1):
            day = n * month + firstFullMoon
        distance = calculateLunarEclipseDistance(day)
        dayMemory = day
        if (distance <= math.degrees(R) + 1):
            day -= 5/24
            eclipseState = 0
            eclipseStateMemory = 0
            maxEclipseState = 0
            contacts = ["-", "-", "-", "-", "-", "-", "-"]
            distanceMemory = 0
            distance = calculateLunarEclipseDistance(day)
            greatest = 0
            magnitude = 0
            while(day <= dayMemory + 5/24):
                eclipseStateMemory = eclipseState
                distanceMemory = distance
                distance = calculateLunarEclipseDistance(day)
                if(distance <= U):
                    eclipseState = 3
                    if(maxEclipseState < eclipseState):
                        maxEclipseState = eclipseState
                elif(distance <= P):
                    eclipseState = 2
                    if(maxEclipseState < eclipseState):
                        maxEclipseState = eclipseState
                elif(distance <= R):
                    eclipseState = 1
                    if(maxEclipseState < eclipseState):
                        maxEclipseState = eclipseState
                elif(distance > R):
                    eclipseState = 0
                if(distance > distanceMemory and greatest == 0):
                    contacts[3] = float(str(round((day - 1/2880) * 10000) / 10000).zfill(4))
                    distance = calculateLunarEclipseDistance(day - 1/2880)
                    if(eclipseState > 1):
                        magnitude = round(((moonAngular - (distance - earthUmbra)) / (2 * moonAngular) * 100) * 10000) / 10000
                    else:
                        magnitude = round(((moonAngular - (distance - earthPenumbra)) / (2 * moonAngular) * 100) * 10000) / 10000
                    greatest = 1
                if(eclipseState > eclipseStateMemory):
                    contacts[int(eclipseState - 1)] = float(str(round((day) * 10000) / 10000).zfill(4))
                elif(eclipseState < eclipseStateMemory):
                    contacts[int(6 - eclipseState)] = float(str(round((day) * 10000) / 10000).zfill(4))
                day += 1/1440
            day = dayMemory
            if(maxEclipseState != 0):
                eclipseType = " "
                if(maxEclipseState == 1):
                    eclipseType = "R"
                elif(maxEclipseState == 2):
                    eclipseType = "P"
                else:
                    eclipseType = "T"
                print(eclipseType + "L: " + str(contacts) + ", FM: " + str(day) + ", D: " + str(round((contacts[6] - contacts[0]) * 10000) / 10000) + ", M: " + str(magnitude))
                i += 1            
    if(sarosmode == 1 and i != 0 and count == i):
        break
    n += 1
if(sarosmode == 1):
    print(str(i) + " Members")
print("Sun Angular Radius: " + str(sunAngular * 180 / math.pi))
print("Moon Angular Radius: " + str(moonAngular * 180 / math.pi))
print("Earth Umbra Angular Radius: " + str(earthUmbra * 180 / math.pi))
print("Earth Penumbra Angular Radius: " + str(earthPenumbra * 180 / math.pi))
print("TS : Total solar eclipse")
print("PS : Partial solar eclipse")
print("UL : Umbral lunar eclipse")
print("PL : Partial lunar eclipse")
print("RL : Penumbral lunar eclipse")
print("Synodic Month : " + str(month))
print("Draconic Month: " + str(draconic))
print("Saros: " + str(j) + " Synodic months")
print("Saros: " + str(saros) + " Days")
