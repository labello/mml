
Na                       = 6.0221412927 * (10**23)  # Avogadro's Number

k = 1.3806488 * 10**-23   # Boltzman's constants kg * m**2
                          #                      ---------
                          #                       K * s**2

SECONDS_TO_FEMPTOSECONDS = 1.0 * 10**15
SECONDS_TO_PICOSECONDS   = 1.0 * 10**12

FEMPTOSECONDS_TO_SECONDS = 1.0 / SECONDS_TO_FEMPTOSECONDS
PICOSECONDS_TO_SECONDS   = 1.0 / SECONDS_TO_PICOSECONDS

AMU_TO_KG                = 1.66053892 * (10**-27)
KG_TO_AMU                = 1.0/AMU_TO_KG

ANGSTROMS_TO_METERS      = 10**-10
METERS_TO_ANGSTROMS      = 10**10 


# k2 is Boltzman's constant in   AMU * ANGSROMS**2
#                              ---------------------
#                               K * FEMPTOSECOND**2  
k2 = k * KG_TO_AMU * METERS_TO_ANGSTROMS**2 * (1/SECONDS_TO_FEMPTOSECONDS)**2 

# k2 is Boltzman's constant in   AMU * ANGSROMS**2
#                              ---------------------
#                               K * PICOSECOND**2  
k3 = k * KG_TO_AMU * METERS_TO_ANGSTROMS**2 * (1/SECONDS_TO_PICOSECONDS)**2
