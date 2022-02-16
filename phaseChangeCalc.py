from sympy.solvers import solve


# --- Lambda Methods ---
def sign(x):
    return (int(abs(x)/x)) if (x != 0) else 1


def diff(x):
    return (x[-1] - x[0]) if (len(x) > 0) else 0


# --- Errors ---
class Error(Exception):
    pass


class NoValidInformation(Error):
    pass


# --- Compounds ---
class Compound(object):
    molar = None
    vapor = None
    fusion = None

    melt = None
    boil = None

    Cp_g = None
    Cp_l = None
    Cp_s = None

    Cp_a = None


class Water(Compound):
    molar = 18.015
    vapor = 40.7
    fusion = 6.02

    melt = 0
    boil = 100

    Cp_g = 1.998
    Cp_l = 4.184
    Cp_s = 2.03


# --- Methods ---
def mCAT(Q="q", M="m", Cp="C", T="T", Print=True):
    solution = solve("{Mass} * {Cap} * {Temp} - {Qenergy}".format(Mass=M, Cap=Cp, Temp=T, Qenergy=Q))
    simplified = [".".join([str(i) for i in str(x).split(".") if int(i) != 0]) for x in solution]

    #print(simplified)

    solution = (simplified) if (len(simplified) != 0 and simplified[0] != '') else ([0])

    if solution[0] != 0 and Print: print("{Qenergy} = {Mass} * {Cap} * {Temp} ({SOLVE})".format(Mass=M, Cap=Cp, Temp=T, Qenergy=Q, SOLVE=float(solution[0])))
    return float(solution[0])


def molVapor(Q="q",mol="m",vapor="v",molar=None,mass="m",Print=True):
    if molar != None:
        mol = mass / molar
    solution = solve("{MOLE} * {VAPOR} * 1000 - {ENERGY}".format(MOLE=mol, VAPOR=vapor, ENERGY=Q))
    solution = [".".join([str(i) for i in str(x).split(".") if int(i) != 0]) for x in solution]

    if Print: print("{ENERGY} = {MOLE:.3f} * {VAPOR} * 1000 ({SOLVE})".format(MOLE=mol, VAPOR=vapor, ENERGY=Q, SOLVE=float(solution[0])))
    return float(solution[0])


def molFusion(Q="q",mol="m",fusion="v",molar=None,mass="m",Print=True):
    if molar != None:
        mol = mass / molar
    solution = solve("{MOLE} * {FUSION} * 1000 - {ENERGY}".format(MOLE=mol, FUSION=fusion, ENERGY=Q))
    solution = [".".join([str(i) for i in str(x).split(".") if int(i) != 0]) for x in solution]

    if Print: print("{ENERGY} = {MOLE:.3f} * {FUSION} * 1000 ({SOLVE})".format(MOLE=mol, FUSION=fusion, ENERGY=Q, SOLVE=float(solution[0])))
    return float(solution[0])


def getDistances(I, F):
    returnVal = []
    Temps = [I, F]

    if (max(Temps) == min(Temps)): return [0,0,0]

    elif (max(Temps) <= 0):
        returnVal = [abs(F - I), 0, 0]
    elif (min(Temps) >= 100):
        returnVal = [0,0,abs(F - I)]
    elif (min(Temps) >= 0 and max(Temps) <= 100):
        returnVal = [0,abs(F - I),0]

    elif (min(Temps) <= 0 and max(Temps) <= 100):
        returnVal = [abs(min(Temps)),max(Temps),0]
    elif (min(Temps) >= 0 and max(Temps) >= 100):
        returnVal = [0, 100-min(Temps), max(Temps)-100]

    elif (min(Temps) <= 0 and max(Temps) >= 100):
        returnVal = [abs(min(Temps)), 100, (max(Temps)-100)]

    return [(i*abs(F - I) / (F - I)) for i in returnVal]


def pCC(compound=Water, mass="m", moles="M", Temp_Init=None, Temp_Final=None, Energy="q", LSconv=False, LGconv=False, Print=True):

    try:

        if compound.Cp_a is None:
            compound.Cp_g = compound.Cp_a
            compound.Cp_l = compound.Cp_a
            compound.Cp_s = compound.Cp_a

        TotalEnergy = 0

        if (mass == "m" and moles == "M" and Energy == "q") or (Temp_Init is None or Temp_Final is None):
            raise NoValidInformation

        elif mass != "m" and moles == "M" and compound.molar is None:
            moles = mass / compound.molar

        elif moles != "M" and mass == "m" and compound.molar is None:
            mass = moles * compound.molar

        Temp_Change = getDistances(Temp_Init, Temp_Final)

        CpE = [
            mCAT(Q=Energy, M=mass, Cp=compound.Cp_g, T=Temp_Change[0],Print=Print),  # x > 100
            mCAT(Q=Energy, M=mass, Cp=compound.Cp_l, T=Temp_Change[1],Print=Print),  # x > 0 and x < 100
            mCAT(Q=Energy, M=mass, Cp=compound.Cp_s, T=Temp_Change[2],Print=Print)   # x <= 0
             ]

        TotalEnergy += sum(CpE)

        if (CpE[1] != 0 and CpE[2] != 0) or (LGconv):
            if Print:
                print("")
            LG = molVapor(Q=Energy, mol=moles, vapor=compound.vapor, Print=Print) * sign(Temp_Final - Temp_Init)
            TotalEnergy += LG
            if Print:
                print("\tLiquid <-> Gas (VAPOR):", LG)

        if (CpE[1] != 0 and CpE[0] != 0) or LSconv:
            if Print:
                print("")
            LS = molFusion(Q=Energy, mol=moles, fusion=compound.fusion, Print=Print) * sign(Temp_Final - Temp_Init)
            TotalEnergy += LS
            if Print:
                print("\tLiquid <-> Solid (FUSION):", LS)

        return TotalEnergy

    except NoValidInformation:
        print("")
        print("pCC: No Valid Inputs were entered (MOLES or MASS or TEMP)")
        print("")