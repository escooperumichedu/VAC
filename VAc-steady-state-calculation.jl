using NLsolve, ComponentArrays

# Nonlinear Large-Scale Process Network for Vinyl Acetate (VAc) Production

# Define nominal inputs
Q_H1 = 5078.69 # kcal/min
Q_H2 = 1461.14 # kcal/min
Q_C3 = 15491.57 # kcal/min
Q_C4 = 7250.42 # kcal/min
Q_C5 = 1881.2 # kcal/min

f_S4 = 12.113916 # kmol/min
f_S5 = 0.47744 # kmol/min
f_S13 = 2.73964 # kmol/min (liquid)
f_S17 = 1.20871 # kmol/min (liquid)
f_S25 = 0.7443 # kmol/min

Q_VAP = 16933.247 # kcal/min
T_RCT_coolant = 133.46 # °C
Ws_COM = 275.64 # kcal/min
T_SEP_coolant = 37.72 # °C

f_S1 = 0.905 # kmol/min
f_S18 = 15.358 # kmol/min (liquid)
f_S22 = 0.8125 # kmol/min (liquid)
f_S3 = 2.1924 # kmol/min (liquid)
f_S27 = 6.556 # kmol/min
f_S30 = 0.00318 # kmol/min

# Steady state flow values given dM/dt = 0 for all unit ops
f_S2 = f_S4 - f_S3
f_S6 = f_S4 # Same flowrate, just heated
f_S7 = f_S6 + f_S5 # Mixing streams as the reactor inlet
f_S8 = f_S7 # hot-stream inlet to HX
f_S9 = f_S8 # cold-stream outlet from HX
f_S10 = f_S9 # Presure letdown to precooler inlet
f_S11 = f_S10 # precooler outlet to reactor inlet

f_S34 = f_S2 - f_S1 # cold-stream outlet from HX
f_S33 = f_S34 # hot-stream inlet from HX

f_S12 = f_S11 - f_S13 # vapor outlet from seperator
f_S14 = f_S12 # Compressed vapor separator outlet
f_S15 = f_S14 # Heated vapor absorber inlet

f_S23 = f_S22 # HAc recycle to absorber inlet
f_S21 = f_S13 + f_S17 # Inlet to the Azeotropic distillation column

f_S19 = f_S17 + f_S18 # Bottom's outlet from the absorber
f_S20 = f_S18 # Cooled feed to absorber

f_S16 = f_S15 + f_S20 + f_S23 - f_S19
f_S28 = f_S16 - f_S27

f_S31 = f_S33 - f_S28
f_S29 = f_S30 + f_S31

f_S32 = f_S27 - f_S29

f_S24 = f_S3 + f_S22 - f_S25

f_S26 = f_S21 - f_S24

# Thermodynamics and physical data
function thermo_data()

    return Dict(
        "O2" => (MW=32.000, SpG=0.50, h_L=2300, a_liq=0.3, b_liq=0, a=0.218, b=0.0001, v=64.178, A=9.2, B=0, C=273),
        "CO2" => (MW=44.010, SpG=1.18, h_L=2429, a_liq=0.60, b_liq=0, a=0.230, b=0, v=37.4, A=7.937, B=0, C=273),
        "C2H4" => (MW=28.052, SpG=0.57, h_L=1260, a_liq=0.60, b_liq=0, a=0.370, b=0.0007, v=49.347, A=9.497, B=-313, C=273),
        "C2H6" => (MW=30.068, SpG=0.57, h_L=1260, a_liq=0.60, b_liq=0, a=0.370, b=0.0007, v=52.866, A=9.497, B=-313, C=273),
        "VAc" => (MW=86.088, SpG=0.85, h_L=8600, a_liq=0.44, b_liq=0.0011, a=0.290, b=0.0006, v=101.564, A=12.6564, B=-2984.45, C=226.66),
        "H2O" => (MW=18.008, SpG=1.00, h_L=10684, a_liq=0.99, b_liq=0.0002, a=0.560, b=-0.0016, v=18.01, A=14.6394, B=-3984.92, C=233.426),
        "HAc" => (MW=60.052, SpG=0.98, h_L=5486, a_liq=0.46, b_liq=0.0012, a=0.520, b=0.0007, v=61.445, A=14.5236, B=-4457.83, C=258.45))

end

# Vapor heat capacity
function cp_vap(T, comp)
    a = comp.a
    b = comp.b
    MW = comp.MW
    return (a + b * T) * MW
end

# Liquid heat capacity
function cp_liq(T, comp)
    a = comp.a_liq
    b = comp.b_liq
    MW = comp.MW
    return (a + b * T) * MW
end

# Vapor enthalpy
function h_vap(T, comp)
    a = comp.a
    b = comp.b
    MW = comp.MW
    h_L = comp.h_L
    return (a * T + 0.5 * b * T^2) * MW + h_L
end

# Liquid enthalpy
function h_liq(T, comp)
    a = comp.a_liq
    b = comp.b_liq
    MW = comp.MW
    return (a * T + 0.5 * b * T^2) * MW
end

# Antoine equation for saturation pressure
function p_sat(T, comp)
    A = comp.A
    B = comp.B
    C = comp.C
    return exp((B / (T + C)) + A)
end

# Example usage of the thermo data dict
thermo = thermo_data()

# Ethylene feed (f_S1) with inert ethane
y_1_C2H6 = 0.001
y_1_C2H4 = 1 - y_1_C2H6
T_1 = 30 # °C
P_1 = 150 # psia

y_5_O2 = 1.0 # Pure oxygen feed (f_S5)
T_5 = 30 # °C
P_5 = 150 # psia

# Assume all vessels are at 50% maximum allowable volume
M_VAP = 4 * 0.5 # m3
M_SEP = 8 * 0.5 # m3
M_TK = 2.83 * 0.5 # m3

# Guesses for the component vector inputs
T_guess = T_1
P_guess = P_1

# Contuct the state variable componennt vector
x0 = ComponentVector(
    VAP=(
        M=M_VAP,
        T=150.0,
        x=(
            O2=0.5,
            CO2=0.1,
            C2H4=0.3,
            C2H6=0.1,
            VAc=0.0,
            H2O=0.0,
            HAc=0.0
        )
    ),
    RCT=(
        T=(
            _1=150.0,
            _2=150.0,
            _3=150.0,
            _4=150.0,
            _5=150.0,
            _6=150.0,
            _7=150.0,
            _8=150.0,
            _9=150.0,
            _10=150.0),
        C=(
            O2=(
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5,
                _9=0.5,
                _10=0.5),
            CO2=(
                _1=0.1,
                _2=0.1,
                _3=0.1,
                _4=0.1,
                _5=0.1,
                _6=0.1,
                _7=0.1,
                _8=0.1,
                _9=0.1,
                _10=0.1),
            C2H4=(
                _1=0.3,
                _2=0.3,
                _3=0.3,
                _4=0.3,
                _5=0.3,
                _6=0.3,
                _7=0.3,
                _8=0.3,
                _9=0.3,
                _10=0.3),
            H2O=(
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5,
                _9=0.5,
                _10=0.5),
            VAc=(
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5,
                _9=0.5,
                _10=0.5),
            HAc=(
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5,
                _9=0.5,
                _10=0.5)
        )
    ),
    SEP=(
        M=M_SEP,
        P=120.0,
        T=(
            liq=100.0,
            vap=120.0),
        x=(
            O2=0.5,
            CO2=0.1,
            C2H4=0.3,
            VAc=0.1,
            H2O=0.0,
            HAc=0.0),
        y=(
            O2=0.5,
            CO2=0.1,
            C2H4=0.3,
            VAc=0.1,
            H2O=0.0,
            HAc=0.0)
    ),
    ABS=(
        M=(
            _b=0.5,
            _1=0.5,
            _2=0.5,
            _3=0.5,
            _4=0.5,
            _5=0.5,
            _6=0.5,
            _7=0.5,
            _8=0.5),
        T=(_b=100.0,
            _1=100.0,
            _2=100.0,
            _3=100.0,
            _4=100.0,
            _5=100.0,
            _6=100.0,
            _7=100.0,
            _8=100.0),
        C=(
            O2=(_b=0.5,
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5),
            CO2=(_b=0.1,
                _1=0.1,
                _2=0.1,
                _3=0.1,
                _4=0.1,
                _5=0.1,
                _6=0.1,
                _7=0.1,
                _8=0.1),
            C2H4=(_b=0.3,
                _1=0.3,
                _2=0.3,
                _3=0.3,
                _4=0.3,
                _5=0.3,
                _6=0.3,
                _7=0.3,
                _8=0.3),
            VAc=(_b=0.1,
                _1=0.1,
                _2=0.1,
                _3=0.1,
                _4=0.1,
                _5=0.1,
                _6=0.1,
                _7=0.1,
                _8=0.1),
            H2O=(_b=0.5,
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5),
            HAc=(_b=0.5,
                _1=0.5,
                _2=0.5,
                _3=0.5,
                _4=0.5,
                _5=0.5,
                _6=0.5,
                _7=0.5,
                _8=0.5)
        )
    ),
    TK=(
        M=M_TK,
        T=80.0,
        x=(VAc=0.1,
            HAc=0.0)
    ),
    COMP=(
        T=130.0,
        P=200.0
    ),
    T=(
        H1=150.0,
        H2=150.0,
        C3=100.0,
        C4=100.0,
        C5=100.0,
        S34=100.0,
        S9=100.0
    )
)


# Construct the parameter (input) component vector
p = ComponentVector(
    S1=(
        f=f_S1,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S2=(
        f=f_S2,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S3=(
        f=f_S3,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S4=(
        f=f_S4,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S5=(
        f=f_S5,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S6=(
        f=f_S6,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S7=(
        f=f_S7,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S8=(
        f=f_S8,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S9=(
        f=f_S9,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S10=(
        f=f_S10,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S11=( # Stream 11 is unique because S10 partial condenes before entering the separator
        f=f_S11,
        T=T_guess,
        P=P_guess,
        x=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        ),
        y=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S12=(
        f=f_S12,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S13=(
        f=f_S13,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S14=(
        f=f_S14,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S15=(
        f=f_S15,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S16=(
        f=f_S16,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S17=(
        f=f_S17,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S18=(
        f=f_S18,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S19=(
        f=f_S19,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S20=(
        f=f_S20,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S21=(
        f=f_S21,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S22=(
        f=f_S22,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S23=(
        f=f_S23,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S24=(
        f=f_S24,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S25=(
        f=f_S25,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S26=(
        f=f_S26,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S27=(
        f=f_S27,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S28=(
        f=f_S28,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S29=(
        f=f_S29,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S30=(
        f=f_S30,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S31=(
        f=f_S31,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S32=(
        f=f_S32,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S33=(
        f=f_S33,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
    S34=(
        f=f_S34,
        T=T_guess,
        P=P_guess,
        z=(
            C2H4=0.1,
            O2=0.5,
            CO2=0.1,
            H2O=0.1,
            VAc=0.1,
            HAc=0.1,
        )
    ),
)

# Wilson parameter and activity coefficient calculations
O2 = thermo["O2"]
CO2 = thermo["CO2"]
C2H4 = thermo["C2H4"]
VAc = thermo["VAc"]
H2O = thermo["H2O"]
HAc = thermo["HAc"]

# Dictionary of Wilson Parameters (aij) of Each Species Pair Existing in the Liquid Phase of the Vinyl Acetate Process
a_ij = Dict(
    ("VAc", "VAc") => 0.0,
    ("VAc", "H2O") => 2266.4,
    ("VAc", "HAc") => 726.7,
    ("H2O", "VAc") => 1384.6,
    ("H2O", "H2O") => 0.0,
    ("H2O", "HAc") => 230.6,
    ("HAc", "VAc") => -136.1,
    ("HAc", "H2O") => 670.7,
    ("HAc", "HAc") => 0.0,
    ("VAc", "O2") => 0.0,
    ("VAc", "CO2") => 0.0,
    ("VAc", "C2H4") => 0.0,
    ("H2O", "O2") => 0.0,
    ("H2O", "CO2") => 0.0,
    ("H2O", "C2H4") => 0.0,
    ("HAc", "O2") => 0.0,
    ("HAc", "CO2") => 0.0,
    ("HAc", "C2H4") => 0.0,
    ("O2", "C2H4") => 0.0,
    ("O2", "HAc") => 0.0,
    ("O2", "H2O") => 0.0,
    ("O2", "VAc") => 0.0,
    ("O2", "O2") => 0.0,
    ("O2", "CO2") => 0.0,
    ("CO2", "C2H4") => 0.0,
    ("CO2", "HAc") => 0.0,
    ("CO2", "H2O") => 0.0,
    ("CO2", "VAc") => 0.0,
    ("CO2", "O2") => 0.0,
    ("CO2", "CO2") => 0.0,
    ("C2H4", "O2") => 0.0,
    ("C2H4", "CO2") => 0.0,
    ("C2H4", "HAc") => 0.0,
    ("C2H4", "H2O") => 0.0,
    ("C2H4", "VAc") => 0.0,
    ("C2H4", "C2H4") => 0.0
)



R_gas = 1.987 # cal/(mol·K)
T_state = 200
x_C2H4 = 0.2
x_O2 = 0.3
x_HAc = 0.1
x_VAc = 0.15
x_CO2 = 0.15
x_H2O = 0.1

i = ["C2H4", "O2", "HAc", "VAc", "CO2", "H2O"]
x = [x_C2H4, x_O2, x_HAc, x_VAc, x_CO2, x_H2O]
j = filter(x -> x != "C2H4", i)
Lambda_C2H4_j = [thermo["C2H4"].v / thermo[j[index]].v * exp(-a_ij[("C2H4", j[index])] / (R_gas * (T_state + 273.15))) for index in 1:length(j)]

j = filter(x -> x != "O2", i)
Lambda_O2_j = [thermo["O2"].v / thermo[j[index]].v * exp(-a_ij[("O2", j[index])] / (R_gas * (T_state + 273.15))) for index in 1:length(j)]

j = filter(x -> x != "HAc", i)
Lambda_HAc_j = [thermo["HAc"].v / thermo[j[index]].v * exp(-a_ij[("HAc", j[index])] / (R_gas * (T_state + 273.15))) for index in 1:length(j)]

j = filter(x -> x != "VAc", i)
Lambda_VAc_j = [thermo["VAc"].v / thermo[j[index]].v * exp(-a_ij[("VAc", j[index])] / (R_gas * (T_state + 273.15))) for index in 1:length(j)]

j = filter(x -> x != "CO2", i)
Lambda_CO2_j = [thermo["CO2"].v / thermo[j[index]].v * exp(-a_ij[("CO2", j[index])] / (R_gas * (T_state + 273.15))) for index in 1:length(j)]

j = filter(x -> x != "H2O", i)
Lambda_H2O_j = [thermo["H2O"].v / thermo[j[index]].v * exp(-a_ij[("H2O", j[index])] / (R_gas * (T_state + 273.15))) for index in 1:length(j)]

Lambda_C2H4_j = vcat(0, Lambda_C2H4_j)
Lambda_O2_j = vcat(Lambda_O2_j[1], 0, Lambda_O2_j[2:end])
Lambda_HAc_j = vcat(Lambda_HAc_j[1:2], 0, Lambda_HAc_j[3:end])
Lambda_VAc_j = vcat(Lambda_VAc_j[1:3], 0, Lambda_VAc_j[4:end])
Lambda_CO2_j = vcat(Lambda_CO2_j[1:4], 0, Lambda_CO2_j[5:end])
Lambda_H2O_j = vcat(Lambda_H2O_j[1:5], 0)
# What I need and am having trouble writing is: ln gamma_i = 1 - ln(sum_j x_j * lambda_ji) - sum_j (x_j * lambda_ij / sum_k x_k * lambda_kj)
# Can you write this out for me? :)
ln_gamma_C2H4 = 1 - log(sum(x[j_index] * Lambda_C2H4_j[j_index] for j_index in 1:length(i))) - sum(x[j_index] * Lambda_C2H4_j[j_index] / sum(x[k_index] * Lambda_C2H4_j[k_index] for k_index in 1:length(i)) for j_index in 1:length(i))
ln_gamma_O2 = 1 - log(sum(x[j_index] * Lambda_O2_j[j_index] for j_index in 1:length(i))) - sum(x[j_index] * Lambda_O2_j[j_index] / sum(x[k_index] * Lambda_O2_j[k_index] for k_index in 1:length(i)) for j_index in 1:length(i))
ln_gamma_HAc = 1 - log(sum(x[j_index] * Lambda_HAc_j[j_index] for j_index in 1:length(i))) - sum(x[j_index] * Lambda_HAc_j[j_index] / sum(x[k_index] * Lambda_HAc_j[k_index] for k_index in 1:length(i)) for j_index in 1:length(i))
ln_gamma_VAc = 1 - log(sum(x[j_index] * Lambda_VAc_j[j_index] for j_index in 1:length(i))) - sum(x[j_index] * Lambda_VAc_j[j_index] / sum(x[k_index] * Lambda_VAc_j[k_index] for k_index in 1:length(i)) for j_index in 1:length(i))
ln_gamma_CO2 = 1 - log(sum(x[j_index] * Lambda_CO2_j[j_index] for j_index in 1:length(i))) - sum(x[j_index] * Lambda_CO2_j[j_index] / sum(x[k_index] * Lambda_CO2_j[k_index] for k_index in 1:length(i)) for j_index in 1:length(i))
ln_gamma_H2O = 1 - log(sum(x[j_index] * Lambda_H2O_j[j_index] for j_index in 1:length(i))) - sum(x[j_index] * Lambda_H2O_j[j_index] / sum(x[k_index] * Lambda_H2O_j[k_index] for k_index in 1:length(i)) for j_index in 1:length(i))

gamma_C2H4 = exp(ln_gamma_C2H4)
gamma_O2 = exp(ln_gamma_O2)
gamma_HAc = exp(ln_gamma_HAc)
gamma_VAc = exp(ln_gamma_VAc)
gamma_CO2 = exp(ln_gamma_CO2)
gamma_H2O = exp(ln_gamma_H2O)

println("Activity Coefficients at T = $T_state °C:")
println("γ_C2H4 = $gamma_C2H4")
println("γ_O2 = $gamma_O2")
println("γ_HAc = $gamma_HAc")
println("γ_VAc = $gamma_VAc")
println("γ_CO2 = $gamma_CO2")
println("γ_H2O = $gamma_H2O")
