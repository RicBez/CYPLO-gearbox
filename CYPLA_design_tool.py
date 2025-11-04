#!/usr/bin/env python3
"""
Planetary-stage design assistant (interactive console tool)

This tool guides the user through the design of the planetary stage of a reducer
(input = sun, output = carrier, fixed ring), given mechanical constraints.

Assumptions and abbreviations:
 - m : module [mm]
 - z_s : number of teeth of the sun
 - z_r : number of teeth of the internal ring (internal gearing)
 - z_p : number of teeth of each planet (computed as (z_r - z_s)/2)
 - D = m * z  : pitch diameter [mm] (used as diameter metric)
 - Outer ring outer diameter estimated as D_ring_outer = m*(z_r + 2)
 - carrier_radius (planet center distance from system center) = m*(z_s + z_p)/2
 - reduction ratio i = 1 + z_r/z_s  (sun -> carrier for sun-driven / ring-fixed planetary)
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ======================================
# ---- UTILITY FUNCTIONS FOR DRAWING ----
# ======================================

def circle(x, y, r, n=200):
    """Return the coordinates of a circle centered at (x, y) with radius r."""
    theta = np.linspace(0, 2 * np.pi, n)
    return x + r * np.cos(theta), y + r * np.sin(theta)

def cycloidal_disk(ecc, roll_r, rollers_num, pitch_R, res=2000):
    """
    Generate the profile of a compact-cam cycloidal disk
    based on geometric relations between pitch circle, rollers, and eccentricity.
    """
    points = []
    theta = np.linspace(0, 2 * np.pi, res)
    for i in theta:
        delta = np.arctan(
            np.sin((1 - rollers_num) * i) /
            ((pitch_R / (ecc * rollers_num)) - np.cos((1 - rollers_num) * i))
        )
        x = (pitch_R * np.cos(i)) \
            - (roll_r * np.cos(i + delta)) \
            - (ecc * np.cos(rollers_num * i))
        y = (-pitch_R * np.sin(i)) \
            + (roll_r * np.sin(i + delta)) \
            + (ecc * np.sin(rollers_num * i))
        points.append((x, y))
    points.append(points[0])
    return np.array(points)

# ======================================
# ---- PLANETARY STAGE DESIGN ----
# ======================================

def planetary_design():
    print("\n=== PLANETARY STAGE DESIGN ===\n")

    # Step 1: maximum radial size
    R_max = float(input(
        "Insert the maximum radius allowed for the overall reducer [mm]: "
    ))

    # Step 2: transmission ratio
    i_plan = float(input(
        "Insert the required transmission ratio for the planetary stage "
        "(allowed range: 1:3 to 1:9). This corresponds to input on the sun, output on the carrier: "
    ))
    if not (3 <= i_plan <= 9):
        print("\nERROR: The planetary ratio must be between 1:3 and 1:9. "
              "Please restart and provide a valid value.")
        return None

    # Step 3: gear module
    m = float(input(
        "Insert the gear module [mm]. This directly affects gear tooth size "
        "and consequently the minimum feasible number of teeth: "
    ))

    # Step 4: compute max number of ring teeth
    Z_ring_user = math.floor(2 * (R_max - 4.5) / m) # foot diameter (2.5 * m) + margin (2*m)
    print(f"\n--> Maximum feasible number of teeth for the ring gear "
          f"(internal gear): {Z_ring_user}")

    # Step 5: approval
    ok = input("Do you want to proceed with this calculated number of teeth? (y/n): ")
    if ok.lower() != "y":
        choice = input(
            "Would you like to restart from step (2/3/4), "
            "or directly specify the ring gear tooth count (type 'ring')? "
        )
        if choice.isdigit():
            return planetary_design()
        elif choice.lower() == "ring":
            Z_ring_user = int(input("Insert the desired number of teeth for the ring gear: "))
            ring_size = (Z_ring_user * m / 2) + (2.5 * m) + 2*m # foot diameter (2.5 * m) + margin (2*m)
            print(f"Resulting engumbrance radius of the ring (thus, of the reducer): {ring_size:.2f} mm")
            #return {"Z_ring": Z_ring_user, "m": m, "i_plan": i_plan}
        else:
            return planetary_design()

    # Step 6: search for gear combinations
    Z_ring = Z_ring_user
    solutions = []
    for Z_sun in range(8, Z_ring - 8, 1):
        Z_carrier = (Z_ring + Z_sun) / 2
        i = 1 + Z_ring / Z_sun
        if abs(i - i_plan) < 0.2:
            r_sun = Z_sun * m / 2
            r_carrier = Z_carrier * m / 2
            # Derived number of teeth per planet
            Z_planet = (Z_ring - Z_sun) / 2
            # Compute all possible planet counts
            possible_N = []
            diff = Z_ring - Z_sun
            for Np in range(3, 8):  # typically between 3 and 7 planets
                if diff % Np == 0:
                    possible_N.append(Np)
            solutions.append((Z_sun, int(Z_ring), int(Z_carrier), i, r_sun, r_carrier, possible_N))

    if not solutions:
        print("\nNo feasible tooth combinations were found for the requested ratio. The solutions without spcified number of planets are not available."
              "Consider adjusting the ratio or module.")
        return None

    print("\nPossible gear combinations (Z_sun, Z_ring, Z_carrier, ratio, r_sun, r_carrier, possible number of equally distanced planets):")
    for idx, sol in enumerate(solutions):
        print(f"{idx+1}: {sol}")

    choice = int(input("\nSelect the index of the preferred solution: ")) - 1
    Z_sun, Z_ring, Z_carrier, i, r_sun, r_carrier, Np = solutions[choice]

    choice = int(input("\nSelect the index of the desired nuber of planets from the vector above: ")) - 1
    num_pla = Np[choice]

    print(f"\n--> Selected planetary stage:\n"
          f"Sun teeth = {Z_sun}, Ring teeth = {Z_ring}, Carrier = {Z_carrier}\n"
          f"Ratio = {i:.2f}, Sun radius = {r_sun:.2f} mm, Carrier radius = {r_carrier:.2f} mm, Number of planets = {num_pla}")

    return {"Z_sun": Z_sun, "Z_ring": Z_ring, "Z_carrier": Z_carrier,
            "m": m, "i": i, "r_sun": r_sun, "r_carrier": r_carrier, "N_planets": num_pla}

# ======================================
# ---- CYCLOIDAL COMPACT-CAM DESIGN ----
# ======================================

def compact_cam_design(D_sun, m):
    print("\n=== COMPACT-CAM CYCLOIDAL DESIGN ===\n")
    
    # --- INPUT SECTION ---
    Dr = float(input("Insert the roller diameter Dr [mm]: "))
    u = float(input("Desired transmission ratio (u): "))
    uth = float(input("Tolerance on the transmission ratio (uth): "))
    
    # Calculation of maximum number of rollers (the approximation has a margin, but can be adjusted)
    def custom_round(x):
        return math.floor(x) if (x - int(x)) < 0.5 else math.ceil(x)
    def compute_maxN(D, Dr):
        return custom_round((D - Dr) * math.pi / Dr)
    def compute_Dmin(N, Dr):
        return N*Dr/math.pi
    
    # Compute and display suggested maximum diameter
    Dme_calc = D_sun - (2.5 * m) - 2 * m - Dr
    print(f"\nSuggested maximum available pitch diameter for the rollers: {Dme_calc:.2f} mm")
    Dme = float(input("Enter desired maximum pitch diameter (must be ≤ suggested value): "))
    NRsup = compute_maxN(Dme, Dr)
    NRinf = 2
    print(f"Maximum number of rollers: {NRsup}")
    # --- DESIGN PROCESS LOOP (restarts from beginning if parameters are changed) ---
    while True:
        # --- DESIGN LOOP (allows re-entry after plotting) ---
        while True:
            # Search for valid combinations for the two stages
            solutions = []
            Nvec = np.arange(NRinf, NRsup + 1)
            for N1 in Nvec:
                for N2 in Nvec:
                    n1 = N1 - 1
                    n2 = N2 - 1
                    sh = n1 * N2
                    if (sh - n2 * N1) == 0:
                        continue
                    ut = sh / (sh - n2 * N1)
                    if abs(ut) <= u + uth and abs(ut) >= u - uth:
                        solutions.append([ut, N1, N2])
            if not solutions:
                print("No solutions found with the provided parameters.")
                return None
            
            # User selection of the solution
            for idx, sol in enumerate(solutions):
                print(f"{idx+1}: Ratio = {sol[0]:.4f}, N1={sol[1]}, N2={sol[2]}")
            choice = int(input("Select the desired solution: ")) - 1
            ut, N1t, N2t = solutions[choice]
            print(f"Selected: N1={N1t}, N2={N2t}, Ratio={ut:.4f}")
            D1min = compute_Dmin(N1t, Dr)
            D2min = compute_Dmin(N2t, Dr)
            
            # --- Pitch circle diameters ---
            while True:
                D1 = float(input("Pitch circle diameter stage 1 D1 [mm]: "))
                D2 = float(input("Pitch circle diameter stage 2 D2 [mm]: "))
                if D1 > Dme or D2 > Dme:
                    print(f"\n Invalid input:")
                    print(f"   Both D1 and D2 must be ≤ {Dme:.2f} mm to respect the geometric constraint.")
                    print("    Please enter valid diameters.\n")
                elif D1 < D1min or D2 < D2min:
                    print(f"\n Invalid input:")
                    print(f"    D1 must be > {D1min:.2f} and D2 must be > {D2min:.2f} mm to respect the geometric constraint.")
                    print("    Please enter valid diameters.\n")
                else:
                    print(" Diameters accepted.\n")
                    break
            
            ecc = float(input("Desired eccentricity [mm]: "))
            
            pitch_R1 = D1 / 2
            pitch_R2 = D2 / 2
            roll_r = Dr / 2
            
            # Function to generate compact-cam profile (two stages)
            def cycloidal_disk_compact(ecc, roll_r, N, pitch_R, res=5000):
                points = []
                theta = np.linspace(0, 2*np.pi, res)
                for i in theta:
                    delta = np.arctan(np.sin((1 - N)*i) / ((pitch_R / (ecc*N)) - np.cos((1 - N)*i)))
                    x = pitch_R * np.cos(i) - roll_r * np.cos(i + delta) - ecc * np.cos(N*i)
                    y = -pitch_R * np.sin(i) + roll_r * np.sin(i + delta) + ecc * np.sin(N*i)
                    points.append((x, y))
                points.append(points[0])
                return np.array(points)
            
            # Profile generation for both stages
            cycloid1 = cycloidal_disk_compact(ecc, roll_r, N1t, pitch_R1)
            cycloid2 = cycloidal_disk_compact(ecc, roll_r, N2t, pitch_R2)
            
            # Plotting of the two stages
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,7))
            ax1.set_aspect('equal'); ax1.set_title("Stage 1 Compact-Cam cycloidal profile")
            ax1.plot(cycloid1[:,0], cycloid1[:,1], color='red')
            for i in range(N1t):
                angle = 2*np.pi*i/N1t
                x = pitch_R1 * np.cos(angle)
                y = pitch_R1 * np.sin(angle)
                xr, yr = circle(x - ecc, y, roll_r)
                ax1.fill(xr, yr, 'yellow')
            
            ax2.set_aspect('equal'); ax2.set_title("Stage 2 Compact-Cam cycloidal profile")
            ax2.plot(cycloid2[:,0], cycloid2[:,1], color='blue')
            for i in range(N2t):
                angle = 2*np.pi*i/N2t
                x = pitch_R2 * np.cos(angle)
                y = pitch_R2 * np.sin(angle)
                xr, yr = circle(x - ecc, y, roll_r)
                ax2.fill(xr, yr, 'yellow')
            
            plt.show()
            
            # --- ASK USER IF HE WANTS TO MODIFY PARAMETERS ---
            print("\nWould you like to modify any design parameters and re-run the simulation?")
            print("Options:")
            print(" 1 - Change the design parameters")
            print(" 0 - Finish and save results")
            
            choice2 = input("Enter your choice: ")
            
            if choice2 == "1":
                print("\nRestarting the design process with new parameters...\n")
                continue  # Torna all'inizio del ciclo esterno per rifare tutto
            
            elif choice2 == "0":
                print("\nDesign process completed.\n")
                break
            
            else:
                print("Invalid choice, exiting.")
                break
        
        # --- RETURN FINAL PARAMETERS ---
        return {"Ratio": ut, "N1": N1t, "D1": D1, "N2": N2t, "D2": D2, "Dr": Dr, "ecc": ecc}





# ======================================
# ---- INTEGRATION PLANETARY + CYCLOIDAL ----
# ======================================


def integrate_design(planet_data, cyclo_data):
    print("\n=== INTEGRATION OF PLANETARY AND COMPACT-CAM STAGES ===\n")

    # === Retrieve planetary data ===
    Zs = planet_data["Z_sun"]
    Zr = planet_data["Z_ring"]
    m = planet_data["m"]
    r_carrier = planet_data["r_carrier"]
    Np = planet_data["N_planets"]
    upl = planet_data["i"]

    # Derived planet gear parameters
    Zp = (Zr - Zs) / 2
    r_sun = (Zs * m) / 2
    r_ring = (Zr * m) / 2
    r_planet = (Zp * m) / 2

    # === Retrieve cycloidal data ===
    ucy = cyclo_data["Ratio"]
    D1 = cyclo_data["D1"]
    D2 = cyclo_data["D2"]
    Dr = cyclo_data["Dr"]
    ecc = cyclo_data["ecc"]
    N1t = cyclo_data["N1"]
    N2t = cyclo_data["N2"]

    pitch_R1 = D1 / 2
    pitch_R2 = D2 / 2
    roll_r = Dr / 2

    utot = ucy*upl;

    # === Helper function ===
    def cycloidal_disk_compact(ecc, roll_r, N, pitch_R, res=5000):
        points = []
        theta = np.linspace(0, 2*np.pi, res)
        for i in theta:
            delta = np.arctan(np.sin((1 - N)*i) / ((pitch_R / (ecc*N)) - np.cos((1 - N)*i)))
            x = pitch_R * np.cos(i) - roll_r * np.cos(i + delta) - ecc * np.cos(N*i)
            y = -pitch_R * np.sin(i) + roll_r * np.sin(i + delta) + ecc * np.sin(N*i)
            points.append((x, y))
        points.append(points[0])
        return np.array(points)

    # === Generate the two cycloidal profiles ===
    cycloid1 = cycloidal_disk_compact(ecc, roll_r, N1t, pitch_R1)
    cycloid2 = cycloidal_disk_compact(ecc, roll_r, N2t, pitch_R2)

    # === Function to plot the planetary geometry (without legenda) ===
    def plot_planetary(ax):
        ax.set_aspect("equal")

        # Ring gear
        xr, yr = circle(0, 0, r_ring)
        ax.plot(xr, yr, "k", linewidth=2)

        # Sun gear
        xs, ys = circle(0, 0, r_sun)
        ax.plot(xs, ys, "b", linewidth=2)

        # Carrier circle
        xc, yc = circle(0, 0, r_carrier)
        ax.plot(xc, yc, "orange", linewidth=2)

        # Planet gears (equally spaced)
        for i in range(Np):
            angle = 2 * np.pi * i / Np
            xp_center = r_carrier * np.cos(angle)
            yp_center = r_carrier * np.sin(angle)
            xp, yp = circle(xp_center, yp_center, r_planet)
            ax.plot(xp, yp, "g", linewidth=1.5)
            ax.fill(xp, yp, "lightgreen", alpha=0.6)

    # === Create figure with two subplots ===
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    # ======== Stage 1 ========
    ax1.set_title("Integration: Planetary + Compact‑Cam Stage 1", fontsize=12)
    plot_planetary(ax1)
    ax1.plot(cycloid1[:, 0] + ecc, cycloid1[:, 1], "r", linewidth=1.8, label="Compact‑Cam Stage 1 Profile")

    xp1, yp1 = circle(0, 0, pitch_R1)
    ax1.plot(xp1, yp1, "r--", linewidth=1, label="Stage 1 Pitch Circle D1")

    for i in range(N1t):
        angle = 2 * np.pi * i / N1t
        x = pitch_R1 * np.cos(angle)
        y = pitch_R1 * np.sin(angle)
        xr, yr = circle(x, y, roll_r)
        ax1.fill(xr, yr, "yellow", alpha=0.6)

    # ======== Stage 2 ========
    ax2.set_title("Integration: Planetary + Compact‑Cam Stage 2", fontsize=12)
    plot_planetary(ax2)
    ax2.plot(cycloid2[:, 0] + ecc, cycloid2[:, 1], "r", linewidth=1.8, label="Compact‑Cam Stage 2 Profile")

    xp2, yp2 = circle(0, 0, pitch_R2)
    ax2.plot(xp2, yp2, "r--", linewidth=1, label="Stage 2 Pitch Circle D2")

    for i in range(N2t):
        angle = 2 * np.pi * i / N2t
        x = pitch_R2 * np.cos(angle)
        y = pitch_R2 * np.sin(angle)
        xr, yr = circle(x, y, roll_r)
        ax2.fill(xr, yr, "yellow", alpha=0.6)

    # === Custom Legends for both axes ===
    legend_elements = [
        Patch(facecolor='none', edgecolor='black', label='Ring Gear (primitive circle)'),
        Patch(facecolor='none', edgecolor='blue', label='Sun Gear (primitive circle)'),
        Patch(facecolor='none', edgecolor='orange', label='Carrier Radius'),
        Patch(facecolor='lightgreen', edgecolor='green', label='Planet Gears (primitive circle)'),
        Line2D([], [], color='r', linewidth=1.8, label='Compact‑Cam Profile'),
        Line2D([], [], color='r', linestyle='--', linewidth=1, label='Pitch Circle'),
        Patch(facecolor='yellow', edgecolor='goldenrod', label='Rollers')
    ]

    ax1.legend(handles=legend_elements, loc="upper right", fontsize=8)
    ax2.legend(handles=legend_elements, loc="upper right", fontsize=8)

    plt.show()

    #Zs = planet_data["Z_sun"]
    #Zr = planet_data["Z_ring"]
    #m = planet_data["m"]
    #r_carrier = planet_data["r_carrier"]
    #Np = planet_data["N_planets"]
    #Zp = (Zr - Zs) / 2

    print(f"\n--> Planetary parameters:\n"
          f"Planetary ratio = {upl}, Sun teeth = {Zs}, Ring teeth = {Zr}, Planet teeth = {Zp}, Radius Carrier = {r_carrier}\n"
          f"Module = {m:.2f}, Number of planets = {Np}")
    
    # === Retrieve cycloidal data ===
    #D1 = cyclo_data["D1"]
    #D2 = cyclo_data["D2"]
    #Dr = cyclo_data["Dr"]
    #ecc = cyclo_data["ecc"]
    #N1t = cyclo_data["N1"]
    #N2t = cyclo_data["N2"]

    #pitch_R1 = D1 / 2
    #pitch_R2 = D2 / 2
    #roll_r = Dr / 2

    print(f"\n--> Cycloidal compact-cam parameters:\n"
          f"Cycloidal Ratio = {ucy:.2f}, Rollers Diameter = {Dr}, Eccentricity = {ecc}, Stage 1 Number of rollers = {N1t}, Stage 1 Pitch Diameter = {D1}\n"
          f"Stage 2 Number of rollers = {N2t}, Stage 2 Pitch Diameter = {D2}")

    print(f"\n--> Total Reduction ratio:\n"
          f"-> {utot:.2f}")

    print("--> Integration completed. Two combined configurations (Stage 1 and Stage 2) have been visualized.\n")




# ======================================
# ---- MAIN ----
# ======================================

if __name__ == "__main__":
    planet = planetary_design()
    if planet:
        ok = input("\nDo you want to proceed with the compact-cam cycloidal stage inside the sun? (y/n): ")
        if ok.lower() == "y":
            cyclo = compact_cam_design(D_sun=planet["r_sun"] * 2, m = planet["m"])
            if cyclo:
                integrate_design(planet, cyclo)

