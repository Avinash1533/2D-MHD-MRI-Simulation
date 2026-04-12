import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. RESEARCH-GRADE SETUP
# ==========================================
print("Initializing Fast 2D MRI Simulation Engine...")
NR, NPHI = 128, 128
R_MIN, R_MAX = 2.0, 12.0
dr = (R_MAX - R_MIN) / NR
dphi = (2 * np.pi) / NPHI

# Grids
R = np.linspace(R_MIN, R_MAX, NR).reshape(NR, 1)
PHI = np.linspace(0, 2 * np.pi, NPHI).reshape(1, NPHI)
R_grid = R * np.ones((NR, NPHI))
PHI_grid = np.ones((NR, NPHI)) * PHI

# FIX: Added a specific array for the radial edges to prevent shape crashing
R_edge = np.linspace(R_MIN, R_MAX, NR + 1).reshape(NR + 1, 1)

# Fluid Base State
rho = np.exp(-((R_grid - 6.0) ** 2) / 2.0) + 0.05
v_phi = np.sqrt(1.0 / R_grid)
v_r = np.zeros((NR, NPHI))

# Perturbation (The spark to trigger MRI)
np.random.seed(42)
v_r += np.random.uniform(-0.05, 0.05, (NR, NPHI)) * (rho > 0.1)
v_phi += np.random.uniform(-0.05, 0.05, (NR, NPHI)) * (rho > 0.1)

# Magnetic Fields
B_r = np.zeros((NR + 1, NPHI))

# FIX: We must seed an in-plane "Toroidal" magnetic field!
# This wraps the magnetic field around the disk.
B_phi = 0.05 * np.ones((NR, NPHI + 1))

# We zero out B_z because it doesn't interact in our flat 2D plane
B_z = np.zeros((NR, NPHI))

# Data arrays for paper graphs
time_data, mass_accretion_rate, magnetic_energy = [], [], []

# ==========================================
# 2. THE COUPLED PHYSICS LOOP
# ==========================================
t = 0.0
t_final = 100.0
dt = 0.002

print("Engine running. Solving vectorized equations...")

# The CFL Safety Factor
C_cfl = 0.3

while t < t_final:
    # --- 0. ADAPTIVE TIME-STEPPING (Strict CFL) ---
    B_r_c = 0.5 * (B_r[:-1, :] + B_r[1:, :])
    B_phi_c = 0.5 * (B_phi[:, :-1] + B_phi[:, 1:])
    B_mag = np.sqrt(B_r_c ** 2 + B_phi_c ** 2 + B_z ** 2)

    # True Alfvén speed
    v_alfven = B_mag / np.sqrt(rho)

    max_v_r = np.max(np.abs(v_r) + v_alfven)
    max_v_phi = np.max(np.abs(v_phi) + v_alfven)

    dt_r = dr / max_v_r
    dt_phi = (R_MIN * dphi) / max_v_phi

    # We MUST respect the strict mathematical limit. No artificial floors!
    dt = C_cfl * min(dt_r, dt_phi)

    # --- A. MAGNETIC UPDATE (Constrained Transport) ---
    E_z_c = (v_r * B_phi_c) - (v_phi * B_r_c)

    E_z = np.zeros((NR + 1, NPHI + 1))
    E_z[1:-1, 1:-1] = 0.25 * (E_z_c[1:, 1:] + E_z_c[:-1, 1:] +
                              E_z_c[1:, :-1] + E_z_c[:-1, :-1])
    E_z[:, 0] = E_z[:, -2]
    E_z[:, -1] = E_z[:, 1]

    B_r[1:-1, :] -= dt * (E_z[1:-1, 1:] - E_z[1:-1, :-1]) / (R_edge[1:-1] * dphi)
    B_phi[:, 1:-1] += dt * (E_z[1:, 1:-1] - E_z[:-1, 1:-1]) / dr

    B_phi[:, 0] = B_phi[:, -2]
    B_phi[:, -1] = B_phi[:, 1]

    # --- B. LORENTZ FORCE & FLUID DYNAMICS ---
    B_r_c = 0.5 * (B_r[:-1, :] + B_r[1:, :])
    B_phi_c = 0.5 * (B_phi[:, :-1] + B_phi[:, 1:])

    dB_phi_dr = np.gradient(B_phi_c, dr, axis=0)
    dB_r_dphi = (np.roll(B_r_c, -1, axis=1) - np.roll(B_r_c, 1, axis=1)) / (2.0 * dphi)
    J_z = dB_phi_dr - (dB_r_dphi / R_grid)

    F_mag_r = J_z * B_phi_c
    F_mag_phi = -J_z * B_r_c

    gravity_centrifugal = (v_phi ** 2 / R_grid) - (1.0 / R_grid ** 2)

    v_r += dt * (gravity_centrifugal + (F_mag_r / rho))
    v_phi += dt * (F_mag_phi / rho)

    # THE CRITICAL FIX: The Vacuum Velocity Cap
    # Prevent runaway acceleration in the low-density regions
    v_r = np.clip(v_r, -10.0, 10.0)
    v_phi = np.clip(v_phi, -10.0, 10.0)

    # --- C. DATA TRACKING ---
    if int(t / dt) % 50 == 0:
        time_data.append(t)
        accretion = np.sum(rho[1, :] * np.abs(np.minimum(0, v_r[1, :])))
        mass_accretion_rate.append(accretion)
        mag_energy = np.sum(B_r_c ** 2 + B_phi_c ** 2 + B_z ** 2)
        magnetic_energy.append(mag_energy)

        if int(t / dt) % 500 == 0:
            print(f"Progress: {t:.2f} / {t_final} seconds | Current dt: {dt:.5f}")

    t += dt

# ==========================================
# 3. RESULTS & VISUALIZATION
# ==========================================
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
plt.style.use('dark_background')

# Graph 1: The Final State
X = R_grid * np.cos(PHI_grid)
Y = R_grid * np.sin(PHI_grid)
c1 = axes[0].pcolormesh(X, Y, v_r, cmap='coolwarm', shading='nearest')
axes[0].set_title("Radial Velocity (Turbulent Flow)")
axes[0].axis('off')
fig.colorbar(c1, ax=axes[0], label="v_r")

# Graph 2: Accretion Proof
axes[1].plot(time_data, mass_accretion_rate, color='#fde725', lw=2)
axes[1].set_title("Mass Accretion Rate over Time")
axes[1].set_xlabel("Time (s)")
axes[1].set_ylabel("Inward Mass Flux")
axes[1].grid(alpha=0.3)

# Graph 3: Magnetic Amplification
axes[2].plot(time_data, magnetic_energy, color='#5ec962', lw=2)
axes[2].set_title("Total Magnetic Energy Growth")
axes[2].set_xlabel("Time (s)")
axes[2].set_ylabel("Magnetic Energy (B^2)")
axes[2].grid(alpha=0.3)

plt.tight_layout()
# ==========================================
# 4. GROWTH RATE (LAMBDA) CALCULATOR
# ==========================================
# 1. Convert your lists to NumPy arrays so the math works!
t_data_np = np.array(time_data)
emag_data_np = np.array(magnetic_energy)

# 2. Set the growth window (Adjust these if your new run spikes earlier/later)
t_start = 20
t_end = 80

# 3. Filter the data for just the "Hockey Stick" part
mask = (t_data_np >= t_start) & (t_data_np <= t_end)
t_growth = t_data_np[mask]
emag_growth = emag_data_np[mask]

# 4. The Math: Fit the log-linear line
log_emag = np.log(emag_growth)
coefficients = np.polyfit(t_growth, log_emag, 1)
lambda_val = coefficients[0]
intercept = coefficients[1]

# 5. Print the result to the PyCharm Terminal
print("\n" + "="*50)
print(f"🔥 CALCULATED GROWTH RATE (Lambda): {lambda_val:.4f}")
print("="*50 + "\n")

# 6. Add the red dashed line to your 3rd graph (axes[2])
fit_line = np.exp(intercept) * np.exp(lambda_val * t_growth)
axes[2].plot(t_growth, fit_line, label=f"Fit ($\lambda$ = {lambda_val:.4f})",
             color='red', linestyle='--', linewidth=3)
axes[2].legend()

# 7. NOW show the final dashboard
plt.tight_layout()
plt.show()