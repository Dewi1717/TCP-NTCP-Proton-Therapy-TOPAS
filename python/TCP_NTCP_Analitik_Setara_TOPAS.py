# -*- coding: utf-8 -*-
"""
Pemodelan analitik TCP-NTCP yang parameternya DISAMAKAN dengan simulasi TOPAS.

Tujuan: agar perbandingan analitik vs TOPAS (tujuan penelitian ke-3) sah.
Satu-satunya perbedaan yang tersisa antara skrip ini dan TOPAS adalah
DISTRIBUSI DOSIS:
  - di sini  : dosis homogen (seluruh volume organ menerima dosis yang sama)
  - di TOPAS : DVH hasil Monte Carlo (dosis heterogen, bergantung geometri)
Sehingga selisih kurva dapat diartikan sebagai efek heterogenitas dosis dan
geometri organ, bukan karena input yang berbeda.

Rumus diambil langsung dari kode sumber TOPAS di folder ini:
  - Koreksi fraksinasi : TsVOutcomeModel::CalculateBED + TsOutcomeModelList
  - TCP Poisson        : TsPoissonModel.cc        (Warkentin et al. 2004)
  - NTCP LKB           : TsLKBModel.cc            (gEUD + probit)
  - NTCP Critical Elem.: TsCriticalElementModel.cc (Niemierko & Goitein 1991)
  - NTCP Critical Vol. : TsCriticalVolumeModel.cc (Stavrev et al. 2001)

Hanya membutuhkan numpy dan matplotlib (tidak perlu scipy).
"""

import csv
import glob
import math
import os

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. PARAMETER (samakan dengan file parameter TOPAS)
# =============================================================================

# Rentang dosis total (Gy) — samakan dengan rentang ScaleFactor di TOPAS
DOSE_MIN, DOSE_MAX, N_POINTS = 0.0, 150.0, 1501

# Fraksinasi (Sc/<scorer>/NumberOfFractions, DosePerFraction, AlphaOverBeta)
# CATATAN: laporan menulis DosePerFraction = 3.0 Gy, file
# "Parameter Simulasi_TOPAS" menulis 2.0 Gy. Isi sesuai file .txt TOPAS
# yang benar-benar dijalankan. Jika tidak di-set, default TOPAS = 2.0 Gy.
N_FRACTIONS = 20
D_REF = 2.0            # Gy, DosePerFraction referensi (default TOPAS)

# Target (TCP) — Sc/Target/...
# Target tidak men-set AlphaOverBeta, sehingga TOPAS memakai default 3 Gy.
# (Untuk tumor biasanya 10 Gy; ubah di sini DAN di TOPAS bila ingin 10.)
TARGET_AB = 3.0
TCD50 = 64.38          # Sc/Target/Poisson/TCD50
GAMMA50_TCP = 0.69     # Sc/Target/Poisson/Gamma50

# OAR1 = ginjal — Sc/OAR1/...
OAR1_AB = 3.0          # Sc/OAR1/AlphaOverBeta
# LKB preset "Kidney" (Burman et al. 1991, ginjal utuh): n, m, TD50
# Cek nilai di LKMModelBurman.txt -> Sc/LKBValues/Kidney
LKB_KIDNEY = dict(n=0.70, m=0.10, TD50=28.0)
# Critical Element (Sc/OAR1/CriticalElement/...)
CE_KIDNEY = dict(gamma=2.0, TD50=28.0)
# Critical Volume (Sc/OAR1/CriticalVolume/...)
CV_KIDNEY = dict(mu_cr=0.22, sigma_mu_cr=0.05, D50FSU=59.0, gamma50FSU=0.4)

# OAR2 = kolon — Sc/OAR2/... (preset "colon", Burman et al. 1991)
OAR2_AB = 3.0
LKB_COLON = dict(n=0.17, m=0.11, TD50=55.0)

# Rasio dosis OAR terhadap dosis target pada model homogen.
# 1.0 = seluruh OAR menerima dosis target (setara kondisi on-axis terburuk).
OAR_DOSE_RATIO = 1.0

# Kriteria therapeutic window (dihitung, bukan digambar manual)
TCP_MIN = 0.50         # TCP minimal
NTCP_MAX = 0.05        # NTCP maksimal yang dapat diterima

# Pembanding: model mekanistik LQ-Poisson dari kode awal (coding Nazwa.txt),
# kini DIFRAKSINASI agar setara (sebelumnya dosis total dianggap 1 fraksi).
SHOW_LQ_REFERENCE = True
LQ_ALPHA, LQ_BETA, LQ_N0 = 0.3, 0.03, 1e7

# Overlay hasil TOPAS (opsional). Letakkan file CSV di folder "topas_csv/":
#   kolom pertama = dosis total (Gy), kolom berikutnya = TCP / NTCP (0-1 atau %)
#   contoh header: Dose,TCP,NTCP_LKB,NTCP_CE,NTCP_CV
TOPAS_CSV_DIR = "topas_csv"

OUTPUT_PNG = "TCP_NTCP_Analitik_Setara_TOPAS.png"
OUTPUT_CSV = "TCP_NTCP_Analitik_Setara_TOPAS.csv"


# =============================================================================
# 2. FUNGSI — replika langsung kode TOPAS (bekerja pada DVH diferensial)
# =============================================================================

def probit(t):
    """TsVOutcomeModel::CalculateProbit."""
    return 0.5 * (1.0 + math.erf(t / math.sqrt(2.0)))


def eqd_ref(dose, ab, n_fx, d_ref):
    """TsVOutcomeModel::CalculateBED (sebenarnya EQD pada d_ref Gy/fraksi).
    dose = dosis total voxel; dosis per fraksi = dose / n_fx."""
    dose = np.asarray(dose, dtype=float)
    return dose * (ab + dose / n_fx) / (ab + d_ref)


def tcp_poisson(dose, vol, TCD50, gamma50):
    """TsPoissonModel::Initialize."""
    v = np.asarray(vol, float) / np.sum(vol)
    arg = np.sum(v * np.exp(2.0 * gamma50 * (1.0 - np.asarray(dose) / TCD50) / math.log(2.0)))
    return 0.5 ** arg


def ntcp_lkb(dose, vol, n, m, TD50):
    """TsLKBModel::Initialize (gEUD dengan eksponen 1/n)."""
    v = np.asarray(vol, float) / np.sum(vol)
    a = 1.0 / n
    geud = np.sum(v * np.asarray(dose, float) ** a) ** (1.0 / a)
    return probit((geud - TD50) / (m * TD50))


def ntcp_critical_element(dose, vol, gamma, TD50):
    """TsCriticalElementModel::Initialize."""
    v = np.asarray(vol, float) / np.sum(vol)
    p1 = 2.0 ** (-np.exp(math.e * gamma * (1.0 - np.asarray(dose) / TD50)))
    return 1.0 - np.prod((1.0 - p1) ** v)


def ntcp_critical_volume(dose, vol, mu_cr, sigma_mu_cr, D50FSU, gamma50FSU):
    """TsCriticalVolumeModel::Initialize."""
    v = np.asarray(vol, float) / np.sum(vol)
    sigma = sigma_mu_cr / -(mu_cr * math.log(mu_cr))
    dose = np.maximum(np.asarray(dose, float), 1e-12)   # hindari log(0)
    args = math.sqrt(2 * math.pi) * gamma50FSU * np.log(dose / D50FSU)
    mud = float(np.sum(v * np.array([probit(a) for a in args])))
    mud = min(max(mud, 1e-15), 1 - 1e-15)
    return probit((-math.log(-math.log(mud)) + math.log(-math.log(mu_cr))) / sigma)


def cumulative_to_differential(dose_bins, cum_volume):
    """TsVOutcomeModel::CumulativeToDifferentialDVH — untuk memakai DVH TOPAS."""
    d = list(dose_bins)
    v = list(cum_volume)
    for i in range(1, len(d)):
        v[i - 1] -= v[i]
    return np.array(d[:-1]), np.array(v[:-1])


def tcp_lq_poisson_fractionated(D, alpha, beta, n0, n_fx):
    """Model mekanistik LQ-Poisson dengan fraksinasi d = D/n_fx."""
    d = D / n_fx
    return np.exp(-n0 * np.exp(-n_fx * (alpha * d + beta * d ** 2)))


# =============================================================================
# 3. KURVA DOSIS HOMOGEN
# =============================================================================

def homogeneous(model, D, ab, **p):
    """Hitung outcome untuk dosis homogen D (satu bin DVH, volume = 1)."""
    out = []
    for dose in D:
        dcorr = eqd_ref([dose], ab, N_FRACTIONS, D_REF)
        out.append(model(dcorr, [1.0], **p))
    return np.array(out)


def first_dose_reaching(D, y, level):
    idx = np.argmax(y >= level)
    return float(D[idx]) if y[idx] >= level else float("nan")


def therapeutic_window(D, tcp, ntcp):
    ok = (tcp >= TCP_MIN) & (ntcp <= NTCP_MAX)
    if not ok.any():
        return float("nan"), float("nan")
    return float(D[ok].min()), float(D[ok].max())


def load_topas_csvs(folder):
    curves = []
    for path in sorted(glob.glob(os.path.join(folder, "*.csv"))):
        with open(path, newline="") as fh:
            rows = list(csv.reader(fh))
        header, data = rows[0], np.array(rows[1:], dtype=float)
        for j in range(1, data.shape[1]):
            y = data[:, j]
            if y.max() > 1.0:          # TOPAS melaporkan dalam %
                y = y / 100.0
            label = f"TOPAS {os.path.splitext(os.path.basename(path))[0]}: {header[j]}"
            curves.append((data[:, 0], y, label))
    return curves


def main():
    D = np.linspace(DOSE_MIN, DOSE_MAX, N_POINTS)
    D_oar = D * OAR_DOSE_RATIO

    tcp = homogeneous(tcp_poisson, D, TARGET_AB, TCD50=TCD50, gamma50=GAMMA50_TCP)
    ntcp = {
        "OAR1 Ginjal - LKB": homogeneous(ntcp_lkb, D_oar, OAR1_AB, **LKB_KIDNEY),
        "OAR1 Ginjal - Critical Element": homogeneous(ntcp_critical_element, D_oar, OAR1_AB, **CE_KIDNEY),
        "OAR1 Ginjal - Critical Volume": homogeneous(ntcp_critical_volume, D_oar, OAR1_AB, **CV_KIDNEY),
        "OAR2 Kolon - LKB": homogeneous(ntcp_lkb, D_oar, OAR2_AB, **LKB_COLON),
    }

    # ---- Tabel ringkasan kuantitatif ----
    print(f"Fraksinasi: {N_FRACTIONS} fraksi, d_ref = {D_REF} Gy, "
          f"rasio dosis OAR/target = {OAR_DOSE_RATIO}")
    print(f"Kriteria window: TCP >= {TCP_MIN:.0%} dan NTCP <= {NTCP_MAX:.0%}\n")
    print(f"TCP Poisson target : TCP50 = {first_dose_reaching(D, tcp, 0.5):6.1f} Gy, "
          f"TCP90 = {first_dose_reaching(D, tcp, 0.9):6.1f} Gy")
    if SHOW_LQ_REFERENCE:
        tcp_lq = tcp_lq_poisson_fractionated(D, LQ_ALPHA, LQ_BETA, LQ_N0, N_FRACTIONS)
        print(f"TCP LQ-Poisson ref : TCP50 = {first_dose_reaching(D, tcp_lq, 0.5):6.1f} Gy, "
              f"TCP90 = {first_dose_reaching(D, tcp_lq, 0.9):6.1f} Gy")
    print()
    print(f"{'Model NTCP':34s} {'NTCP5':>7s} {'NTCP50':>7s}   Therapeutic window")
    rows_summary = []
    for name, y in ntcp.items():
        n5, n50 = first_dose_reaching(D, y, 0.05), first_dose_reaching(D, y, 0.5)
        w0, w1 = therapeutic_window(D, tcp, y)
        wtxt = "tidak ada" if math.isnan(w0) else f"{w0:.1f} - {w1:.1f} Gy (lebar {w1 - w0:.1f} Gy)"
        print(f"{name:34s} {n5:7.1f} {n50:7.1f}   {wtxt}")
        rows_summary.append((name, n5, n50, w0, w1))

    # ---- Simpan kurva ke CSV ----
    with open(OUTPUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Dose_Gy", "TCP_Poisson"] + [f"NTCP_{k}" for k in ntcp])
        for i in range(len(D)):
            w.writerow([f"{D[i]:.2f}", f"{tcp[i]:.6f}"] + [f"{ntcp[k][i]:.6f}" for k in ntcp])

    # ---- Plot ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    topas = load_topas_csvs(TOPAS_CSV_DIR) if os.path.isdir(TOPAS_CSV_DIR) else []
    for ax, (name, y), (_, _, _, w0, w1) in zip(axes.flat, ntcp.items(), rows_summary):
        ax.plot(D, tcp, color="tab:blue", lw=2, label="TCP Poisson (analitik homogen)")
        if SHOW_LQ_REFERENCE:
            ax.plot(D, tcp_lq, color="tab:blue", ls=":", lw=1.2, label="TCP LQ-Poisson (ref., difraksinasi)")
        ax.plot(D, y, color="tab:red", lw=2, label=f"NTCP {name.split(' - ')[-1]} (analitik homogen)")
        if not math.isnan(w0):
            ax.axvspan(w0, w1, color="gray", alpha=0.25,
                       label=f"Window {w0:.0f}-{w1:.0f} Gy (dihitung)")
        for x, yy, lab in topas:
            ax.plot(x, yy, "o--", ms=3, lw=1, label=lab)
        ax.set_title(name)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc="center right")
    for ax in axes[1]:
        ax.set_xlabel(f"Dosis total (Gy), {N_FRACTIONS} fraksi")
    for ax in axes[:, 0]:
        ax.set_ylabel("Probabilitas")
    fig.suptitle("TCP-NTCP analitik (dosis homogen) dengan parameter setara TOPAS")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300)
    print(f"\nGrafik disimpan: {OUTPUT_PNG}\nData kurva disimpan: {OUTPUT_CSV}")
    plt.show()


if __name__ == "__main__":
    main()
