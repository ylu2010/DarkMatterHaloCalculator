import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, sys, glob, time, cStringIO
import parameter as par
import cosmology_lib as Cosmo
reload(Cosmo)

def compute(m, z):
    rvir = rvir_mvir(m, z=z)
    vvir = vvir_mvir(m, z=z)
    tvir = tvir_vvir(vvir)
    plot = plot_halo_histories(z, m, 0.3)
    #plot = plot_test(z, m, 0.3)
    return rvir, vvir, tvir, plot

def rvir_mvir(mvir, z=0.0):
    h = par.Hubble/100.
    rhovir_z = Cosmo.delta_vir(z) * Cosmo.rho_crit(z)*h*h
    rvir = (3.*mvir / 4 /np.pi/rhovir_z)**(1./3) * 1e3 # in unit of kpc
    return rvir

def vvir_mvir(mvir, z=0.0):
    h = par.Hubble/100.
    rhovir_z = Cosmo.delta_vir(z) * Cosmo.rho_crit(z)*h*h
    vvir = np.sqrt(Cosmo.Grav) * (4./3*np.pi*rhovir_z)**(1./6) * mvir**(1./3);
    return vvir

def tvir_vvir(vvir):
    tvir = 35.9 * vvir * vvir
    return tvir

def mah_wechsler(z, g=0.3):
    return np.exp(-g*z)

def plot_halo_histories(z0, m0, g, resolution=100):
    """Return filename of plot of the MAH."""
    z = np.linspace(0, 20, resolution+1)
    mnorm = m0 / mah_wechsler(z0, g)
    mah = mah_wechsler(z, g) * mnorm
    mah0 = mah[0]

    rvirh = rvir_mvir(mah, z)
    rvirh0 = rvirh[0]
    rvirhmin = np.min(rvirh)

    vvirh = vvir_mvir(mah, z)
    vvirh0 = vvirh[0]
    vvirhmin = np.min(vvirh)

    tvirh = tvir_vvir(vvirh)
    tvirh0 = tvirh[0]
    tvirhmin = np.min(tvirh)

    # start plotting
    try:
        # Store App Engine's modified stdout so we can restore it later
        gae_stdout = sys.stdout

        plt.clf()
        plt.figure()  # needed to avoid adding curves in plot
        font = {'family' : 'serif',
            'weight' : 'normal',
            'size'   : 16}
        matplotlib.rc('font', **font)
        nrow = 2
        ncol = 2
        fig, ax = plt.subplots(nrow, ncol, sharex=False, sharey=False)
        fig.set_size_inches(7.,6)

        # Panel 00
        ir = 0
        ic = 0
        xran = [1, 20]
        yran = [mah0*0.002, mah0*2]
        ax[ir,ic].plot(z+1, mah, '-', color='blue')
        ax[ir,ic].plot([z0+1,z0+1], yran, ':', color = 'red')
        #ax.plot(nstu_vec, d_alpha_nattr10, 'o-', color='blue')
        #plt.plot(z, u)
        ax[ir,ic].set_xscale('log')
        ax[ir,ic].set_yscale('log')
        #ax[ic].set_xlim([1, 600])
        #ax[ic].set_ylim(0,0.1)
        ax[ir,ic].set_xlim(xran)
        ax[ir,ic].set_ylim(yran)
        #ax.set_xlabel(u'$N_{stu}$')
        ax[ir,ic].set_ylabel(r'$M_{\rm vir}$'+' [$M_{\odot}$]')
        ax[ir,ic].set_xlabel(u'$1+z$')
        #ax[ic].legend(loc=0, frameon=False, fontsize=16)

        # Panel 01
        ir = 0
        ic = 1
        xran = [1, 20]
        yran = [rvirhmin/3, rvirh0*3]
        ax[ir,ic].plot(z+1, rvirh, '-', color='blue')
        ax[ir,ic].plot([z0+1,z0+1], yran, ':', color = 'red')
        #ax.plot(nstu_vec, d_alpha_nattr10, 'o-', color='blue')
        #plt.plot(z, u)
        ax[ir,ic].set_xscale('log')
        ax[ir,ic].set_yscale('log')
        #ax[ic].set_xlim([1, 600])
        #ax[ic].set_ylim(0,0.1)
        ax[ir,ic].set_xlim(xran)
        ax[ir,ic].set_ylim(yran)
        #ax.set_xlabel(u'$N_{stu}$')
        ax[ir,ic].set_ylabel(r'$R_{\rm vir}$'+' [$kpc$]')
        ax[ir,ic].set_xlabel(u'$1+z$')
        #ax[ic].legend(loc=0, frameon=False, fontsize=16)

        # Panel 10
        ir = 1
        ic = 0
        xran = [1, 20]
        yran = [vvirhmin/3, vvirh0*3]
        ax[ir,ic].plot(z+1, vvirh, '-', color='blue')
        ax[ir,ic].plot([z0+1,z0+1], yran, ':', color = 'red')
        #ax.plot(nstu_vec, d_alpha_nattr10, 'o-', color='blue')
        #plt.plot(z, u)
        ax[ir,ic].set_xscale('log')
        ax[ir,ic].set_yscale('log')
        #ax[ic].set_xlim([1, 600])
        #ax[ic].set_ylim(0,0.1)
        ax[ir,ic].set_xlim(xran)
        ax[ir,ic].set_ylim(yran)
        #ax.set_xlabel(u'$N_{stu}$')
        ax[ir,ic].set_ylabel(r'$V_{\rm c}$'+' [$km\,s^{-1}$]')
        ax[ir,ic].set_xlabel(u'$1+z$')
        #ax[ic].legend(loc=0, frameon=False, fontsize=16)

        # Panel 11
        ir = 1
        ic = 1
        xran = [1, 20]
        yran = [tvirhmin/3, tvirh0*3]
        ax[ir,ic].plot(z+1, tvirh, '-', color='blue')
        ax[ir,ic].plot([z0+1,z0+1], yran, ':', color = 'red')
        #ax.plot(nstu_vec, d_alpha_nattr10, 'o-', color='blue')
        #plt.plot(z, u)
        ax[ir,ic].set_xscale('log')
        ax[ir,ic].set_yscale('log')
        #ax[ic].set_xlim([1, 600])
        #ax[ic].set_ylim(0,0.1)
        ax[ir,ic].set_xlim(xran)
        ax[ir,ic].set_ylim(yran)
        #ax.set_xlabel(u'$N_{stu}$')
        ax[ir,ic].set_ylabel(r'$T_{\rm vir}$'+' [$K$]')
        ax[ir,ic].set_xlabel(u'$1+z$')
        #ax[ic].legend(loc=0, frameon=False, fontsize=16)

        #plt.title('MAH')
        fig.tight_layout()
        # Redirect stdout to a StringIO object
        new_stdout = cStringIO.StringIO()
        sys.stdout = new_stdout

        #plt.savefig(new_stdout, format="png")
        plt.savefig(new_stdout)
        #plt.clf()
        plotfile = """data:image/png;base64,%s""" % new_stdout.getvalue().encode("base64").strip()
        sys.stdout = gae_stdout
        return plotfile
    finally:
        plt.clf()

def plot_test(z, m, g):
    try:
        # Store App Engine's modified stdout so we can restore it later
        gae_stdout = sys.stdout

        plt.title("Dynamic PNG")
        plt.plot([0,1,2,3,4], [2,3,4,3,2], 'o-')

        # Redirect stdout to a StringIO object
        new_stdout = cStringIO.StringIO()
        sys.stdout = new_stdout

        #plt.savefig(new_stdout, format="png")
        plt.savefig(new_stdout)
        #plt.clf()
        plotfile = """data:image/png;base64,%s""" % new_stdout.getvalue().encode("base64").strip()
        sys.stdout = gae_stdout
        return plotfile
    finally:
        plt.clf()
