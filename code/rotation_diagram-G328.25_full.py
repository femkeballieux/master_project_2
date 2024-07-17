import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.optimize

#Lines and regions
path = '/net/vdesk/data2/bach1/ballieux/master_project_2/data/high_mass_data/G328.25/'
lines = ['CO21','CO32']
regions= ['redshifted_outflow','blueshifted_outflow']
pixel_region_size_21 = [2.787e3, 5.938e3] #The region size in pixels
pixel_region_size_32 = [1.0759e4, 1.5812e4] #Important to do them seperately, as some regions get cut off
CO21_pixel_scale= (0.11 * 0.11) #arcsec^2 per pixel
CO32_pixel_scale= (0.056 * 0.056) #arcsec^2 per pixel


def column_density(integrated_flux_jansky, area_region_arcsec):
    """
    Enter the integrated flux in Jy km / s
    The radius in arcsecond of the emitted region 
    Adding the target is only for printing statements
    Adding the line decides the Einstein coefficient

    Returns the column density in [cm^-2]
    """

    integrated_flux_erg = integrated_flux_jansky * 10 ** (-18) # [erg cm^-1 s^-1]

    # beam_arcsec = np.pi * radius_arcsec ** 2 / (4 * np.log(2)) # in arcseconds^2
    area_steradians = area_region_arcsec * (1/3600**2) * (2 * np.pi / 360)**2 #In steradians

    if line =='CO21':
        A= 10**(-6.1605) #In units of s^-1 
    if line =='CO32':
        A= 10**(-5.6026) #In units of s^-1

    h = 6.625*10**(- 27) # [erg s]
    c = 3 * 10**10       #[cm/s]

    column_density = 4 * np.pi * integrated_flux_erg / (A * area_steradians * h * c) #[cm^-2]
    return column_density


def f(x,a,b):
    """
    straight line function
    """
    return a * x + b


# def gaussian(x, a, b, c):
#     """
#     Gaussian for the fit
#     """
#     return a * np.exp(-0.5 * ((x - b) / c) ** 2)



for i, region in enumerate(regions): 
    print('Region:', region)
    output_column_density = [] #To store the column densities
    for line in lines:
        #Import data
        data=np.genfromtxt(path+'integrated_spectra/'+region+'_'+line+'.txt')

        # Here the right columns are selected
        x_data = np.flip(data[:, 3])
        y_data = np.flip(data[:, 4])

        plt.step(x_data, y_data)
        plt.xlabel('velocity [km/s]')
        plt.ylabel('Flux density [Jy]')
        plt.title(line+' '+region)
        plt.savefig(path+'integrated_spectra/spectrum'+ line+'_' + region +'.pdf', bbox_inches='tight')
        plt.close()

        #use trapz to integrate
        integral = np.trapz(y_data, x_data) 
        e_integral = 0.1 * integral #TODO: uncertainties?? Did not deal with any of it now. 

        # print(f"The integrated line intensity is {integral} [Jy km / s] with an uncertainty of {e_integral} for {line}") 

        if line == 'CO21':
            area_region_arcsec = pixel_region_size_21[i] * CO21_pixel_scale #in square arcseconds

        if line == 'CO32':
            area_region_arcsec = pixel_region_size_32[i] * CO32_pixel_scale #in square arcseconds
        print(area_region_arcsec)
        coldens = column_density(integral, area_region_arcsec)

        #Now divide by the degeneracies
        if line == 'CO21':
            degeneracy = 5
        elif line == 'CO32': 
            degeneracy = 7


        normalized_column_density= (coldens/degeneracy)
        output_column_density.append(normalized_column_density)

    #The energies as found from the database, upper energies
    E_1 = 11.5350 * 1.438 #K
    E_2 = 23.0695 * 1.438 #K
    
    energies = np.array([E_1, E_2])

    #The column densities
    ln_E_g = np.log(output_column_density)
    E_ln_E_g = 0.1 #From prop of errors, this is absolute uncertainty

    #Plotting it
    plt.errorbar(energies, ln_E_g, yerr = E_ln_E_g, fmt='o', color='rebeccapurple')

    a = (ln_E_g[1] - ln_E_g[0]) / (energies[1] - energies[0])
    b = ln_E_g[0] - a * energies[0]
    E_a = np.abs(np.sqrt(2) * E_ln_E_g  / (energies[1] - energies[0]) )
    E_b = np.sqrt( E_ln_E_g ** 2 + (E_a * energies[0]) ** 2)
    plt.plot(energies, f(energies, a, b), linestyle='--', color='black')
    plt.ylabel(r'$\ln(\frac{N_u}{g_u})$')
    plt.xlabel(r'$E_u$ [K] ')
    # plt.legend()

    T = (-1/a)
    e_T = T ** 2 * np.abs(E_a )
    print('T of region in Kelvin:', T, 'pm', e_T )
    plt.text(17, ln_E_g[0] + 0.05, 'T = '+str(np.round(T,0))+ ' $\pm$ ' + str(np.round(e_T, 0))+' K') #position is first 2 nummbers

    #Calculate the partition function
    T_range = [2.725, 5.000, 9.375, 18.75, 37.5, 75, 150, 225, 300, 500, 1000, 2000]
    Q_range = [1.4053, 2.1824, 3.7435, 7.1223, 13.8965, 27.4545, 54.5814, 81.7184, 108.8651, 181.3025, 362.6910, 726.7430]

    #Calculate the column density
    Q = np.interp(T, T_range, Q_range)
    N = np.exp(b) * Q #cm^-2
    e_N = np.abs(N * E_b) #follows from prop of errors

    print('AGN column density total', N,'pm', e_N)
    print("")

    plt.text(17, ln_E_g[1]+0.05, 'N = '+str(np.round(N,-16))+ ' $\pm$ ' + str(np.round(e_N, -16))+' $cm^{-2}$' ) #position is first 2 nummbers
    plt.title(region)
    plt.savefig(path+'rotation_diagrams/rotation_diagram'+region+'.png', bbox_inches='tight')
    plt.close()