import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.optimize

#Lines and regions
lines = ['2-1','3-2']
regions= ['1.1', '2.1', '2.2', '2.3', '4.1', '5.1', '5.2' ]
pixel_region_size_21 = [2.685e3, 4.835e3, 4.292e3, 1.5e2, 2.1353e4, 2.795e3, 1.017e4] #The region size in pixels
pixel_region_size_32 = [3.397e3, 6.088e3, 5.415e3, 1.91e2, 2.1556e4, 3.527e3, 1.249e4] #Important to do them seperately, as some regions get cut off
CO21_pixel_scale= (0.1 * 0.1) #arcsec^2 per pixel
CO32_pixel_scale= (0.089 * 0.089) #arcsec^2 per pixel



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

    if line =='2-1':
        A= 10**(-6.1605) #In units of s^-1 
    if line =='3-2':
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


def gaussian(x, a, b, c):
    """
    Gaussian for the fit
    """
    return a * np.exp(-0.5 * ((x - b) / c) ** 2)



for i, region in enumerate(regions): 
    print('Region:', region)
    output_column_density = [] #To store the column densities
    for line in lines:
        #Import data
        data=np.genfromtxt('/net/vdesk/data2/bach1/ballieux/master_project_2/data/high_mass_data/G338.93/integrated_spectra/'+ line+'_region'+region+'.txt')

        # Here the right columns are selected
        x_data = np.flip(data[:, 3])
        y_data = np.flip(data[:, 4])

        plt.step(x_data, y_data)
        plt.xlabel('velocity [km/s]')
        plt.ylabel('Flux density [Jy]')
        plt.title(line+' '+region)
        plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_2/data/high_mass_data/G338.93/integrated_spectra/spectrum'+ line+'_' + region +'.pdf', bbox_inches='tight')
        plt.close()

        #use trapz to integrate
        integral = np.trapz(y_data, x_data) 
        e_integral = 0.1 * integral #TODO: uncertainties?? Did not deal with any of it now. 

        # print(f"The integrated line intensity is {integral} [Jy km / s] with an uncertainty of {e_integral} for {line}") 

        if line == '2-1':
            area_region_arcsec = pixel_region_size_21[i] * CO21_pixel_scale #in square arcseconds

        if line == '3-2':
            area_region_arcsec = pixel_region_size_32[i] * CO32_pixel_scale #in square arcseconds

        coldens = column_density(integral, area_region_arcsec)

        #Now divide by the degeneracies
        if line == '2-1':
            degeneracy = 5
        elif line == '3-2': 
            degeneracy = 7


        normalized_column_density= (coldens/degeneracy)
        output_column_density.append(normalized_column_density)

    #The energies as found from the database
    E_1 = 11.5350 * 1.438 #K
    E_2 = 23.0695 * 1.438 #K
    energies = np.array([E_1, E_2])
    # print(output_column_density)

    ln_E_g= np.log(output_column_density)
    #Plotting it, uncertainties follow from prop of errors
    plt.errorbar(energies, ln_E_g, yerr = 0.1, fmt='o', color='rebeccapurple')


    popt, pcov = scipy.optimize.curve_fit(f, energies, ln_E_g, sigma = 0.1, absolute_sigma=True )
    # print('AGN',popt, np.sqrt(np.diag(pcov)))
    plt.plot(energies, f(energies, popt[0], popt[1]), linestyle='--', color='black')


    plt.ylabel(r'$\ln(\frac{N_u}{g_u})$')
    plt.xlabel(r'$E_u$ [K] ')
    # plt.legend()



    T = (-1/popt[0])
    e_T = (T ** 2 * np.sqrt(np.diag(pcov))[0])
    print('T of region in Kelvin:', (-1/popt[0]), 'pm', e_T )
    plt.text(17, ln_E_g[0]+0.05, 'T = '+str(np.round(T,3))+' K')

    #Calculate the partition function
    T_range = [2.725, 5.000, 9.375, 18.75, 37.5, 75, 150, 225, 300, 500, 1000, 2000]
    Q_range = [1.4053, 2.1824, 3.7435, 7.1223, 13.8965, 27.4545, 54.5814, 81.7184, 108.8651, 181.3025, 362.6910, 726.7430]

    #Calculate the column density
    Q = np.interp(T, T_range, Q_range)
    N = np.exp(popt[1]) * Q # I have no idea what units this is in. 
    e_N = N * np.sqrt(np.diag(pcov))[1]
    print('AGN column density total', N,'pm', e_N)
    print("")
    plt.text(17, ln_E_g[1]+0.05, 'N = '+str(np.round(N,1))+ ' cm^-2')

    plt.savefig('/net/vdesk/data2/bach1/ballieux/master_project_2/data/high_mass_data/G338.93/Rotation_diagrams/rotation_diagram'+region+'.pdf', bbox_inches='tight')
    plt.close()