import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math


# Some useful functions:
def witherr(x, err, dig):
    """
    Computes floats with significant digits only

    :param x: experimental value
    :param err: uncertainty of experimental value
    :param dig: number of significant digits that should be displayes
    :return: Returns value, uncertainty, and order of magnitude
    """

    xpower = math.floor(math.log10(abs(x)))
    errpower = math.floor(math.log10(err))

    sign = xpower - errpower + dig - 1

    x = round(x * 10 ** (-xpower), sign)
    err = round(err * 10 ** (-xpower), sign)

    return x, err, xpower


def R2(x, y, fit, param, yerr=np.zeros(1)):
    """
    Computes the Coefficient of determinatin R^2
    The values are weighted by their uncertainty

    :param x: list of independent x values
    :param y: list of experimental y values
    :param yerr: list of uncertainties of experimental y values
    :param fit: fit function
    :param param: fit function parameters
    """

    yfit = fit(x, *param)
    if np.isin(0, yerr):
        ss_res = np.sum((y - yfit)**2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)

    else:
        ymean = np.sum(y / yerr**2) / np.sum(yerr**(-2))
        ss_res = np.sum(((y - yfit) / yerr)**2)
        ss_tot = np.sum(((y - ymean) / yerr)**2)

    return float(1 - (ss_res / ss_tot))


# Some constants:
F = 96485.3329  # C/mol      Faraday constant
R = 8.3144598   # J/(mol K)   gas constant
T = 296.15  # K              estimated room temperature
Terr = 2    # K              estimated uncertainty of temperature
Verr = 0.002    # V         uncertainty of reference electroded


# My fitting function:
def lin(x, m, b):
    return 10 ** (m * x + b)


# Importing sample data:
table = pd.read_excel('samples.xlsx')

# Now I loop over all samples:
for row in range(len(table)):
    if table['done'][row] != 1:
        print('\nsample', table['sample'][row], '...')


# Importing data and setting working directory:
        sample = table['sample'][row]
        pH = table['ph'][row]
        Ru = table['Ru'][row]
        comment = table['comments'][row]
        if str(comment) in ['nan', 'NAN', 'Nan']:
            comment = ''

        folder = 'nofolder'
        for file in os.listdir('expdata/' + sample):
            if 'Tafel' in file:
                folder = 'expdata/' + sample + '/' + file + '/'
            elif 'tafel' in file:
                folder = 'expdata/' + sample + '/' + file + '/'

        if folder == 'nofolder':
            print('no tafel folder found.')
            break

# Organizing data into dataframes and lists:
        print('reading data')

        CV = pd.read_csv(folder + 'CV.DTA', sep='	', skiprows=53, header=None)
        CV = CV.drop([0, 1, 5, 7, 8, 9], axis=1)
        CV.columns = ['t', 'Vf', 'Im', 'Sig']
        CV['Vf'] = pd.to_numeric(CV['Vf'], errors='coerce')
        CV = CV.dropna().astype('float').reset_index()

        CV_final = pd.read_csv(folder + 'CV_final.DTA', sep='	', skiprows=53, header=None)
        CV_final = CV_final.drop([0, 1, 5, 7, 8, 9], axis=1)
        CV_final.columns = ['t', 'Vf', 'Im', 'Sig']
        CV_final['Vf'] = pd.to_numeric(CV_final['Vf'], errors='coerce')
        CV_final = CV_final.dropna().astype('float').reset_index()

        chronop_1 = pd.read_csv(folder + 'CHRONOP_1.DTA', sep='	', skiprows=51, header=None)
        chronop_1 = chronop_1.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_1.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_2 = pd.read_csv(folder + 'CHRONOP_2.DTA', sep='	', skiprows=51, header=None)
        chronop_2 = chronop_2.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_2.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_3 = pd.read_csv(folder + 'CHRONOP_3.DTA', sep='	', skiprows=51, header=None)
        chronop_3 = chronop_3.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_3.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_4 = pd.read_csv(folder + 'CHRONOP_4.DTA', sep='	', skiprows=51, header=None)
        chronop_4 = chronop_4.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_4.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_5 = pd.read_csv(folder + 'CHRONOP_5.DTA', sep='	', skiprows=51, header=None)
        chronop_5 = chronop_5.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_5.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_6 = pd.read_csv(folder + 'CHRONOP_6.DTA', sep='	', skiprows=51, header=None)
        chronop_6 = chronop_6.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_6.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_7 = pd.read_csv(folder + 'CHRONOP_7.DTA', sep='	', skiprows=51, header=None)
        chronop_7 = chronop_7.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_7.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_8 = pd.read_csv(folder + 'CHRONOP_8.DTA', sep='	', skiprows=51, header=None)
        chronop_8 = chronop_8.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_8.columns = ['t', 'Vf', 'Im', 'Sig']

        chronop_9 = pd.read_csv(folder + 'CHRONOP_9.DTA', sep='	', skiprows=51, header=None)
        chronop_9 = chronop_9.drop([0, 1, 5, 7, 8, 9], axis=1)
        chronop_9.columns = ['t', 'Vf', 'Im', 'Sig']

        steps = [chronop_1, chronop_2, chronop_3, chronop_4, chronop_5, chronop_6, chronop_7, chronop_8, chronop_9]


# A comparison of CV curves before and after electrolysis:
        print("comparing CV's")

        # slicing the CV's into up and down parts:
        intervals_before = [0, ]  # start points for 'up', 'down', 'up', 'down' etc.
        for i in range(2, len(CV)):
            if (CV['Vf'][i] - CV['Vf'][i - 1]) * (CV['Vf'][i - 1] - CV['Vf'][i - 2]) < 0:
                intervals_before.append(i - 1)

        intervals_after = [0, ]
        for i in range(2, len(CV_final)):
            if (CV_final['Vf'][i] - CV_final['Vf'][i - 1]) * (CV_final['Vf'][i - 1] - CV_final['Vf'][i - 2]) < 0:
                intervals_after.append(i - 1)


# Potential correction (vs NHE and open circuit potential):
        CV['Vf'] = CV['Vf'] + 0.2 - CV['Im'] * Ru
        CV_final['Vf'] = CV_final['Vf'] + 0.2 - CV_final['Im'] * Ru


# conversion A -> mA
        CV['Im'] *= 1000
        CV_final['Im'] *= 1000

        CV_before = CV[intervals_before[2]:intervals_before[4]]
        CV_after = CV_final[intervals_after[2]:]


# Plots for a comparison of CV before and after
        E0 = 1.229 - 0.059 * pH

        fig, ax = plt.subplots()

        ax.plot(CV_before['Vf'], CV_before['Im'], label='before')
        ax.plot(CV_after['Vf'], CV_after['Im'], label='after')
        ax.set_xlabel('E  / V vs NHE')
        ax.set_ylabel(r'j  / A cm$^{{-2}}$')
        #if comment == '':
        #    ax.set_title('CV (pH {})'.format(pH))
        #else:
        #    ax.set_title('CV (pH {} - {})'.format(pH, comment))
        ax.legend()

        plt.axhline(0, color='grey', lw=1)
        ax.axvline(E0, color='grey', ls='--', lw=1)
        #plt.grid()
        plt.savefig('plots/CV_{}.png'.format(sample), dpi=300)
        # plt.show(fig)
        plt.close(fig)


# Plots for the different steps
        print('plotting E/t curves for current steps')

        fig, ax = plt.subplots(3, 3, figsize=(10, 8))

        counter = 0
        for i in range(3):
            for j in range(3):
                counter += 1
                counts = [counter * 2 - 1, counter * 2]
                ax[i, j].plot(steps[i * 3 + j]['t'], steps[i * 3 + j]['Vf'])
                if counter < 9:
                    ax[i, j].set_title('CP Steps {} and {}'.format(counts[0], counts[1]))
                elif counter == 9:
                    ax[i, j].set_title('CP Step {}'.format(counts[0]))
                ax[i, j].set_xlabel('t   / s')
                ax[i, j].set_ylabel('E   / V vs Ag/AgCl')

        plt.tight_layout()
        plt.savefig('plots/steps_{}.png'.format(sample), dpi=300)
        # plt.show()
        plt.close(fig)


# For data analysis, only the last 1000 values of each of the steps are used:
        print('generating tafel plot')

        lastvalues = []

        for step in steps:
            for i in range(5, len(step) - 1):
                if abs(step['Sig'][i] - step['Sig'][i + 1]) > 0.001:
                    last = i
                    lastvalues.append([step['Vf'][last - 1000: last], step['Im'][last - 1000: last]])
                    break

            lastvalues.append([step['Vf'][-1000:], step['Im'][-1000:]])

        lastvalues = np.array(lastvalues)


# The points in the tafel plot are calculated as a mean of these last values:
        tafel = []

        for i in range(len(lastvalues)):
            tafel.append([[lastvalues[i, 0].mean(), lastvalues[i, 0].std()],
                          [lastvalues[i, 1].mean(), lastvalues[i, 1].std()]])

        tafel = np.array(tafel)  # [[V, deltaV], [I, deltaI]]


# Potential correction
        # eta = E(AgCl) + 0.2 V - Ru * I - (1.229 V - 0.059 * pH)
        for i in range(len(tafel)):
            tafel[i, 0, 0] = tafel[i, 0, 0] + 0.2 - tafel[i, 1, 0] * Ru - 1.229 + 0.059 * pH

            # errors:
            tafel[i, 0, 1] = np.sqrt(tafel[i, 0, 1] ** 2 + (tafel[i, 1, 1] * Ru) ** 2)


# Reading fit points
        fitpoints = pd.read_csv('scripts/fitpoints.csv', header=0, index_col=0)

        if sample in fitpoints.columns:
            startup, stopup, startdown, stopdown = fitpoints[sample]

        else:
            # Simple plot of data points on tafel plot to select points for fit:
            plt.errorbar(tafel[:, 0, 0], tafel[:, 1, 0], xerr=tafel[:, 0, 1], yerr=tafel[:, 1, 1],
                         ls='', marker='.')
            plt.yscale('log')
            plt.title('Select points for fit:'.format(sample))
            plt.xlabel(r'$\eta$   / V')
            plt.ylabel(r'j   / A cm$^2$)')
            plt.show()

            # Selecting data points for fit:
            startup = int(input('\nfit upwards from point: ')) - 1
            stopup = int(input("until point: "))
            startdown = int(input('fit downwards from point: ')) - 1
            stopdown = int(input("until point: "))

            fitpoints = fitpoints.join(pd.DataFrame([startup, stopup, startdown, stopdown], columns=[sample]))
            fitpoints.to_csv('scripts/fitpoints.csv')

            plt.close(fig)


# Fitting:
        print('fitting, processing and saving data')

        pointsup = tafel[startup:stopup]
        pointsdown = tafel[startdown:stopdown]
        pointslist = [pointsdown, pointsup]

        slopes = []
        j0s = []

        for i in range(2):
            param, covar = curve_fit(lin, pointslist[i][:, 0, 0], pointslist[i][:, 1, 0], sigma=pointslist[i][:, 1, 1])

            m, b = param
            merr, berr = np.sqrt(np.diag(covar))


# Calculating j0, alpha, R^2, Chi^2, tafel slope in mV/dec:
            # Errors are calculated according to Gauss.
            determination = R2(pointslist[i][:, 0, 0], pointslist[i][:, 1, 0], lin, param, yerr=pointslist[i][:, 1, 1])

            merr = merr / m  # to relative uncertainty
            m = 1000 / m
            merr = m * merr  # back to absolute uncertainty

            j0 = 10 ** b
            j0err = j0 * np.log(10) * berr

# Saving slopes, errors and R^2:
            slopes.append([m, merr, determination])
            j0s.append([j0, j0err])


# Computing average slope:
        if stopup - startup < 3 or stopdown - startdown < 3: # only if there are too few points to calculate std
            average = (slopes[0][0] + slopes[1][0]) / 2
            averageerr = float('nan')
        else:
            average = (slopes[0][0] / slopes[0][1]**2 + slopes[1][0] / slopes[1][1]**2) /\
                      (1 / slopes[0][1]**2 + 1 / slopes[1][1]**2)
            averageerr = 1 / math.sqrt(1 / slopes[0][1]**2 + 1 / slopes[1][1]**2)

        j0 = (j0s[0][0] / j0s[0][1]**2 + j0s[1][0] / j0s[1][1]**2) / (1 / j0s[0][1]**2 + 1 / j0s[1][1]**2)
        j0err = 1 / math.sqrt(1 / j0s[0][1]**2 + 1 / j0s[1][1]**2)


# Plotting Tafel plot:
        x = np.linspace(pointslist[i][:, 0, 0].min() - 0.03, pointslist[i][:, 0, 0].max() + 0.02, 1000)
        fit = lin(x, *param)

        fig, ax = plt.subplots()

        ax.errorbar(tafel[:, 0, 0], tafel[:, 1, 0], xerr=tafel[:, 0, 1], yerr=tafel[:, 1, 1],
                    ls='', marker='.', label='all points')  # all points
        ax.errorbar(pointslist[1][:, 0, 0], pointslist[1][:, 1, 0], xerr=pointslist[1][:, 0, 1],
                    yerr=pointslist[1][:, 1, 1], ls='', marker='o', label='points used for fit')
        ax.plot(x, fit, color='grey')

        # Plotting values as text:
        if abs(m) < 0.01:
            m, merr, mpower = witherr(m, merr, 1)
            ax.text(0.05, 0.70, r'slope = ({} $\pm$ {}) $\cdot 10^{{{}}}$  mV/dec'.format(m, merr, mpower),
                    transform=ax.transAxes)
        else:
            if merr < 1:
                significant = int(- math.floor(math.log10(abs(merr))))
                m = round(m, significant)
                merr = round(merr, significant)
            else:
                m = round(m)
                merr = round(merr)
            ax.text(0.05, 0.70, r'slope = ({} $\pm$ {})  mV/dec'.format(m, merr), transform=ax.transAxes)

        #determination = round(determination, 3)
        #ax.text(0.05, 0.62, r'R$^2$ = {}'.format(determination), transform=ax.transAxes)

        ax.set_yscale('log')
        #if comment == '':
        #    ax.set_title('Tafel Plot (pH {})'.format(pH))
        #else:
        #    ax.set_title('Tafel Plot (pH {} - {})'.format(pH, comment))
        ax.set_xlabel(r'$\eta$   / V')
        ax.set_ylabel(r'j   / A cm$^{{-2}}$')
        plt.legend()
        plt.savefig('plots/tafelplot_{}.png'.format(sample), dpi=300)
        # plt.show()
        plt.close(fig)


# Saving data to table:
        table.loc[row, 'slope up'] = slopes[1][0]
        table.loc[row, 'slope up err'] = slopes[1][1]
        table.loc[row, 'R^2 up'] = slopes[1][2]
        table.loc[row, 'slope down'] = slopes[0][0]
        table.loc[row, 'slope down err'] = slopes[0][1]
        table.loc[row, 'R^2 down'] = slopes[0][2]
        table.loc[row, 'average slope'] = average
        table.loc[row, 'average err'] = averageerr
        table.loc[row, 'j0'] = j0
        table.loc[row, 'j0 err'] = j0err


        table.loc[row, 'done'] = 1


# Saving up points used for fit from Tafel plot
        tafelpoints = pd.read_csv('scripts/tafelpoints.csv', header=0, index_col=0)
        if sample not in tafelpoints.columns:
            #points = np.transpose(np.array([pointsup[:, 0, 0], pointsup[:, 1, 0]]))
            tafelpoints = tafelpoints.append(pd.DataFrame([pointsup[:, 0, 0], pointsup[:, 1, 0]],
                                                          index=[sample + 'x', sample + 'y']))
            tafelpoints.to_csv('scripts/tafelpoints.csv')


# saving table
        table.to_excel('samples.xlsx', engine='xlsxwriter', index=False)
