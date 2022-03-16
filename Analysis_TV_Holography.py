# -*- coding: utf-8 -*-
"""
Spyder Editor
This is the code for TV_Hologpraphy experiment 

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import easygui
import csv
from operator import itemgetter
import os.path
from os import listdir
from scipy.optimize import curve_fit


# Specific constants for each data set
OFFSET = [500,3000]
HEIGHT = 12.5
DISTANCE = 150
MID_POINT = 1500
FIT_FUNCTION = 1
LASER_WAVELENGTH = 632.8 #nm
COLOUR_FACTOR = 1 # for red 1 ; for green 1.7

def case(voltage):
    
    global DISTANCE
    global MID_POINT
    global COLOUR_FACTOR
    
    MID_POINT = 1100
    
    V1 = int(voltage[0].replace('V',''))
    V2 = int(voltage[1].replace('V',''))
    diff = np.abs(V1 - V2)
    
    if diff < 8:
        DISTANCE = 850/ COLOUR_FACTOR
        
    elif 8 <= diff <= 12:
        DISTANCE = 430/ COLOUR_FACTOR
        
    elif 12 < diff <= 16:
        DISTANCE = 280/ COLOUR_FACTOR
        
    elif 16 < diff < 20:
        DISTANCE = 230/ COLOUR_FACTOR
        
    elif 20 <= diff < 23:
        DISTANCE = 180/ COLOUR_FACTOR
        
    elif 23 <= diff < 28:
        DISTANCE = 150/ COLOUR_FACTOR
        
    elif 28 <= diff < 32:
        DISTANCE = 130/ COLOUR_FACTOR
        
    elif 32 <= diff < 36:
        DISTANCE = 110/ COLOUR_FACTOR
        
    elif 36 <= diff < 40:
        DISTANCE = 100/ COLOUR_FACTOR
    
    elif 40 <= diff < 44:
        DISTANCE = 90/ COLOUR_FACTOR
        
    elif 60 < diff <= 65:
        DISTANCE = 70/ COLOUR_FACTOR
        
    elif 65 < diff < 70:
        DISTANCE = 60/ COLOUR_FACTOR
        
    elif 70 <= diff < 80:
        DISTANCE = 45/ COLOUR_FACTOR
        
    elif 80 <= diff <= 100:
        DISTANCE = 30/ COLOUR_FACTOR
        
    elif 100 < diff <= 150:
        DISTANCE = 25/ COLOUR_FACTOR
        
    else:
        DISTANCE = 80/ COLOUR_FACTOR


def pixel_mm(pix):
    return pix * 0.0043  # mm


def import_csv(path):
    """
    

    Returns:
        x - 1st column
        y - 2nd column
    -------
    None.

    """
    
    x = []
    y = []
    with open(str(path)) as file_name:
        file_read = csv.reader(file_name)
        
        for row in file_read:
            if not row:
                row = [0,0]
            x.append(int(row[0]))
            y.append(float(row[1]))
    x.pop()
    y.pop()
    return np.array(x),np.array(y)
    

def filter_data(data,sep_fac=10):
    
    y_array = data[1]
    x_data = []
    y_data = []
    
    for i in range(int(len(y_array)/sep_fac)-1):
        x_data.append(i*sep_fac + int(sep_fac/2))
        bb = 0
        for j in range(sep_fac):
            bb += y_array[i*sep_fac + int(sep_fac/2)]
            
        y_data.append(bb/sep_fac)
        
    return np.array(x_data), np.array(y_data)

def plot_data(x, y, peaks, voltages = ['give', 'title']):
    
    
    plt.errorbar(x, y,fmt=',')
    
    peak_value = []
    peak_position = []
    for i in peaks:
        peak_value.append(y[i])
        peak_position.append(x[i])
    
    plt.plot(peak_position, peak_value, "x")
    plt.axis('off')
    plt.title(voltages[0] + '_' + voltages[1])
    plt.show()
    
    
def find_maxima(data,sep_fac=10):
    global DISTANCE
    global HEIGHT
    global OFFSET
    
    peaks,_ = find_peaks(data[1],height = HEIGHT,distance = DISTANCE/sep_fac)
    peak = []
    
    for i in peaks:
        if i > OFFSET[0]/sep_fac and i < OFFSET[1]/sep_fac:
            peak.append(i)
        
    return peak

def possibe_minima(peaks):
    
    minima = []
    for i in range(len(peaks)-1):
        minima.append(int((peaks[i]+peaks[i+1])/2))
    
    return minima

def find_spacing(peaks):
    
    spacing = []
    for i in range(len(peaks)-1):
        spacing.append(np.abs(pixel_mm(peaks[i])-pixel_mm(peaks[i+1])))
        
    return spacing

def find_amplitude(maxima,minima,y_data):
    
    amp =[]
    for i in range(len(minima)):
        amp.append(y_data[maxima[i]]-y_data[minima[i]])
        amp.append(y_data[maxima[i+1]]-y_data[minima[i]])
    return amp

def find_points_around(peak, n=1, sepfac=10):
    ### Look for point near to MID_POINT, then look for the nearest point to the near point
    global MID_POINT
    
    peaks = []
    a = [[np.abs(x - MID_POINT/sepfac), x] for x in peak]
    a = sorted(a, key=itemgetter(0))
    new_point = a[0][1]
    del(a)
    
    a = [[np.abs(x - new_point), x] for x in peak]
    a = sorted(a, key=itemgetter(0))
    for i in range(n):
        peaks.append(a[i][1])
    
    peaks = sorted(peaks)
    return peaks

def do_analysis(data,n=3,sepfac=10):
    
    all_peaks = find_maxima(data)
    x_values = data[0]
    peaks = find_points_around(all_peaks,n,sepfac) # n maximas near MID_POINT
    spacing = [element * sepfac for element in find_spacing(peaks)] #distance between maximas
    min_peaks = possibe_minima(peaks) # Minimas 
    minima_peaks = []
    
    for i in range(len(min_peaks)):
        
        a = [[np.abs(x - min_peaks[i]), x] for x in x_values]
        
        a = sorted(a, key=lambda x: x[0])
        minima_peaks.append(int(a[0][1]))
        
    amp = find_amplitude(peaks,minima_peaks,data[1]) # amplitude
    
    return amp, spacing

def analyze(results):
    spacing = []
    amp = []
    voltages = []
    
    mean_amp_list = []
    mean_spacing_list = []
    voltage_diff_list = []
    
    for i in range(len(results)):
        amp.append(results[i][0][0])
        spacing.append(results[i][0][1])
        voltages.append(results[i][1])
    
    for i in range(len(results)):

        #amp1 = max(amp[i])
        #amp[i].remove(amp1)
        #amp2 = max(amp[i])

        #mean_amp = (amp1 + amp2)/2
        try:
            mean_amp = (max(amp[i])+min(amp[i]))/2
        except:
            mean_amp = np.mean(amp[i])
            
        mean_spacing = np.mean(spacing[i])
        
        ### Image - Reference Image
        volt_diff = np.abs(int(voltages[i][1].replace('V','')) - int(voltages[i][0].replace('V','')))
        
        mean_amp_list.append(mean_amp)
        mean_spacing_list.append(mean_spacing)
        voltage_diff_list.append(volt_diff)
        
    return mean_amp_list, mean_spacing_list, voltage_diff_list
        

def angle_of_PZT(spacing):
    global LASER_WAVELENGTH
    wavelen = LASER_WAVELENGTH * 10**(-6) #mm
    angles = []
    
    for sp in spacing:
        angles.append(wavelen / sp)
        
    return angles

def result_plot(result):
    values = result
    
    #spacing against voltage
    plt.plot(values[2],values[1],'o')
    plt.title('spacing')
    plt.show()
    # amplitude against voltage
    plt.plot(values[2],values[0],'o')
    plt.title('amplitude')
    plt.show()
    # angle against voltage
    plt.plot(values[2], angle_of_PZT(values[1]),'o')
    plt.show()
    

def plot_error(x,y,x_error = None, y_error = None, label =''):
    plt.errorbar(x,y,xerr = x_error,yerr=y_error,fmt='.', label=label) #xerr = x_error,yerr=y_error,

def analyse_multiple(PATH,sepfac = 8):
    global COLOUR_FACTOR
    global LASER_WAVELENGTH
    
    #PATH = "D:\\TVHolography\\NewFolder\\photos\\Analyse"
    input_path = PATH
    filenames = listdir(input_path)
    folder_list = [ filename for filename in filenames]
    all_results = []
    for file in folder_list: # what day
        new_path = input_path + "\\" + file
        new_files = [filename for filename in listdir(new_path)]

        for file2 in new_files: # files in that day
            final_path = new_path + '\\' + file2
            colour_laser = file2.split('_')[0]
            
            if colour_laser == 'greenlaser':
                COLOUR_FACTOR = 1.5
                #LASER_WAVELENGTH = 532
            else:
                COLOUR_FACTOR = 1
                #LASER_WAVELENGTH = 632.8
                
            label = file + ' - ' + file2
            res = main_analysis(final_path, sep_fac = sepfac)
            all_results.append([res,label])
            
    return all_results #2D array
    

def main_analysis(PATH,sep_fac = 6):
    
    results = []
    
    filenames = listdir(PATH+'\\results')
    csv_files_list = [ filename for filename in filenames if filename.endswith( ".csv" ) ]
    
    for file in csv_files_list:
        new_path = PATH +'\\results' +'\\'+ file
        
        voltage_string = [file.split('_')[0],file.split('_')[2]] # 0 and 2
        case(voltage_string)
        data = filter_data(import_csv(new_path),sep_fac)


        try:
            all_peaks = find_maxima(data,sep_fac)
            peaks = find_points_around(all_peaks,n=3, sepfac= sep_fac) # n maximas near MID_POINT
            ana = do_analysis(data,3,sep_fac)
            results.append([ana,voltage_string])
            #plot_data(data[0],data[1],peaks,voltage_string)
            del ana
        except:
            try:
                all_peaks = find_maxima(data,sep_fac)
                peaks = find_points_around(all_peaks,n=2, sepfac= sep_fac) # n maximas near MID_POINT
                ana = do_analysis(data,2,sep_fac)
                results.append([ana,voltage_string])
                #plot_data(data[0],data[1],peaks,voltage_string)
                del ana
            except:
                print('exception')
                
    return results # 2D array ([amplitude,spacing], voltage_string)


def max_range(tan_angle):
    L = 37.5 #mm
    distance = []
    for i in range(len(tan_angle)):
        distance.append(L*tan_angle[i]*10**3)
    return distance #nm

def y_error_formula(tan,x,sepfac=12):
    global LASER_WAVELENGTH
    x_err = 2.15*10**(-3) + pixel_mm((float(sepfac/2))) #mm
    tan_error = []
    for i in range(len(x)):
        tan_err = np.sqrt((0.1/LASER_WAVELENGTH)**2 + (x_err/x[i])**2 )*tan[i]
        tan_error.append(tan_err)
        
    return tan_error

def distance_error(tan,x,sepfac=12):
    distance_error = []
    for i in range(len(tan)):
        distance = max_range(tan)
        y_e = y_error_formula(tan,x,sepfac)
        y_error = np.sqrt((y_e[i]/tan[i])**2 + (0.1/37.5)**2)*distance[i]
        distance_error.append(y_error)
    return distance_error

def finish_plot(title='', x_label='', y_label =''):
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(bbox_to_anchor=(1.0, 0.6))
    plt.title(title)
    plt.show()
    

def find_param(xdata,ydata):
    param,cov = np.polyfit(xdata,ydata,1, cov=True)
    return param,cov

def sub_main(input_path,sepfac = 12):
    
    SEPARATION_FACTOR = sepfac
    main_path = input_path
    all_results = analyse_multiple(main_path,SEPARATION_FACTOR)
    y_err = 2.15*10**(-3) + pixel_mm((float(SEPARATION_FACTOR/2)))
    results = [item[0] for item in all_results]
    labels = [item[1] for item in all_results]
    
    volt_data = []
    angle_data = []
    for i in range(len(results)):
        values = analyze(results[i])
        
        plot_error(values[2],values[1],x_error = 0.05, y_error = y_err,label = labels[i])
        
    xdata = np.linspace(0, 100, 1000)
    plt.plot(xdata,0*xdata + 6.9,'r',linestyle='dashed') #New parameter
    finish_plot('spacing between fringes', 'Voltage [V]', 'spacing [mm]')
    
    for i in range(len(results)):
        values = analyze(results[i])
        plot_error(values[2],values[0],label = labels[i])
    finish_plot('amplitude')
    
    for i in range(len(results)):
        values = analyze(results[i])
        volt_data.extend(values[2])
        angle_data.extend(max_range(angle_of_PZT(values[1])))
        plot_error(values[2],max_range(angle_of_PZT(values[1])),x_error = 0.05, y_error = distance_error(angle_of_PZT(values[1]), values[2],SEPARATION_FACTOR),label = labels[i])
    
    """
    lin_param = find_param(volt_data,angle_data)[0]
    xdata = np.linspace(0, 150, 1000)
    
    plt.plot(xdata,2.438*10**(-5)*xdata,'k') #Fixed parameter
    plt.plot(xdata,lin_param[0]*xdata,'k',linestyle='dashed') #New parameter
    
    plt.plot(xdata,0*xdata+ 0.0027,'g',linestyle='dashed') #New parameter
    plt.plot(xdata,0*xdata+ 0.0033,'r',linestyle='dashed') #New parameter
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0, 0.003, "a = "+str(float('%.4g' % lin_param[0])),bbox=props)
    #plt.text(0, 0.0012, "a = "+str(float('%.4g' % lin_param[0])),bbox=props)
    finish_plot("distance PZT moved from initial position",'Voltage [V]',r'maximal distance [mm]')
    """
    lin_param = find_param(volt_data,angle_data)[0]
    xdata = np.linspace(0, 150, 1000)
    
    plt.plot(xdata,2.438*10**(-5)*37.5*1000*xdata,'k') #Fixed parameter
    plt.plot(xdata,lin_param[0]*xdata,'k',linestyle='dashed') #New parameter
    
    plt.plot(xdata,0*xdata+ 0.0026*37.5*1000,'g',linestyle='dashed') #New parameter
    plt.plot(xdata,0*xdata+ 0.0032*37.5*1000,'r',linestyle='dashed') #New parameter
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0, 80, "a = "+str(float('%.4g' % lin_param[0])),bbox=props)
    #plt.text(0, 0.0012, "a = "+str(float('%.4g' % lin_param[0])),bbox=props)
    finish_plot("distance PZT moved from initial position",'Voltage [V]',r'maximal distance [$\mu$m]')
   
    for i in range(len(results)):
        values = analyze(results[i])
        volt_data.extend(values[2])
        angle_data.extend(max_range(angle_of_PZT(values[1])))
        plt.errorbar(values[2], np.zeros(len(values[1])), xerr = 0.05, yerr=distance_error(angle_of_PZT(values[1]), values[2],SEPARATION_FACTOR),fmt='o')
    plt.show()
    #finish_plot("hysteresis",'Voltage [V]',r'tan(${\Theta}$)')

def main():
    #PATH = os.path.dirname(easygui.fileopenbox(filetypes='*.csv'))
    #results = main_analysis(PATH)
    SEPARATION_FACTOR = 12
    main_path = "D:\\TVHolography\\NewFolder\\photos\\Analyse"
    sub_main(main_path,SEPARATION_FACTOR)
    
    SEPARATION_FACTOR = 12
    main_path = "D:\\TVHolography\\NewFolder\\photos\\Analyse - Copy"
    sub_main(main_path,SEPARATION_FACTOR)
    
    SEPARATION_FACTOR = 10
    main_path = "D:\\TVHolography\\New folder"
    sub_main(main_path,SEPARATION_FACTOR)
    
    # hysteresis
    SEPARATION_FACTOR = 10
    main_path = "D:\\TVHolography\\NewFolder\\photos\\hysteresis"
    sub_main(main_path,SEPARATION_FACTOR)
    
    # Green laser
    SEPARATION_FACTOR = 10
    main_path = "D:\\TVHolography\\NewFolder\\photos\\Analyse_Green"
    sub_main(main_path,SEPARATION_FACTOR)
    
    
    
    
    ### list 'results' is 2D list, contains [ [[amp], [separation]], [reference voltage, voltage]]
    ### so it [ values(2D) ,  voltages(1D)]
    #values = analyze(results)
    #result_plot(values)
    


main()
x_err = 2.15*10**(-3) + pixel_mm((float(12/2))) #mm
print(x_err)
tan_err = np.sqrt((0.1/532)**2 + (x_err/6.88)**2 )*0.000077
print(tan_err)